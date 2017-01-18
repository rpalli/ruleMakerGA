import utils as utils
import liu_networks as liu
#import stuff
from random import random
from random import shuffle
from random import gauss
import numpy as numpy
from math import floor, ceil, log
import operator
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
import networkx as nx
import csv
import itertools as itertool
import copy as copy
import gc as gc
import re as re

class modelClass:
     def __init__(self,graph, sss):         
		repeat=True
		while(repeat):
			repeat=False
			for node in graph.nodes():
				if node in graph.successors(node):
					graph.remove_edge(node,node)
					repeat=True
		
		evaluateNodes=[] #list of nodes that need to be compared to the steady state values (sss)
		individualParse=[] #list of triples that contain the start and end of the noder order value and operator flag list for each node
		inputOrderList=[] #list of lists of possible node orders for each node
		inputOrderInvertList=[] # keeps track of which incoming nodes for each node need to be inverted
		possibilityNumList=[] # keeps track of length of above inputOrderList for each node
		earlyEvalNodes=[] # nodes that don't have initial value and need to be re-evaluated early on in the simulation
		
		nodeList=graph.nodes()#define the node list simply as the nodes in the graph. 
		nodeDict={} #identifies names of nodes with their index in the node list
		for i in range(0,len(nodeList)):
			nodeDict[nodeList[i]]=i #constructs the node dict so we can easily look up nodes
		
		counter=int(0) #keeps track of where we are in the generic individual
		initValueList=[] #starting states for nodes
		for j in range(0,len(sss)): #construct empty lists to stick values in later for intiial value list
				initValueList.append([])
		

		for i in range(0,len(nodeList)):
			preds=graph.predecessors(nodeList[i]) # get predecessors of node. 
			if len(preds)>15: #handle case where there are too many predecessors by truncation
				preds=preds[1:15]
			for j in range(0,len(preds)):
				preds[j]=nodeDict[preds[j]]
			
			# the followign section constructs a list of possible node orders
			withNones = zip(preds, itertool.repeat('empty'))
			possibilities=list(itertool.product(*withNones))
			for j in range(0,len(possibilities)):
				possibilities[j]=list(possibilities[j])
				while 'empty' in possibilities[j]:
					possibilities[j].remove('empty')
				while [] in possibilities[j]:
					possibilities[j].remove([])
			while [] in possibilities:
				possibilities.remove([])
			
			#store activities of nodes in correspondence with node order. 
			activities=[] #list to store activities of nodes (a vs i)
			for sequence in possibilities:
				activity=[]
				for node in sequence:
					if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
						activity.append(False)
					else:
						activity.append(True)
				activities.append(activity)
			inputOrderList.append(possibilities)
			inputOrderInvertList.append(activities)
			possibilityNumList.append(len(possibilities))
			#deal with nodes with no incoming or only one incoming first
			if len(possibilities)==0:
				individualParse.append([counter,counter,counter])
			elif len(possibilities)==1:
				individualParse.append([counter,counter,counter])
				counter=int(counter+1)
			else:
				logNum=ceil(log(len(possibilities))/log(2))#determine number of bits necessary to store information of which node order is the one for that individual
				individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(preds)-1)]) #set up limits: lognum-1 for list then len(preds)-1 for the activity flags
				counter=int(counter+logNum+len(preds)-1) #uptick the counter so we know where to start on the next node
			#set up lists of initial values and values that need to be evaluated early
			for j in range(0,len(sss)):
				ss=sss[j]
				if  nodeList[i] in sss[0].keys():
					initValueList[j].append(ss[nodeList[i]])
					evaluateNodes.append(i)
				else:
					initValueList[j].append(0.5)
					earlyEvalNodes.append(i)
		self.size=counter
		self.evaluateNodes=evaluateNodes #list of nodes that need to be compared to the steady state values (sss)
		self.individualParse=individualParse #list of triples that contain the start and end of the noder order value and operator flag list for each node
		self.inputOrderList=inputOrderList #list of lists of possible node orders for each node
		self.inputOrderInvertList=inputOrderInvertList # keeps track of which incoming nodes for each node need to be inverted
		self.possibilityNumList=possibilityNumList # keeps track of length of above inputOrderList for each node
		self.earlyEvalNodes=earlyEvalNodes
		self.nodeList=nodeList#define the node list simply as the nodes in the graph. 
		self.nodeDict=nodeDict
		self.initValueList=initValueList

class paramClass:
	def __init__(self):     
		self.simSteps=100 # number of steps each individual is run when evaluating
		self.generations=15 # generations to run
		self.popSize=100 #size of population
		self.mu= 100#individuals selected
		self.lambd= 200#children produced
		self.bitFlipProb=.1 # prob of flipping bits inside mutation
		self.crossoverProb=.2 # prob of crossing over a particular parent
		self.mutationProb=.2 # prob of mutating a particular parent
		self.async=False # run in asynchronous mode
		self.iters=100 #number of simulations to try in asynchronous mode
		self.complexPenalty=False #penalize models with increased complexity
		self.genSteps=10000 # steps to find steady state with fake data
		self.sigmaNetwork=0
		self.sigmaNode=0
		self.hofSize=10
		self.cells=100

def fuzzyAnd(num1, num2, index1, index2,corrMat):
	return min(num1,num2)
def fuzzyOr(num1, num2, index1, index2,corrMat):
	return max(num1, num2)
def naiveAnd(num1, num2, index1, index2,corrMat):
	return num1*num2	
def naiveOr(num1, num2, index1, index2,corrMat):
	return (num1+num2-(num1*num2))
def propOrE1(num1, num2, index1, index2,corrMat):
	if num1==1 or num2==1:
		return 1
	else:
		return max(0,min(1,(num1+num2-(num1*corrMat[index2][index1]+num2*corrMat[index1][index2])/2)))
def propAndE1(num1,num2,index1,index2,corrMat):
	return max(0,min(1,((num1*corrMat[index2][index1])/2+(num2*corrMat[index1][index2])/2)))

def propOrE2(A,B,AgivenB, BgivenA):
	return max(0,min(1,A+B-(B*AgivenB+A*BgivenA)/2))

def propAndE2(A,B,AgivenB, BgivenA):
	return max(0,min(1,(A*BgivenA+B*AgivenB)/2))

def Inv(x, inverter): #inverts if necessary then applies hill fun
	if inverter:
		return (1-x)
	else:
		return (x)

class simulatorClass:
	def __init__(self,simTyping):
		self.simType=simTyping
		if simTyping=='propE1':
			self.And=propAndE1
			self.Or=propOrE1
			self.corrMat={}
		if simTyping=='propE1=2':
			self.And=propAndE2
			self.Or=propOrE2
			self.corrMat={}
		if simTyping=='fuzzy':
			self.And=fuzzyAnd
			self.Or=fuzzyOr
			self.corrMat=0
		if simTyping=='propNaive':
			self.And=naiveAnd
			self.Or=naiveOr
			self.corrMat=0

#run a simulation given a starting state then add noise if noise is TRUE
def runModel(individual, params, model, simulator, initValues, noise):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out.  
	#list of triples which gives the total list of nodes. 
	counter=0;
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	seq=range(0,len(model.individualParse))
	for node in range(0,len(model.earlyEvalNodes)):
		newValue[model.earlyEvalNodes[node]]=updateNet(model.earlyEvalNodes[node],newValue,individual,individualParse[seq[i]], model,simulator)
	for step in range(0,params.simSteps):
		oldValue=list(newValue)
		if params.async:
			shuffle(seq)
		for i in range(0,len(model.individualParse)):	
			if params.async:
				temp=updateNet(seq[i],newValue,individual, model.individualParse[seq[i]],  model,simulator)
			else:
				temp=updateNet(seq[i],oldValue,individual, model.individualParse[seq[i]], model,simulator)

			if not noise:
				newValue[seq[i]]=temp
			else:
				newValue[seq[i]]=temp+params.sigmaNode*(1-temp/(.1+temp)-temp/10)*random()
		simData.append(list(newValue))
	avg= [0 for x in range(0,len(newValue))]
	stdev= [0 for x in range(0,len(newValue))]
	for step in range(0,len(simData)):
		for element in range(0,len(avg)):
			simData[step][element]=simData[step][element]+params.sigmaNetwork*random()
	for step in range(len(simData)-10,len(simData)):
		for element in range(0,len(avg)):
			avg[element]=avg[element]+simData[step][element]
	for element in range(0,len(avg)):
		avg[element]=avg[element]/10
	#for step in range(len(simData)-10,len(simData)):
		#for element in range(0,len(avg)):
			#stdev[element]=stdev[element]+(simData[step][element]-avg[element])^2
	#return (avg, stdev)
	return avg

#run a simulation and average it over iters trials
def averageResultModelSim(individual, params, model, simulator, initValues, iters):
	sum=[0 for x in range(0,len(initValues))]
	for i in range(0,iters):
		avg=runModel( individual, params, model, simulator, initValues,False)
		for j in range(0,len(sum)):
			sum[j]=sum[j]+avg[j]
	avgs=list(sum)
	for i in range(0,len(sum)):
		avgs[i]=sum[i]/float(iters)
	return avgs

def updateNet(currentNode,oldValue,individual, triple, model,simulator):
	inputOrder=model.inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
	inputOrderInvert=model.inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if model.possibilityNumList[currentNode]>0:
		logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
		inputOrder=list(inputOrder[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]]) # determine what order of inputs is
		inputOrderInvert=inputOrderInvert[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
		if len(inputOrder)==0:
			value=oldValue[currentNode] #if no inputs, maintain value
			return value
		elif len(inputOrder)==1:
			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
			value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
			return value
		else:
			#update nodes with more than one input
			# update and then or
			upstreamVals=[]
			for upstream in range(0,len(inputOrder)):
				upstreamVals.append(Inv(oldValue[inputOrder[upstream]],inputOrderInvert[upstream]))
			counter =0
			# print(upstreamVals)
			# print(logicOperatorFlags)
			while counter < len(logicOperatorFlags) and counter+1<len(inputOrder):
				if logicOperatorFlags[counter]==0:
					tempVal=simulator.And(upstreamVals[counter],upstreamVals[counter+1],inputOrder[counter], inputOrder[counter+1],simulator.corrMat)
					inputOrder.pop(counter)
					inputOrder.pop(counter)
					logicOperatorFlags.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.insert(counter,tempVal)
				else:
					counter=counter+1
				# print(upstreamVals)

			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
			while len(upstreamVals)>1:
				tempVal=simulator.Or(upstreamVals.pop(0),upstreamVals.pop(0),inputOrder.pop(0),inputOrder.pop(0),simulator.corrMat)
				upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
			return upstreamVals[0]
	else:
		#returns savme value if now inputs
		return oldValue[currentNode]						
def updateNetOLD(currentNode,oldValue,individual, triple, model,simulator):
	inputOrder=model.inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
	inputOrderInvert=model.inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if model.possibilityNumList[currentNode]>0:
		logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
		inputOrder=list(inputOrder[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]]) # determine what order of inputs is
		inputOrderInvert=inputOrderInvert[utils.bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
		if len(inputOrder)==0:
			value=oldValue[currentNode] #if no inputs, maintain value
			return value
		elif len(inputOrder)==1:
			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
			value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
			return value
		else:
			#update nodes with more than one input
			# update and then or
			upstreamVals=[]
			for upstream in range(0,len(inputOrder)):
				upstreamVals.append(Inv(oldValue[inputOrder[upstream]],inputOrderInvert[upstream]))
			counter =0
			# print(upstreamVals)
			# print(logicOperatorFlags)
			while counter < len(logicOperatorFlags) and counter+1<len(inputOrder):
				if logicOperatorFlags[counter]==0:
					tempVal=simulator.And(upstreamVals[counter],upstreamVals[counter+1],inputOrder[counter], inputOrder[counter+1],simulator.corrMat)
					inputOrder.pop(counter)
					inputOrder.pop(counter)
					logicOperatorFlags.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.insert(counter,tempVal)
				else:
					counter=counter+1
				# print(upstreamVals)

			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
			while len(upstreamVals)>1:
				tempVal=simulator.Or(upstreamVals.pop(0),upstreamVals.pop(0),inputOrder.pop(0),inputOrder.pop(0),simulator.corrMat)
				upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
			return upstreamVals[0]
	else:
		#returns savme value if now inputs
		return oldValue[currentNode]						

def updateNode(currentNode,oldValue,individual, model,simulator):
	
	oldTriple=model.individualParse[currentNode]
	triple=[]
	triple.append(0)
	triple.append(oldTriple[1]-oldTriple[0])
	triple.append(oldTriple[2]-oldTriple[0])
	retinue = updateNet(currentNode,oldValue,individual, triple,model,simulator)
	return retinue
					

