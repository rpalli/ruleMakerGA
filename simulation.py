#import from other parts of ruleMaker
import utils as utils
#import python modules
from random import random, shuffle, randint
import operator
import networkx as nx
import itertools as itertool
import scipy.stats as regress
import numpy as np

class modelClass:
	def __init__(self,graph, sss): 
		#remove self loops from the graph
		for node in graph.nodes():
			repeat=True
			while(repeat):
				repeat=False
				if node in graph.successors(node):
					graph.remove_edge(node,node)
					repeat=True
		#set up empty lists and dicts for later
		individualParse=[] # list of the number of shadow and nodes that contribute to each node, in order by index num
		andNodeList=[] #a list of the shadow nodes that represent and relations between incoming edge
		andNodeInvertList=[] # keeps track of which incoming nodes for each node need to be inverted
		andLenList=[] # keeps track of how many nodes are coming into each shadow AND node
		nodeList=graph.nodes()#define the node list simply as the nodes in the graph. 
		nodeDict={} #identifies names of nodes with their index in the node list- provide name, get index
		for i in range(0,len(nodeList)):
			nodeDict[nodeList[i]]=i #constructs the node dict so we can easily look up nodes
		counter=int(0) #keeps track of where we are in the generic individual
		initValueList=[] #starting states for nodes
		for j in range(0,len(sss)): #construct empty lists to stick values in later for intiial value list
			initValueList.append([])
		
		#find all possible combinations of upstream contributors for each node. These become the shadow And nodes
		for i in range(0,len(nodeList)):
			preds=graph.predecessors(nodeList[i]) # get predecessors of node. 
			while len(preds)>3: #handle case where there are too many predecessors by truncation
				preds.pop(randint(0, len(preds)-1))
			for j in range(0,len(preds)):
				preds[j]=nodeDict[preds[j]]
			# the followign section constructs a list of possible node orders
			# this is accomblished by finding all possible subsets of the list of predecessor nodes
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

			# create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
			activities=[] #list to store activities of nodes (a vs i)
			for sequence in possibilities:
				activity=[]
				for node in sequence:
					if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
						activity.append(False)
					else:
						activity.append(True)
				activities.append(activity)
			andNodeList.append(possibilities)
			andNodeInvertList.append(activities)
			andLenList.append(len(possibilities))
			
			# construct the list of lengths of possibilties for each node, add to the counter that keeps track of how many bits are necessary
			individualParse.append(counter)
			counter=counter+len(possibilities)
		self.size=counter
		individualParse.append(counter)
		self.individualParse=individualParse #index of start value of current node on the individual
		self.andNodeList=andNodeList # shadow and node inputs
		self.andNodeInvertList=andNodeInvertList # keeps track of which incoming nodes for each node need to be inverted
		self.andLenList=andLenList # keeps track of length of above inputOrderList for each node
		self.nodeList=nodeList #define the node list simply as the nodes in the graph. 
		self.nodeDict=nodeDict #identifies names of nodes with their index in the node list.. provide name, get index
		self.initValueList=initValueList #puts an empty and correctly structured initValueList together for later population. 
		
class paramClass:
	def __init__(self,):    
		self.bioReplicates=1
		self.cells=1000
		self.samples=10
		self.generations=5 # generations to run
		self.popSize=8 #size of population
		self.mu= 8 #individuals selected
		self.lambd= 8 #children produced
		self.iters=100 #number of simulations to try in asynchronous mode
		self.genSteps=100 # steps to find steady state with fake data
		self.simSteps=100 # number of steps each individual is run when evaluating
		self.crossoverProb=.3 # prob of3crossing over a particular parent
		self.mutationProb=.7 # prob of mutating a particular parent
		self.rewire=False
		self.async=False # run in asynchronous mode
		self.adaptive=True
		self.verbose=True
		self.bitFlipProb=.1 # prob of flipping bits inside mutation
		self.sigmaNetwork=0
		self.sigmaNode=0
		self.hofSize=10
		self.trials=1
		self.IC=0 #tells the information criterion... 0- no criterion; 1- AIC; 2- BIC
		
def updateBool(currentNode,oldValue,nodeIndividual, model):
	# we update node by updating shadow and nodes then combining them to update or nodes. 
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	#update nodes with more than one input
	# first deal with case of simple logic without need of linear regression
	counter =0
	orset=[]
	# go through list of possible shadow and nodes to see which ones actually contribute
	for andindex in range(len(nodeIndividual)):
		if nodeIndividual[andindex]==1:
			# if a shadow and contributes, compute its value using its upstream nodes
			# calculate value of first then use and to append rest in list of predecessors
			newval=oldValue[andNodes[andindex][0]]^andNodeInvertList[andindex][0]
			for addnode in range(1,len(andNodes[andindex])):
				newval=newval&(oldValue[andNodes[andindex][addnode]]^andNodeInvertList[andindex][addnode])
			orset.append(newval)
	#combine the shadow and nodes with or operations
	newval=orset.pop()
	for val in orset:
		newval= newval | val
	return newval

def updateFuzzy(currentNode,oldValue,nodeIndividual, model):
	# we update node by updating shadow and nodes then combining them to update or nodes. 
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	#update nodes with more than one input
	# first deal with case of simple logic without need of linear regression
	counter = 0
	orset=[]
	# go through list of possible shadow and nodes to see which ones actually contribute
	for andindex in range(len(nodeIndividual)):
		if nodeIndividual[andindex]==1:
			# if a shadow and contributes, compute its value using its upstream nodes
			# calculate value of first then use and to append rest in list of predecessors
			if andNodeInvertList[andindex][0]:
				newval=1-oldValue[andNodes[andindex][0]]
			else:
				newval=oldValue[andNodes[andindex][0]]
			for addnode in range(1,len(andNodes[andindex])):
				if oldValue[andNodes[andindex][addnode]]:
					newval=min(newval,1-andNodeInvertList[andindex][addnode])
				else:
					newval=min(newval,andNodeInvertList[andindex][addnode])
			orset.append(newval)
	#combine the shadow and nodes with or operations
	newval=orset.pop()
	for val in orset:
		newval= max(newval , val)
	return newval

#run a simulation given a starting state
def runBool(individual, model,  simSteps, initValues, params, knockouts, knockins):
	if params.async:
		return asyncBool(individual, model, simSteps, initValues, params.iters, knockouts, knockins)
	else:
		return syncBool(individual, model, simSteps, initValues, knockouts, knockins)

def runFuzzy(individual, model,  simSteps, initValues, params, knockouts, knockins):
	if params.async:
		return asyncFuzzy(individual, model,  initValues, params.iters, knockouts, knockins)
	else:
		return syncFuzzy(individual, model, initValues, knockouts, knockins)

def syncBool(individual, model, simSteps, initValues, knockouts, knockins):
	# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
	# set up data storage for simulation, add step 0
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))

	#iterate over number of steps necessary
	for step in range(0,simSteps):
		oldValue=list(newValue)
		for i in range(0,len(model.nodeList)):
			#find start and finish for each node to update from the individualParse list
			if i in knockouts:
				temp=0
			elif i in knockins:
				temp=1
			elif model.andLenList[i]==1:
				temp=   oldValue[i]^model.andNodeInvertList[i][0][0]
			elif model.andLenList[i]==0:
				temp=oldValue[i]
			else:
				if i==len(model.nodeList)-1:
					end= model.size
				else:
					end=model.individualParse[i+1]	 
				temp=updateBool(i,oldValue,individual[model.individualParse[i]:end], model)
			newValue[i]=temp
		simData.append(list(newValue))
	avgs= [[] for x in range(0,len(newValue))]
	for element in range(0,len(avgs)):
		avgs[element]=[simData[step][element] for step in range(len(simData)-10,len(simData)) ]
	avg= [1.0*np.sum(element)/len(element) for element in avgs]
	return avg

def syncFuzzy(individual, model, simSteps, initValues, knockouts, knockins):
	# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
	# set up data storage for simulation, add step 0
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	#iterate over number of steps necessary
	for step in range(0,simSteps):
		oldValue=list(newValue)
		for i in range(0,len(model.nodeList)):
			#find start and finish for each node to update from the individualParse list
			if i in knockouts:
				temp=0
			elif i in knockins:
				temp=1
			elif model.andLenList[i]==1:
				temp=   oldValue[i]^model.andNodeInvertList[i][0][0]
			elif model.andLenList[i]==0:
				temp=oldValue[i]
			else:
				if i==len(model.nodeList)-1:
					end= model.size
				else:
					end=model.individualParse[i+1]	 
				temp=updateFuzzy(i,oldValue,individual[model.individualParse[i]:end], model)
			newValue[i]=temp
		simData.append(list(newValue))
	avgs= [[] for x in range(0,len(newValue))]
	for element in range(0,len(avgs)):
		avgs[element]=[simData[step][element] for step in range(len(simData)-10,len(simData)) ]
	avg= [1.*np.sum(element)/len(element) for element in avgs]
	return avgs

#run asyncrhnonoyus simulation and average it over iters trials
def asyncBool(individual, model, simSteps, initValues, iters, knockouts, knockins):
	sum1=[0 for x in range(0,len(initValues))]
	# run iterations with different orderings
	for i in range(0,params.iters):
		# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
		# set up data storage for simulation, add step 0
		newValue=list(initValues)
		# totalNodes=0... turn this on for information criterion

		# set up the sequence of nodes to be updated
		seq=range(0,len(model.nodeList))
		#iterate over number of steps necessary
		for step in range(0,simSteps):
			oldValue=list(newValue)
			#shuffle- async 
			shuffle(seq)
			for seq[i] in range(0,len(model.nodeList)):
				#find start and finish for each node to update from the individualParse list
				if seq[i] in knockouts:
					temp=0
				elif seq[i] in knockins:
					temp=1
				elif model.andLenList[seq[i]]==1:
					if model.andNodeInvertList[seq[i]][0][0]:
						temp=   1- oldValue[seq[i]]
					else:
						temp=  oldValue[seq[i]]
				elif model.andLenList[seq[i]]==0:
					temp=oldValue[seq[i]]
				else:
					if seq[i]==len(model.nodeList)-1:
						end= model.size
					else:
						end=model.individualParse[i+1]	 
					temp=updateBool(seq[i],newValue,individual[model.individualParse[seq[i]]:end],  model)
				newValue[seq[i]]=temp
		for j in range(0,len(sum1)):
			sum1[j]=sum1[j]+newValue[j]
	avgs=list(sum1)
	for i in range(0,len(sum1)):
		avgs[i]=sum1[i]/float(params.iters)
	return avgs

#run asyncrhnonoyus simulation and average it over iters trials
def asyncFuzzy(individual, model, simSteps, initValues, iters, knockouts, knockins):
	sum1=[0 for x in range(0,len(initValues))]
	# run iterations with different orderings
	for i in range(0,params.iters):
		# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
		# set up data storage for simulation, add step 0
		newValue=list(initValues)
		# totalNodes=0... turn this on for information criterion

		# set up the sequence of nodes to be updated
		seq=range(0,len(model.nodeList))
		#iterate over number of steps necessary
		for step in range(0,simSteps):
			oldValue=list(newValue)
			#shuffle- async 
			shuffle(seq)
			for i in range(0,len(model.nodeList)):
				#find start and finish for each node to update from the individualParse list
				if seq[i] in knockouts:
					temp=0
				elif seq[i] in knockins:
					temp=1
				elif model.andLenList[seq[i]]==1:
					if model.andNodeInvertList[seq[i]][0][0]:
						temp=   1- oldValue[seq[i]]
					else:
						temp=  oldValue[seq[i]]
				elif model.andLenList[seq[i]]==0:
					temp=oldValue[seq[i]]
				else:
					if seq[i]==len(model.nodeList)-1:
						end= model.size
					else:
						end=model.individualParse[seq[i]+1]	 
					temp=updateFuzzy(seq[i],newValue,individual[model.individualParse[seq[i]]:end],  model)
				newValue[seq[i]]=temp
		for j in range(0,len(sum1)):
			sum1[j]=sum1[j]+newValue[j]
	avgs=list(sum1)
	for i in range(0,len(sum1)):
		avgs[i]=sum1[i]/float(params.iters)
	return avgs

# init value generator for EBNs
def genEBNInitValues(individual, model,sampleProbs):
	#return [True if (random()<sampleProbs[node]) else False for node in range(0,len(sampleProbs))]
	initValues=[False for x in range(0,len(model.nodeList))]
	for node in range(0,len(sampleProbs)):
		if random()<sampleProbs[node]:
			initValues[node]=True
	return initValues

# main EBN simulation code. Runs an EBN
def EBNbool(individual, model, cells, initProbs, params, KOs, KIs, iteratorDict):
	cellArray=[]
	simSteps=3*len(model.nodeList)
	if not str(individual) in iteratorDict:
		iteratorDict[str(individual)]={}
	valueDict=iteratorDict[str(individual)]

	for j in range(0,cells):
		# shuffle nodes to be initially called.... 
		#simulations that are truly random starting states should have all nodes added to this list
		#get initial values for all nodes
		initValues=genEBNInitValues(individual, model,initProbs)
		if str(initValues) in valueDict.keys():
			vals=valueDict[str(initValues)]
		else:
			# run Boolean simulation with initial values and append
			vals=runBool(individual, model,simSteps, initValues, params, KOs, KIs)
			valueDict[str(initValues)]=vals
		cellArray.append(vals)
	return [1.*np.sum(col) / float(cells) for col in zip(*cellArray)]


def EBNfuzzy(individual, model, cells, initProbs, params, KOs, KIs, iteratorDict):
	cellArray=[]
	simSteps=3*len(model.nodeList)
	if not str(individual) in iteratorDict:
		iteratorDict[str(individual)]={}
	valueDict=iteratorDict[str(individual)]

	for j in range(0,cells):
		# shuffle nodes to be initially called.... 
		#simulations that are truly random starting states should have all nodes added to this list
		#get initial values for all nodes
		initValues=genEBNInitFuzzies(individual, model,initProbs)
		if str(initValues) in valueDict.keys():
			vals=valueDict[str(initValues)]
		else:
			# run Boolean simulation with initial values and append
			vals=runFuzzy(individual, model,simSteps, initValues, params, KOs, KIs)
			valueDict[initValues]=vals
		cellArray.append(vals)
	return [1.*np.sum(col) / float(cells) for col in zip(*cellArray)]