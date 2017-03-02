#import from other parts of ruleMaker
import utils as utils
#import python modules
from random import random
from random import shuffle
import operator
import networkx as nx
import itertools as itertool
import scipy.stats as regress

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
		evaluateNodes=[] #list of nodes that need to be compared to the steady state values (sss)
		individualParse=[] # list of the number of shadow and nodes that contribute to each node, in order by index num
		earlyEvalNodes=[] # nodes that don't have initial value and need to be re-evaluated early on in the simulation.. 
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
			if len(preds)>15: #handle case where there are too many predecessors by truncation
				preds=preds[1:15]
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
			if  nodeList[i] in sss[0].keys():
				evaluateNodes.append(i)
			else:
				earlyEvalNodes.append(i)
		self.size=counter
		individualParse.append(counter)
		self.evaluateNodes=evaluateNodes #list of nodes that need to be compared to the steady state values (sss)
		self.individualParse=individualParse #index of start value of current node on the individual
		self.andNodeList=andNodeList # shadow add node inputs
		self.andNodeInvertList=andNodeInvertList # keeps track of which incoming nodes for each node need to be inverted
		self.andLenList=andLenList # keeps track of length of above inputOrderList for each node
		self.earlyEvalNodes=earlyEvalNodes
		self.nodeList=nodeList #define the node list simply as the nodes in the graph. 
		self.nodeDict=nodeDict #identifies names of nodes with their index in the node list.. provide name, get index
		self.initValueList=initValueList #puts an empty and correctly structured initValueList together for later population. 

		
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
		self.genSteps=100 # steps to find steady state with fake data
		self.sigmaNetwork=0
		self.sigmaNode=0
		self.hofSize=10
		self.cells=1000
		self.samples=4
		self.trials=10
		self.IC=3 #tells the information criterion... 0- no criterion; 1- AIC; 2- BIC

def fuzzyAnd(num1, num2):
	return min(num1,num2)
def fuzzyOr(num1, num2):
	return max(num1, num2)
def naiveAnd(num1, num2):
	return num1*num2	
def naiveOr(num1, num2):
	return (num1+num2-(num1*num2))
def propOr(A,B,AgivenB, BgivenA):
	return max(0,min(1,A+B-((B*AgivenB+A*BgivenA)/2)))
def propAnd(A,B,AgivenB, BgivenA):
	return max(0,min(1,(A*BgivenA+B*AgivenB)/2))
def Inv(x, inverter): #inverts if necessary then applies hill fun
	if inverter:
		return (1-x)
	else:
		return (x)

class simulatorClass:
	#class for simulations.... the simType gives the formalism to be used
	# if statements below map 'and' and 'or' to the correct functions for the simtype
	def __init__(self,simTyping):
		self.simType=simTyping
		self.trainingData=[]
		if simTyping=='prop':
			self.And=propAnd
			self.Or=propOr
			self.corrMat={}
			self.switch=1			
			self.simSteps=1
		if simTyping=='fuzzy':
			self.And=fuzzyAnd
			self.Or=fuzzyOr
			self.corrMat=0
			self.switch=0
			params=paramClass()
			self.simSteps= params.simSteps
		if simTyping=='propNaive':
			self.And=naiveAnd
			self.Or=naiveOr
			self.corrMat=0
			self.switch=0
			params=paramClass()
			self.simSteps= params.simSteps
	def train(self, model):
		#self.trainingData
		self.andTrainer=[]
		self.orDicts=[]
		# saved as [m,b,m,b],....[m, b, m, b], [train]
		for currentNode in range(0,len(model.nodeList)):
			newdict={}
			self.orDicts.append(newdict) # add a blank dictionary for the or regressions on each node...
			andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
			andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
			if model.andLenList[currentNode]==0:
				self.andTrainer.append([]) #if no inputs, no and nodes
			elif len(andNodes)==1 : 
				#if only one input, still don't need shadow and nodes
				self.andTrainer.append([])
			else:
				#if more than one input, we need to train and nodes
				tempTrainer=[] # store all and nodes for a particular node

				#go through and nodes and produce solutions for each
				for andindex in range(0,len(andNodes)):
					temp=[] #store values for a particular and node
					if not len(andNodes[andindex])==1:
						#do iterative calculation on AND nodes storing as you go.... 
						for nodeIn in range(len(andNodes[andindex])):
							#set up initial vars with first node in the set of nodes to be connnected by AND
								if nodeIn==0:
									node1index=andNodes[andindex][0]
									node1=[]
									for sample in self.trainingData:
										node1.append(sample[node1index])
									tempdata=node1
								else:
									# update training data and value using AND rule
									node1index=andNodes[andindex][nodeIn]
									node1=[]
									for sample in self.trainingData:
										node1.append(sample[node1index])
									slope1, intercept1, r_value, p_value, std_err = regress.linregress(tempdata,node1)
									slope2, intercept2, r_value, p_value, std_err = regress.linregress(node1,tempdata)
									#store relevant linreg results
									temp.append([slope1,intercept1, slope2, intercept2])
									for sample in range(len(self.trainingData)):
										tempdata[sample]=self.And(self.trainingData[sample][node1index],tempdata[sample],slope1*tempdata[sample]+intercept1,slope2*self.trainingData[sample][node1index]+intercept2)
						# append the tempdata to do or calculations....
						temp.append(tempdata)
					tempTrainer.append(temp)
				self.andTrainer.append(tempTrainer)

def updateNode(currentNode,oldValue,nodeIndividual, model,simulator):
	# we update node by updating shadow and nodes then combining them to update or nodes. 
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if model.andLenList[currentNode]==0:
		return oldValue[currentNode] #if no inputs, maintain value
	elif len(andNodes)==1 or sum(nodeIndividual)==0: 
		#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
		if nodeIndividual[0]==1:
			#if only one input, then set to that number
			value=Inv(oldValue[andNodes[0][0]],andNodeInvertList[0][0])
		else:
			value=oldValue[currentNode] #if no inputs, maintain value
		return value
	else:
		#update nodes with more than one input

		# first deal with case of simple logic without need of linear regression
		if simulator.switch==0:
			counter =0
			orset=[]
			# go through list of possible shadow and nodes to see which ones actually contribute
			for andindex in range(len(nodeIndividual)):
				if nodeIndividual[andindex]==1:
					# if a shadow and contributes, compute its value using its upstream nodes
					# calculate value of first then use and to append rest in list of predecessors
					newval=Inv(oldValue[andNodes[andindex][0]],andNodeInvertList[andindex][0])
					for addnode in range(1,len(andNodes[andindex])):
						newval=simulator.And(newval,Inv(andNodes[andindex][addnode],andNodeInvertList[andindex][addnode]))
					orset.append(newval)
			#combine the shadow and nodes with or operations
			newval=orset.pop()
			for val in orset:
				newval=simulator.Or(newval,val)
			return newval
		else:
			#now we have the case of a proportional boolean network to solve. 
			ortraininglist=[]
			orset=[]
			incNodes=[] # and nodes included
			orDict=simulator.orDicts[currentNode]
			for andindex in range(len(nodeIndividual)):
				if nodeIndividual[andindex]==1:
					incNodes.append(andindex)
					# if a shadow and contributes, compute its value using its upstream nodes
					# first for special case of only a single node as pred to shadow and node
					if len(andNodes[andindex])==1:
						node1index=andNodes[andindex][0]
						# add value of that pred to list of shadow and node values
						orset.append(oldValue[andNodes[andindex][0]])
						node1=[]
						# compute value of this shadow and across the training data
						for sample in simulator.trainingData:
							node1.append(sample[node1index])
						ortraininglist.append(node1)
					else:
						#calculate values across AND nodes
						for nodeIn in range(len(andNodes[andindex])):
							#set up initial vars with first node in the set of nodes to be connnected by AND
							
							if nodeIn==0:
								node1index=andNodes[andindex][0]
								node1=[]
								for sample in simulator.trainingData:
									node1.append(sample[node1index])
								temptraining=node1
								tempVal=oldValue[node1index]
							else:
								# update training data and value using AND rule
								node1index=andNodes[andindex][nodeIn]
								slope1=simulator.andTrainer[currentNode][andindex][nodeIn][0]
								intercept1=simulator.andTrainer[currentNode][andindex][nodeIn][1]
								slope2=simulator.andTrainer[currentNode][andindex][nodeIn][2]
								intercept2=simulator.andTrainer[currentNode][andindex][nodeIn][3]
								tempVal=simulator.And(oldValue[node1index],tempVal,slope1*tempVal+intercept1,slope2*oldValue[node1index]+intercept2)
							# append calculated values of training and test 
							ortraininglist.append(simulator.andTrainer[currentNode][andindex][-1])
							orset.append(tempVal)
			# combine or nodes below
			currentrain=[]
			for i in range(len(orset)):
				if i==0:
					#set up initial or value
					currentval=orset[0]
					if len(orset)>1:
						currentrain=ortraininglist[0]
				else:
					# update with OR values
					if tuple(incNodes[0:i]) in orDict.keys():
						slope1=orDict[tuple(incNodes[0:i])][0]
						intercept1=orDict[tuple(incNodes[0:i])][1]
						slope2=orDict[tuple(incNodes[0:i])][2]
						intercept2=orDict[tuple(incNodes[0:i])][3]
						currentrain=orDict[tuple(incNodes[0:i])][4]
					else:
						slope1, intercept1, r_value, p_value, std_err = regress.linregress(currentrain,ortraininglist[i])
						slope2, intercept2, r_value, p_value, std_err = regress.linregress(ortraininglist[i],currentrain)
						for j in range(len(ortraininglist[i])):
							currentrain[j]=simulator.Or(currentrain[j],ortraininglist[i][j],slope1*ortraininglist[i][j]+intercept1,slope2*currentrain[j]+intercept2)
						orDict[tuple(incNodes[0:i])]=[slope1,intercept1,slope2, intercept2,currentrain]
					currentval=simulator.Or(currentval,orset[i],slope1*orset[i]+intercept1,slope2*currentval+intercept2)
			return currentval

#run a simulation given a starting state
def runModel(individual, model, simulator, initValues):
	# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
	params=paramClass()
	# set up data storage for simulation, add step 0
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	totalNodes=0

	# set up the sequence of nodes to be updated
	seq=range(0,len(model.nodeList))
	for node in range(0,len(model.earlyEvalNodes)):
		newValue[model.earlyEvalNodes[node]]=updateNode(model.earlyEvalNodes[node],newValue,individual,individualParse[seq[i]], model,simulator)
	
	#iterate over number of steps necessary
	for step in range(0,simulator.simSteps):
		oldValue=list(newValue)
		if params.async: #shuffle if async
			shuffle(seq)
		for i in range(0,len(model.nodeList)):
			#find start and finish for each node to update from the individualParse list
			if seq[i]==len(model.nodeList)-1:
				end= model.size
			else:
				end=model.individualParse[seq[i]+1]	 
			#sum up number of nodes in rule
			if i==0:
				for bit in range(model.individualParse[seq[i]],end):
					if individual[bit]==1:
						totalNodes=totalNodes+len(model.andNodeInvertList[seq[i]][bit-model.individualParse[seq[i]]])
			#update based on sync or async assumptions
			if params.async:
				temp=updateNode(seq[i],newValue,individual[model.individualParse[seq[i]]:end],  model,simulator)
			else:
				temp=updateNode(seq[i],oldValue,individual[model.individualParse[seq[i]]:end], model,simulator)
			newValue[seq[i]]=temp

		simData.append(list(newValue))
	if simulator.switch==0:
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
		return avg, totalNodes
	else:
		return simData.pop(), totalNodes

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