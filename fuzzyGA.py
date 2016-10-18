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
from scipy.stats import logistic
import re
import urllib2 
import fuzzyNetworkConstructor as constructor
import csv
import itertools as itertool
import copy as copy
import cPickle as pickle

def returnSimParams():
	simSteps=100
	generations=100
	popSize=100
	p=.5
	h= 3
	hillOn=False
	bitFlipProb=.1
	crossoverProb=.2
	mutationProb=.2
	return simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb

def setupGAparams(graph, ss):
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	
	evaluateNodes=[] 
	individualParse=[] 
	nodeOrderList=[] 
	nodeOrderInvertList=[]
	possibilityNumList=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 100
	counter=0
	initValues=[]
	for i in range(0,len(nodeList)):
		preds=graph.predecessors(nodeList[i])
		if len(preds)>15:
			preds=preds[1:15]
		for j in range(0,len(preds)):
			preds[j]=nodeDict[preds[j]]
		withNones = zip(preds, itertool.repeat('empty'))
		
		possibilities=list(itertool.product(*withNones))
		activities=[]
		for j in range(0,len(possibilities)):
			possibilities[j]=list(possibilities[j])
			while 'empty' in possibilities[j]:
				possibilities[j].remove('empty')
			while [] in possibilities[j]:
				possibilities[j].remove([])
		while [] in possibilities:
			possibilities.remove([])
		for sequence in possibilities:
			activity=[]
			for node in sequence:
			#print(nodeList[node[0]])
			#print(nodeList[i])
			#print(graph.edges())
				
				if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
					activity.append(False)
					#print('a')
				else:
					activity.append(True)
					#print('i')
			activities.append(activity)
		nodeOrderList.append(possibilities)
		nodeOrderInvertList.append(activities)
		possibilityNumList.append(len(possibilities))
		if len(possibilities)==0:
			logNum=0
		else:
			logNum=ceil(log(len(possibilities))/log(2))
		individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(possibilities)-1)])
		counter=counter+logNum+len(possibilities)
		if  nodeList[i] in ss.keys():
			initValues.append(ss[nodeList[i]])
			evaluateNodes.append(i)
		else:
			initValues.append(0.5)
	return counter, individualParse, nodeList, nodeOrderList, initValues, evaluateNodes, possibilityNumList, nodeOrderInvertList
	
	
def genRandBits(individualLength):
	arr = numpy.random.randint(2, size=(individualLength,))
	return list(arr)
	
def loadFatimaData(filename,tempDict):
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]!='#':
			kline=line.split(' ')
			tempDict[str.upper(kline[0][3:])]=logistic.cdf(float(kline[1]))

			
def loadFpkms(filename):
	with open(filename) as csvfile:
		data=[]
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			data.append(row)
		return data

def sortFpkms(data,ss):
	for i in range(1,len(data)):
		maxdata=0
		for j in range(1,len(data[i])):
			if float(data[i][j])>maxdata:
				maxdata=float(data[i][j])
		if maxdata==0:
			maxdata=1
		ss[str.upper(data[i][0])]=float(data[i][1])/maxdata
	
def hill(x,h,p,hillOn):
	if hillOn:
		return ((1+h**p)*x**p)/(h**p+x**p)		
	else:
		return x
		
def hillInv(x, inverter,h,p,hillOn):
	if inverter:
		return hill(1-x, h,p,hillOn)
	else:
		return hill(x, h,p,hillOn)
		
def fuzzyAnd(num1, num2):
	return min(num1,num2)
			
def fuzzyOr(num1, num2):
	return max(num1, num2)

def bit2int(bitlist):
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out

def fuzzyUpdate(currentNode,oldValue,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn):
	triple=individualParse[currentNode]
	nodeOrder=nodeOrderList[currentNode]
	nodeOrderInvert=nodeOrderInvertList[currentNode]
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=individual[triple[1]:triple[2]]
		nodeOrder=nodeOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		nodeOrderInvert=nodeOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		if len(nodeOrder)==0:
			value=oldValue[currentNode]
		elif len(nodeOrder)==1:
			#print(nodeOrder[0])
			value=hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0],h,p,hillOn)
		else:
			counter =0
			if logicOperatorFlags[0]==0:
				value=fuzzyAnd(hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0],h,p,hillOn),hillInv(oldValue[nodeOrder[1]],nodeOrderInvert[1],h,p,hillOn))
			else:
				value=fuzzyOr(hillInv(oldValue[nodeOrder[0]],nodeOrderInvert[0],h,p,hillOn),hillInv(oldValue[nodeOrder[1]],nodeOrderInvert[1],h,p,hillOn))
			for i in range(2,len(nodeOrder)):
				if logicOperatorFlags[i]==0:
					value=fuzzyAnd(value,hillInv(oldValue[nodeOrder[i]],nodeOrderInvert[i],h,p,hillOn))
				else:
					value=fuzzyOr(value,hillInv(oldValue[nodeOrder[i]],nodeOrderInvert[i],h,p,hillOn))
		return value
	else:
		# print('origvalue')
		# print(oldValue[currentNode])
		return oldValue[currentNode]						
				
def runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out.  
	#list of triples which gives the total list of nodes. 
	counter=0;
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	seq=range(0,len(individualParse))
	nodeNoise=[]
	for node in range(0,len(individualParse)):
		nodeNoise.append(sigmaNode*random())
	for step in range(0,steps):
		shuffle(seq)
		for i in range(0,len(individualParse)):	
			newValue[seq[i]]=nodeNoise[seq[i]]+fuzzyUpdate(seq[i],newValue,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn)
		simData.append(list(newValue))
	avg= [0 for x in range(0,len(newValue))]
	stdev= [0 for x in range(0,len(newValue))]
	for step in range(0,len(simData)):
		for element in range(0,len(avg)):
			simData[step][element]=simData[step][element]+sigmaNetwork*random()
	for step in range(len(simData)-10,len(simData)):
		for element in range(0,len(avg)):
			avg[element]=avg[element]+simData[step][element]
	for element in range(0,len(avg)):
		avg[element]=avg[element]/10
	for step in range(len(simData)-10,len(simData)):
		for element in range(0,len(avg)):
			stdev[element]=stdev[element]+abs(simData[step][element]-avg[element])
		
	return (avg, stdev)
			
#this is the evaluation function. Need to create a graph then recreate everything with it.
def evaluate(individual, ss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn):
	RME=0
	boolValues, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,0,0)	
	for i in range(0, len(evaluateNodes)):
		RME=RME+(boolValues[evaluateNodes[i]]-ss[nodeList[evaluateNodes[i]]])**2
	return RME,

def simplifyNetwork(graph, ss):
#network simplification algorithm. 
# # 1. remove nodes which are neither part of input nor have input to them...
# # 2. remove straigth paths. 
# # 3. 

#  collapse straight lines
	print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1 and (len(graph.successors(x))==1)and (not (x in ss.keys())))]
	for rm in removeNodeList:
		before=graph.predecessors(rm)[0]
		after=graph.successors(rm)[0]
		edge1=graph.get_edge_data(before,rm)['signal']
		edge2=graph.get_edge_data(rm,after)['signal']
		inhCount=0
		if edge1=='i':
			inhCount=inhCount+1
		if edge2=='i':
			inhCount=inhCount+1
		if inhCount==1:
			graph.add_edge(before,after,signal='i')
		else:
			graph.add_edge(before,after,signal='a')
		graph.remove_node(rm)
		
		
# remove nodes with no predecessors or value in given steady state data
	print(len(graph.nodes()))
	newNodes = [x for x in graph.nodes() if not (len(graph.predecessors(x))==0 and (not (x in ss.keys())))]
	graph=graph.subgraph(newNodes)
	
	print(len(graph.nodes()))

def buildFatimaNetwork():
	dataFileName='C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv'
	
	#two dicts for the models
	nodeUpdateDict={}
	ss={} 
	data=loadFpkms(dataFileName)
	sortFpkms(data,ss)
	#print(ss.keys())
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)	
	
	currentfile='IL1b_pathways.txt'
	inputfile = open(currentfile, 'r')
	line = inputfile.read()
	codelist=re.findall('ko\d\d\d\d\d',line)	
	print(codelist)
	constructor.uploadKEGGcodes(codelist, graph, KEGGdict)
	for node in graph.nodes():
		if node in graph.successors(node):
			graph.remove_edge(node,node)
	nodeList=graph.nodes()
	
	simplifyNetwork(graph,ss)
	return graph, ss
	
def printSubgraphs(graph):
	nodeList=graph.nodes()
	subgraphnodes=set()
	for node in nodeList:
		if len(graph.predecessors(node))>15:
			subgraphnodes.add(node)
			for element in graph.predecessors(node):
				subgraphnodes.add(element)
			
			print('new node')
			print(node)
			print(len(graph.predecessors(node)))
			print(graph.predecessors(node))

			
def buildToolbox(individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues,simSteps,h,p,hillOn):
	# # #setup toolbox
	toolbox = base.Toolbox()

	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)

	# make a fitness minimization function
	creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
	# create a class of individuals that are lists from networkx
	creator.create("Individual", list, fitness=creator.FitnessMin)


	# how to create aliases for your individuals... 
	# toolbox.register("attr_float", random.random)
	# need this alias to create new graphs... but we should just be using the first one.... 

	toolbox.register("genRandomBitString", genRandBits, individualLength=individualLength)
	toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.individual)
	
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)
	hof = tools.HallOfFame(5, similar=numpy.array_equal)
	
	# finish registering the toolbox functions
	toolbox.register("mate", tools.cxTwoPoint)
	toolbox.register("mutate", tools.mutFlipBit, indpb=bitFlipProb)
	toolbox.register("select", tools.selNSGA2)
	toolbox.register("evaluate", evaluate, ss=ss, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, nodeOrderList=nodeOrderList, nodeOrderInvertList=nodeOrderInvertList, possibilityNumList=possibilityNumList, initValues=initValues, steps=simSteps, h=h,p=p,hillOn=hillOn)
	toolbox.register("similar", numpy.array_equal)
	return toolbox, hof, stats

def runFatimaGA():
	graph, ss = buildFatimaNetwork()
	simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb = returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValues, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	toolbox, hof, stats=buildToolbox(individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues,simSteps,h,p,hillOn)
	population=toolbox.population(n=popSize)
	algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)

	
def runFatimaSim(networkNoises, nodeNoises, trials):
	graph, ss=buildFatimaNetwork()
	simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb ,mutationProb= returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValues, evaluateNodes, possibilityNumList, nodeOrderInvertList=setupGAparams(graph, ss)
	#open file to save, use noises in name
	#save graph seperately
	toolbox, hof, stats=buildToolbox(individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues,simSteps,h,p,hillOn)
	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(toolbox, hof, stats, graph, trials, ss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValues,simSteps, h,p,hillOn,nodeNoise,networkNoise,popSize ,crossoverProb,mutationProb,generations )

		
		
def runFatimaSimTester(networkNoises, nodeNoises, trials):
	ss={}
	ss['zero']=1
	ss['one']=0
	ss['two']=0
	print(ss.keys())
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='a')
	graph.add_edge('one','two', signal='i')	
	simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb ,mutationProb= returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValues, evaluateNodes, possibilityNumList, nodeOrderInvertList=setupGAparams(graph, ss)
	#open file to save, use noises in name
	#save graph seperately
	toolbox, hof, stats=buildToolbox(individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues,simSteps,h,p,hillOn)
	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(toolbox, hof, stats, graph, trials, ss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValues,simSteps, h,p,hillOn,nodeNoise,networkNoise,popSize ,crossoverProb,mutationProb,generations )
		
			
def fuzzySimulator(toolbox, hof, stats, graph, trials, ss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValues, steps, h,p,hillOn,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations ):
	avgSSEs=[]
	hofs=[]
	truthIndividuals=[]
	truthAvgs=[]
	truthStdevs=[]
	for i in range(0,trials):
		newSS=copy.deepcopy(ss)
		individual=genRandBits(individualLength)
		boolValues, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,nodeNoise,networkNoise)
		for j in range(0,len(evaluateNodes)):
			newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
		toolbox.register("evaluate", evaluate, ss=newSS, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, nodeOrderList=nodeOrderList, nodeOrderInvertList=nodeOrderInvertList, possibilityNumList=possibilityNumList, initValues=boolValues, steps=steps, h=h,p=p,hillOn=hillOn)
		#toolbox.register("evaluate", evaluate, ss=ss, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, nodeOrderList=nodeOrderList, nodeOrderInvertList=nodeOrderInvertList, possibilityNumList=possibilityNumList, initValues=boolValues, steps=steps, h=h,p=p,hillOn=hillOn)
		population=toolbox.population(n=popSize)
		output=algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True)
		hofs.append(hof)
		truthIndividuals.append(individual)
		truthAvgs.append(boolValues)
		truthStdevs.append(stdev)
	prefix='fatimaTest_noise_'+str(nodeNoise)+'_networkNoise_'+str(networkNoise)
	outfile = open(prefix+'_hof.pkl', 'wb')
	pickle.dump(hof, outfile)
	outfile.close()
	print(hof)
	outfile = open(prefix+'_truth_individual.pkl', 'wb')
	pickle.dump(truthIndividuals, outfile)
	outfile.close()
	print(truthIndividuals)
	outfile = open(prefix+'_truth_avgs.pkl', 'wb')
	pickle.dump(truthAvgs, outfile)
	outfile.close()
	print(truthAvgs)
	outfile = open(prefix+'_truth_stdevs.pkl', 'wb')
	pickle.dump(truthStdevs, outfile)
	outfile.close()
	print(truthStdevs)
	
def testFuzzySim():
	filename='inputDataFatima.txt'
	KEGGfileName='ko04060.xml'
		
	#two dicts for the models
	nodeUpdateDict={}
	ss={}
	for i in range(1,5):
		ss[str(i)]=0
	ss['zero']=1
	ss['one']=0
	ss['two']=0
	print(ss.keys())
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='a')
	graph.add_edge('one','two', signal='i')	
	individualLength, individualParse, nodeList, nodeOrderList, initValues, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb=returnSimParams()
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	avg, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	
	
	sigmaNode=0
	sigmaNetwork=0
	individual=[0,1,0,0,0,0]
	print(individual)
	avg, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	sigmaNode=0
	sigmaNetwork=0
	individual=[1,0,0,0,0,0]
	print(individual)
	avg, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValues, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	
	
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='i')
	graph.add_edge('one','two', signal='i')	
	individualLength=setupGAparams(graph, ss)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, individualParse, nodeList, initValues))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, individualParse, nodeList, initValues))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, individualParse, nodeList, initValues))
	
	ss['zero']=1
	ss['one']=1
	ss['two']=1
	individualLength=setupGAparams(graph, ss)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, individualParse, nodeList, initValues))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, individualParse, nodeList, initValues))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, individualParse, nodeList, initValues))
	
if __name__ == '__main__':
	# ss={}
	# data=loadFpkms('C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv')
	# sortFpkms(data,ss)
	# print(ss)
	#runFatimaSim([.001,.005],[.001,.002],1)
	runFatimaSim([.001],[.001],1)
	#testFuzzySim()