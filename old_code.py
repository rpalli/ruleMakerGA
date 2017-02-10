#from scipy.stats import logistic
#import re
#from scoop import futures
#import pyximport; pyximport.install()
#import evaluator
#import cPickle as pickle


from random import random
import networkx as nx

def hill(x,h,p,hillOn): #applies hill function if called for
	if hillOn:
		return ((1+h**p)*x**p)/(h**p+x**p)		
	else:
		return x
	
def logRegPrepNet(graph, ss):
	newNodes = [x for x in graph.nodes() if  (not (x in ss.keys()))]
	for node in newNodes:
		befores=graph.predecessors(node)
		afters=graph.successors(node)
		for before in befores:
			for after in afters:
				edge1=graph.get_edge_data(before,node)['signal']
				edge2=graph.get_edge_data(node,after)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(before,after,signal='i')
				else:
					graph.add_edge(before,after,signal='a')
		graph.remove_node(node)
	print(len(graph.nodes()))



def buildFatimaNetwork(): #build a network from Fatimas data and from IL1 pathways in KEGG
	dataFileName='C:/Users/Rohith/Desktop/Rohith_data/ruleMaker_GA/Data/fatima_fpkms.csv'
	#dataFileName='/home/rohith/Documents/fatima_fpkms.csv'
	#two dicts for the models
	nodeUpdateDict={}
	data=loadFpkms(dataFileName)
	sss=sortFpkms(data)
	gc.collect()
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
	
	graph=simplifyNetwork(graph,sss[0])
	return graph, sss

			

def randomizeInputs(sss,samples): # randomzes extant inputs for steady state simulations
	newSS=[]
	for i in range(0,samples):
		newSS.append({})
	for key in sss[0].keys():
		for i in range(0,samples):
			newSS[i][key]=random()
	return newSS

#above are old functions and imports
	
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

def testSimConvergence(async, iters, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList, nodeNoise,networkNoise,earlyEvalNodes, genSteps):
	newSSS=[]
	trueNodeList=[]
	for i in range(0,trials):
		individual=genRandBits(individualLength)
		newSSS=[]
		for k in range(0,len(sss)):
			if async:
				boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[k], genSteps ,nodeNoise,networkNoise,earlyEvalNodes)
			else:
				boolValues=runFuzzyModel(individual, async,individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[k], genSteps,nodeNoise,networkNoise,earlyEvalNodes)
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(evaluateNodes)):
				newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
			newSSS.append(newSS)
		#update eval nodes
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		for values in newInitValueList:
			print(values)

# def fuzzyUpdate(currentNode,oldValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
# 	triple=individualParse[currentNode] #find locations of bits we need in individual
# 	inputOrder=inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
# 	inputOrderInvert=inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
# 	if possibilityNumList[currentNode]>0:
# 		logicOperatorFlags=individual[triple[1]:triple[2]] # find demarcations on individual for bits we need
# 		inputOrder=inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # determine what order of inputs is
# 		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
# 		if len(inputOrder)==0:
# 			value=oldValue[currentNode] #if no inputs, maintain value
# 		elif len(inputOrder)==1:
# 			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
# 			if individual[triple[0]]==1:
# 				value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
# 			else:
# 				value=oldValue[currentNode]
# 		else:
# 			#update nodes with more than one input
# 			counter =0
# 			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
# 			if logicOperatorFlags[0]==0:
# 				value=fuzzyAnd(Inv(oldValue[inputOrder[0]],inputOrderInvert[0]),Inv(oldValue[inputOrder[1]],inputOrderInvert[1]))
# 			else:
# 				value=fuzzyOr(Inv(oldValue[inputOrder[0]],inputOrderInvert[0]),Inv(oldValue[inputOrder[1]],inputOrderInvert[1]))
# 			for i in range(2,len(inputOrder)):
# 				#combines subsequent inputs
# 				if logicOperatorFlags[i-1]==0:
# 					value=fuzzyAnd(value,Inv(oldValue[inputOrder[i]],inputOrderInvert[i]))
# 				else:
# 					value=fuzzyOr(value,Inv(oldValue[inputOrder[i]],inputOrderInvert[i]))
# 		return value
# 	else:
# 		#returns savme value if now inputs
# 		return oldValue[currentNode]						
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
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb=returnSimParams()
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))

	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	
	
	sigmaNode=0
	sigmaNetwork=0
	individual=[0,1,0,0,0,0]
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	sigmaNode=0
	sigmaNetwork=0
	individual=[1,0,0,0,0,0]
	print(individual)
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	
	
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='i')
	graph.add_edge('one','two', signal='i')	
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	ss['zero']=1
	ss['one']=1
	ss['two']=1
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))


def testFuzzySimTimeCourse():

	csvfile=open('output2.csv', 'w')
	writer = csv.writer(csvfile)
    
	simCount=0	
	#two dicts for the models
	nodeUpdateDict={}
	ss={}
	sss=[]
	for i in range(1,5):
		ss[str(i)]=0
	ss['zero']=1
	ss['one']=0
	ss['two']=0
	ss['three']=0
	print(ss.keys())
	sss.append(ss)
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='a')
	graph.add_edge('one','two', signal='i')	
	graph.add_edge('two','three', signal='a')
	graph.add_edge('three','one', signal='a')	

	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes= setupGAparams(graph, sss)
	async, iters, simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb=returnSimParams()
	print(individualLength)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async,individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

	print("with_noise .01")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)





	
	sigmaNode=0
	sigmaNetwork=0
	individual=[0,1,0,0,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)


	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,1,0,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,0,1,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	


	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,0,0,1,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	




	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,0,0,0,1,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

#above are old testing functions, below are current testing functions

def LiuNetwork1Builder():
	graph = nx.DiGraph()
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='i')	
	graph.add_edge('f','k', signal='i')
	graph.add_edge('a','c', signal='a')
	graph.add_edge('b','d', signal='a')
	graph.add_edge('c','f', signal='a')	
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')

	return graph
def LiuSSS1Builder():
	sss=[]
	ss={}
	ss['a']=1
	ss['b']=0
	ss['c']=0
	ss['d']=0
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=1
	ss['c']=0
	ss['d']=0
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=0
	ss['c']=1
	ss['d']=0
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=0
	ss['c']=0
	ss['d']=1
	ss['f']=0
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	sss.append(ss)
	ss={}
	ss['a']=0
	ss['b']=0
	ss['c']=0
	ss['d']=0
	ss['f']=1
	ss['g']=0
	ss['h']=0
	ss['j']=0
	ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=1
	# ss['h']=0
	# ss['j']=0
	# ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=1
	# ss['j']=0
	# ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=0
	# ss['j']=1
	# ss['k']=0
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=0
	# ss['j']=0
	# ss['k']=1
	# sss.append(ss)
	# ss={}
	# ss['a']=1
	# ss['b']=1
	# ss['c']=1
	# ss['d']=1
	# ss['f']=1
	# ss['g']=1
	# ss['h']=1
	# ss['j']=1
	# ss['k']=1
	# sss.append(ss)
	# ss={}
	# ss['a']=0
	# ss['b']=0
	# ss['c']=0
	# ss['d']=0
	# ss['f']=0
	# ss['g']=0
	# ss['h']=0
	# ss['j']=0
	# ss['k']=0
	# sss.append(ss)
	return sss
def LiuNetwork3Builder():
	graph = nx.DiGraph()
	
	graph.add_edge('EGFR','RAS', signal='a')
	graph.add_edge('EGFR','PI3K', signal='a')
	graph.add_edge('TNFR','FADD', signal='a')	
	graph.add_edge('TNFR','PI3K', signal='a')
	graph.add_edge('TNFR','RIP1', signal='a')	
	graph.add_edge('DNA_DAMAGE','ATM', signal='a')
	graph.add_edge('RAS','ERK', signal='a')	
	graph.add_edge('RAS','JNK', signal='a')
	graph.add_edge('FADD','CASP8', signal='a')	
	graph.add_edge('PI3K','AKT', signal='a')
	graph.add_edge('RIP1','P38', signal='a')
	graph.add_edge('RIP1','JNK', signal='a')
	graph.add_edge('ATM','CHK', signal='a')	
	graph.add_edge('ATM','P53', signal='a')	
	graph.add_edge('ERK','RAS', signal='a')
	graph.add_edge('ERK','MYC', signal='a')	
	graph.add_edge('ERK','BID', signal='i')
	graph.add_edge('CASP8','BID', signal='a')	
	graph.add_edge('CASP8','CASP3', signal='a')
	graph.add_edge('CASP8','RIP1', signal='i')
	graph.add_edge('JNK','P53', signal='a')	
	graph.add_edge('JNK','FOXO1', signal='a')
	graph.add_edge('AKT','BID', signal='i')	
	graph.add_edge('AKT','CASP8', signal='i')
	graph.add_edge('AKT','BIM', signal='i')
	graph.add_edge('AKT','4EBP1', signal='i')	
	graph.add_edge('AKT','S6K', signal='a')
	graph.add_edge('AKT','P53', signal='i')	
	graph.add_edge('AKT','FOXO1', signal='i')
	graph.add_edge('P38','CDC25', signal='a')	
	graph.add_edge('P38','FOXO1', signal='a')
	graph.add_edge('P38','P53', signal='a')
	graph.add_edge('CHK','CDC25', signal='i')	
	graph.add_edge('CHK','P53', signal='a')	
	graph.add_edge('CDC25','CDK', signal='a')
	graph.add_edge('P53','PUMA', signal='a')
	graph.add_edge('P53','FADD', signal='a')
	graph.add_edge('P53','ATM', signal='a')	
	graph.add_edge('P53','CHK', signal='i')
	graph.add_edge('P53','CDK', signal='i')
	graph.add_edge('P53','CYCLIN', signal='i')	
	graph.add_edge('P53','BIM', signal='a')
	graph.add_edge('FOXO1','P27', signal='a')
	graph.add_edge('BIM','BAX', signal='a')	
	graph.add_edge('BID','BAX', signal='a')
	graph.add_edge('MYC','PROLIFERATION', signal='a')	
	graph.add_edge('PUMA','BAX', signal='a')
	graph.add_edge('S6K','S6', signal='a')	
	graph.add_edge('S6K','PI3K', signal='a')
	graph.add_edge('P27','CYCLIN', signal='i')
	graph.add_edge('P27','CELL_CYCLE', signal='a')	
	graph.add_edge('BAX','SMAC', signal='a')
	graph.add_edge('SMAC','CASP3', signal='a')	
	graph.add_edge('4EBP1','PROLIFERATION', signal='i')
	graph.add_edge('S6','PROLIFERATION', signal='a')	
	graph.add_edge('CYCLIN','CDK', signal='a')
	graph.add_edge('CYCLIN','PROLIFERATION', signal='a')
	graph.add_edge('CDK','PROLIFERATION', signal='a')	
	graph.add_edge('CHK','CELL_CYCLE', signal='a')
	graph.add_edge('BAX','CASP9', signal='a')	
	graph.add_edge('CASP9','CASP3', signal='a')
	graph.add_edge('CASP3','APOPTOSIS', signal='a')	
	return graph

def LiuNetwork3SimTest(inputNum,trialNum):
	graph = LiuNetwork3Builder()
	#runNoiseSimTest(graph, sss, [.001,.002,.005,.01,.02,.05,.1,.2,.5], [0], 10)
	sss1=synthesizeInputs(graph,inputNum)
	runNoiseSimTest(graph, sss1, [0], [0], trialNum)

def LiuNetwork1SimTest(inputNum,trialNum):
	#two dicts for the models
	#sss=LiuSSS1Builder()
	graph = LiuNetwork1Builder()
	#runNoiseSimTest(graph, sss, [.001,.002,.005,.01,.02,.05,.1,.2,.5], [0], 10)
	sss1=synthesizeInputs(graph,inputNum)
	runNoiseSimTest(graph, sss1, [.001], [0], trialNum)

def LiuNetwork2SimTest(inputNum,trialNum):
	#two dicts for the models
	sss, graph=LiuNetwork2Builder()
	sss1=synthesizeInputs(graph,inputNum)
	runNoiseSimTest(graph, sss1, [0], [0], trialNum)

def LiuNetwork2Builder():
	graph = nx.DiGraph()
	graph.add_edge('igf1','ras', signal='a')
	graph.add_edge('tgfa','ras', signal='a')
	graph.add_edge('igf1','pi3k', signal='a')
	graph.add_edge('tgfa','pi3k', signal='a')
	graph.add_edge('ras','pi3k', signal='a')
	graph.add_edge('ras','map3k1', signal='a')
	graph.add_edge('ras','mek12', signal='a')
	graph.add_edge('tnfa','pi3k', signal='a')
	graph.add_edge('tnfa','jnk12', signal='a')
	graph.add_edge('tnfa','map3k1', signal='a')
	graph.add_edge('tnfa','map3k7', signal='a')
	graph.add_edge('tnfa','mkk4', signal='a')
	graph.add_edge('il1a','map3k7', signal='a')
	graph.add_edge('il1a','map3k1', signal='a')
	graph.add_edge('map3k7','ikk', signal='a')
	graph.add_edge('map3k7','mkk4', signal='a')
	graph.add_edge('map3k7','p38', signal='a')
	graph.add_edge('map3k7','hsp27', signal='a')
	graph.add_edge('map3k1','ikk', signal='a')
	graph.add_edge('map3k1','jnk12', signal='a')
	graph.add_edge('map3k1','mkk4', signal='a')
	graph.add_edge('pi3k','map3k1', signal='a')
	graph.add_edge('pi3k','akt', signal='a')
	graph.add_edge('pi3k','mek12', signal='a')
	graph.add_edge('akt','mek12', signal='a')
	graph.add_edge('akt','ikk', signal='a')
	graph.add_edge('mek12','erk12', signal='a')
	graph.add_edge('ikk','ikb', signal='a')
	graph.add_edge('mkk4','jnk12', signal='a')
	graph.add_edge('mkk4','p38', signal='a')
	graph.add_edge('erk12','hsp27', signal='a')
	graph.add_edge('p38','hsp27', signal='a')

	# listy=['igf1','tgfa','ras','tnfa','il1a','map3k7','pi3k','map3k1','akt','mek12','ikk','mkk4','erk12','jk12','ikb','p38','hsp27']
	# sss=[]
	# for i in range(0,10):
	# 	sss.append({})
	# for element in listy:
	# 	for i in range(0,10):
	# 		sss[i][element]=random()
	return  graph

def runNoiseSimTest(graph, sss, networkNoises, nodeNoises, trials):
	async, iters, steps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb,mu,lambd, complexPenalty, genSteps= returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss)
	#open file to save, use noises in name
	#save graph seperately
	toolbox, stats=buildToolbox(individualLength, bitFlipProb, len(sss))
	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList, steps, h,p,hillOn,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps)




def runFatimaGA():
	#very simple function to just run GA from fatima data with no frills. just see output from a single call to eaMuCommaLambda
	graph, ss = buildFatimaNetwork()
	iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb = returnSimParams()
	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList= setupGAparams(graph, ss)
	toolbox, hof, stats=buildToolbox(iters,individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList,simSteps, )
	toolbox.register("evaluate", evaluate, iters=iters,ss=ss,  evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValueList=initValues, steps=simSteps, h=h,p=p,hillOn=hillOn)
	population=toolbox.population(n=popSize)
	algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)

def bootstrapRules(graph,sss,reSampleNum,sampleSize, counter, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList):
	ruleList=[]
	for i in range(0,reSampleNum):
		sss2=numpy.random.choice(sss, size=sampleSize, replace=True, p=None)
		ruleList.append(propBoolModel( individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList,graph,sss2))
	ruleDict=[]

	# for i in range(0,len(nodeList)):
	# 	ruleDict.append({})
	# for trialNum in range(0,len(ruleList)):
	# 	for nodeNum in range(0,len(ruleList[trialNum])):
	# 		ruleDict[nodeNum][ruleList[trialNum][nodeNum]]=0
	# for trialNum in range(0,len(ruleList)):
	# 	for nodeNum in range(0,len(ruleList[trialNum])):
	# 		ruleDict[nodeNum][ruleList[trialNum][nodeNum]]=ruleDict[nodeNum][ruleList[trialNum][nodeNum]]+1
	return ruleList

def drawGraph3(representative, filename, KEGGdict):
	rep=representative.copy()
	dictionary={}
	names=nx.get_node_attributes(representative, 'name')
	for n in rep.nodes():
		if len(names[n].split())==1:
			if names[n] in KEGGdict.keys():
				dictionary[n]=KEGGdict[names[n]]
				#print(rep.node[n]['name'])
		else :
			translated=''
			for word in names[n].split():
				word1=word.lstrip('ko:')
				word1=word1.lstrip('gl:')
				if word1 in KEGGdict.keys():
					translated=translated+KEGGdict[word1]+'-'
				else:
					translated=translated+word1+'-'
			dictionary[n]=translated
	repar= nx.relabel_nodes(rep,dictionary)
	#print(repar.nodes())
	#print(repar.edges())
	B=agraph.to_agraph(repar)        # convert to a graphviz graph\
	B.layout()            # neato layout
	B.draw(filename)       # write postscript in k5.ps with neato layout

def drawGraph2(representative, filename):
	B=agraph.to_agraph(representative)        # convert to a graphviz graph\
	B.layout()            # neato layout
	B.draw(filename)       # write postscript in k5.ps with neato layout

def readKEGG(lines, graph, KEGGdict):
	#the network contained in a KEGG file to the graph. 
	grouplist=[] #sometimes the network stores a group of things all activated by a single signal as a "group". We need to rewire the arrows to go to each individual component of the group so we save a list of these groups and a dictionary between id numbers and lists of elements that are part of that group
	groups={}
	i=0 #i controls our movement through the lines of code. 
	dict={} #identifies internal id numbers with names (much like the code earlier does for groups to id numbers of elements of the group)
	while i< len(lines):
		line=lines[i]
		# print line
		#add nodes
		if "<entry" in line:		
			nameline=line.split('name="')
			namestring=nameline[1]
			nameline=namestring.split('"')
			name=str.lstrip(nameline[0],'ko:')
			name=str.lstrip(name,'gl:')
			name=name.replace(' ','-')
			name=name.replace('ko:','')
			name=name.replace(':','-')
			
			# print(name)
			
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			#build a dictionary internal to this file
			idline=line.split('id="')
			idstring=idline[1]
			idline=idstring.split('"')
			id=idline[0]
			
			
			namelist=[]
			if '-' in name: 
				namelines=name.split('-')
				for x in namelines:
					if x in KEGGdict.keys():
						namelist.append(KEGGdict[x])
					else:
						namelist.append(x)
				name=namelist[0]+'-'
				for x in range(1,len(namelist)):
					name=name+namelist[x]+'-'
				name=name[:-1:]
			if name in KEGGdict.keys():
				dict[id]=KEGGdict[name]
			else:
				dict[id]=name
			#print(dict[id])
			#if the node is a group, we create the data structure to deconvolute the group later
			if type=='group':
				i=i+1
				newgroup=[]
				grouplist.append(id)
				while(not '</entry>' in lines[i]):
					if 'component' in lines[i]:
						newcompline=lines[i].split('id="')
						newcompstring=newcompline[1]
						newcompline=newcompstring.split('"')
						newcomp=newcompline[0]
						newgroup.append(newcomp)
					i=i+1
				groups[id]=newgroup
			else:
				graph.add_node(dict[id],{'name':dict[id],'type':type})
				#print(name)
				#print(id)
				i=i+1
		#add edges from KEGG file
		elif "<relation" in line:
			color='black'
			signal='a'
			subtypes=[]
			
			entryline=line.split('entry1="')
			entrystring=entryline[1]
			entryline=entrystring.split('"')
			entry1=entryline[0]
			
			entryline=line.split('entry2="')
			entrystring=entryline[1]
			entryline=entrystring.split('"')
			entry2=entryline[0]
			

			
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			while(not '</relation>' in lines[i]):
				if 'subtype' in lines[i]:
					nameline1=lines[i].split('name="')
					namestring1=nameline1[1]
					nameline1=namestring1.split('"')
					subtypes.append(nameline1[0])
				i=i+1
			
			
			#color and activation assignament based on the type of interaction in KEGG
			if ('activation' in subtypes) or ('expression' in subtypes):
				color='green'
				signal='a'
			elif 'inhibition' in subtypes:
				color='red'
				signal='i'
			elif ('binding/association' in subtypes) or('compound' in subtypes):
				color='purple'
				signal='a'
			elif 'phosphorylation' in subtypes:
				color='orange'
				signal='a'
			elif 'dephosphorylation' in subtypes:
				color='pink'
				signal='i'
			elif 'indirect effect' in subtypes:
				color='cyan'
				signal='a'
			elif 'dissociation' in subtypes:
				color='yellow'
				signal='i'
			elif 'ubiquitination' in subtypes:
				color='cyan'
				signal='i'
			else:
				print('color not detected. Signal assigned to activation arbitrarily')
				print(subtypes)
				signal='a'
			
			#this code deconvolutes the "groups" that KEGG creates 
			#each element of the group gets an edge to or from it instead of the group
			if entry1 in grouplist:
				for current1 in groups[entry1]:
					node1=dict[current1]
					if entry2 in grouplist:
						for current2 in groups[entry2]:
							node2=dict[current2]
							graph.add_edge(node1, node2, color=color, subtype=name, type=type, signal=signal )
					else:
						node2=dict[entry2]
						graph.add_edge(node1, node2, color=color,  subtype=name, type=type, signal=signal )
			elif entry2 in grouplist:
				node1=dict[entry1]
				for current in groups[entry2]:
					if current in grouplist:
						for current1 in groups[current]:
							if current in dict.keys():
								node2=dict[current]
								graph.add_edge(node1,node2, color=color,  subtype=name, type=type, signal=signal )
							else: 
								print('groups on groups on groups')
						#print('groups on groups')
					elif current in dict.keys():
						node2=dict[current]
						graph.add_edge(node1,node2, color=color, subtype=name, type=type, signal=signal )
					else:
						print('not in the keys')
						print(current)
			else:
				node1=dict[entry1]
				node2=dict[entry2]
				
				graph.add_edge(node1,node2, color=color, subtype=name, type=type, signal=signal )
		i=i+1
def read_biopax(lines, graph):
	soup = BeautifulSoup(''.join(lines), 'xml')

	##old code, potentially useless
	# pathways = soup.find_all('pathway')
	# for pathway in pathways:
	# 	name = pathway.find('NAME').string
	# 	synonyms = [s.string for s in pathway.find_all('SYNONYMS')]
	# 	component_ids = [c.get('rdf:resource')[1:] for c in pathway.find_all('PATHWAY-COMPONENTS')]
	# 	print name
	# 	for c_id in component_ids:
	# 		component = soup.find(attrs={'rdf:ID': c_id})
	# 		print component.prettify()
	##end old code


	#treating edge classes as nodes for now, will simplify later
	#constuct a BeautifulSoup4 parser object with the libxml backend
	#input in lines is given as an array, the join function converts to 1 long string
	#On windows this fails, try 'lxml' instead of 'lxml-xml'
	soup = BeautifulSoup(''.join(lines), 'xml')

	#in BioPAX, both physical objects (proteins, genes, etc.) and interactions (activation, 
	#inhibition, etc.) are represented as NODES.  The edges in the system represent membership 
	#or participation in these interactions.  Accordingly, we parse all objects into nodes below.

	#for every biopax class name (defined on lines 29-31)

## Rohith's OLD version, kept for history
# def parseKEGGdict(filename):
# 	#makes a dictionary to convert ko numbers from KEGG into real gene names
# 	#this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two.
# 	dict={}
# 	inputfile = open(filename, 'r')
# 	lines = inputfile.readlines()
# 	for line in lines:
# 		if line[0]=='D':
# 			kline=line.split('      ')
# 			kstring=kline[1]
# 			kline=kstring.split('  ')
# 			k=kline[0]
# 			nameline=line.replace('D      ', 'D')
# 			nameline=nameline.split('  ')
# 			namestring=nameline[1]
# 			nameline=namestring.split(';')
# 			name=nameline[0]
# 			dict[k]=name
# 	return dict

def runFatimaSim(networkNoises, nodeNoises, trials):
	#simulate data on Fatima network. Output how well it goes. 
	samples=20
	graph, sss=buildFatimaNetwork()
	sss1=synthesizeInputs(graph,samples)
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb ,mutationProb, mu,lambd= returnSimParams()
	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss1)
	print("individualLength")
	print(individualLength)
	#open file to save, use noises in name
	#save graph seperately
	
	toolbox, stats=buildToolbox(individualLength, bitFlipProb, samples)

	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss1,individualLength,evaluateNodes, individualParse,nodeList,inputOrderList,inputOrderInvertList,possibilityNumList,initValueList,simSteps ,nodeNoise,networkNoise,popSize ,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd)
		
		

#if __name__ == '__main__':
	# ss={}
	#data=loadFpkms('C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv')
	#sortFpkms(data,ss)
	# print(ss)
	#teseter.LiuNetwork2SimTest()
	#runFatimaSim([.01,.05],[.01,.02],1)
	#runFatimaSim([0],[0.001],10)
	#tester.LiuNetwork2SimTest(20,1)
	# inputNum=50
	# counters=[]
	# graphStart=tester.LiuNetwork2Builder()
	# sss=synthesizeInputs(graphStart,inputNum)
	# graph=simplifyNetwork(graphStart,sss[0])
	# nx.write_graphml(graph,'Liu2.graphml')
	# for i in range(0,1):
	# 	sss=synthesizeInputs(graph,inputNum)
	# 	counters.append(testBootstrap(5,10,graph,sss))
	# print(counters)

def setupBootstrapParams(graph):
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	
	evaluateNodes=[] 
	individualParse=[] 
	inputOrderList=[] 
	inputOrderInvertList=[]
	possibilityNumList=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 100
	counter=0
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
		inputOrderList.append(possibilities)
		inputOrderInvertList.append(activities)
		possibilityNumList.append(len(possibilities))
		if len(possibilities)==0:
			logNum=0
		else:
			logNum=ceil(log(len(possibilities))/log(2))
		individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(possibilities)-1)])
		counter=counter+logNum+len(possibilities)
		# print(counter)
		# print(individualParse)
	return counter, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList

def testBootstrap(reSampleNum,sampleSize, graph,sss):

	# graph=nx.DiGraph()
	# graph.add_edge('1','3',signal='a')
	# graph.add_edge('2','3',signal='a')

	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes= setupGAparams(graph, sss)
	nodeNoise=0
	networkNoise=0
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps=returnSimParams()
	#for i in range(0,trials):
	individual=genRandBits(individualLength)
	# print(individualLength)
	# print(individual)
	
	truths=[]
	#print(writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
	for i in range(0,len(nodeList)):
		truths.append(writeFuzzyNode(i,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))		
	newSSS=[]
	for k in range(0,len(sss)):
		if async:
			boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps,nodeNoise,networkNoise,earlyEvalNodes)
		else:
			boolValues=runFuzzyModel(individual, async,individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps ,nodeNoise,networkNoise,earlyEvalNodes)
		newSS=copy.deepcopy(sss[k])
		for j in range(0,len(evaluateNodes)):
			newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
		newSSS.append(newSS)
		#update eval nodes
		#print(evaluate(individual, async, iters, newSSS, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList, steps, ,earlyEvalNodes,complexPenalty))
	
	initValueList=[]
	for j in range(0,len(newSSS)):
		initValueList.append([])
	for i in range(0,len(nodeList)):
		for j in range(0,len(newSSS)):
			ss=newSSS[j]
			if  nodeList[i] in newSSS[0].keys():
				initValueList[j].append(ss[nodeList[i]])
			else:
				initValueList[j].append(0.5)
				
	# print(individual)
	for i in range(0,len(nodeList)):
		# print('new Node')
		# print(individual)
		# print(individualParse[i])
		deviation=0
		# for steadyStateNum in range(0,len(newSSS)):
		# 	derivedVal=propUpdate(i,initValueList[steadyStateNum],individual[individualParse[i][0]:individualParse[i][2]], individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList, )
		# 	deviation=deviation+(derivedVal-newSSS[steadyStateNum][nodeList[i]])**2
		# print(deviation)

	ruleList=bootstrapRules(graph,newSSS,reSampleNum,sampleSize, individualLength, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList)
	wrongCount=0
	correctCount=0
	for item in ruleList:
		for i in range(0,len(nodeList)):
			if possibilityNumList[i]>1:
				if truths[i]==writeFuzzyNode(i,item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
					correctCount=correctCount+1
				else:
					wrongCount=wrongCount+1
					print(truths[i])
					print(writeFuzzyNode(i,item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))

	return wrongCount, correctCount
		# print(item)
		# print(writeFuzzyModel(item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))



def fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,inputOrderList,inputOrderInvertList,possibilityNumList,initValueList, steps ,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps):
	hofSize=10
	avgSSEs=[]
	hofs=[]
	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	truthCounter=[0 for number in xrange(hofSize)]
	trueNodeList=[]
	for i in range(0,trials):
		individual=genRandBits(individualLength)
		newSSS=[]
		for k in range(0,len(sss)):
			if async:
				boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps,nodeNoise,networkNoise,earlyEvalNodes)
			else:
				boolValues=runFuzzyModel(individual, async,individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps ,nodeNoise,networkNoise,earlyEvalNodes)
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(evaluateNodes)):
				newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
			newSSS.append(newSS)
		#update eval nodes
		print(evaluate(individual, async, iters, newSSS, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList, steps,earlyEvalNodes,complexPenalty))
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		#print(len(sss))
		toolbox.register("evaluate", evaluate, iters=iters, async=async, sss=newSSS, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValueList=newInitValueList, steps=steps, earlyEvalNodes=earlyEvalNodes, complexPenalty=complexPenalty)
		#toolbox.register("evaluate", evaluate, ss=ss, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValues=boolValues, steps=steps, h=h,p=p,hillOn=hillOn)
		population=toolbox.population(n=popSize)
		hof = tools.HallOfFame(hofSize, similar=numpy.array_equal)
		output = tools.Logbook()
		output=algo.eaMuCommaLambda(population, toolbox, mu=mu, lambda_=lambd, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)
		truth=(writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
		
		for j in range(0,10):
			bestRun=(writeFuzzyModel(hof[j], individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
			if truth==bestRun:
				truthCounter[j]+=1
				break
			
			elif j==0:
				truthlines=truth.split('\n')
				newlines=bestRun.split('\n')
				trueNodes=0
				for k in range(0,len(truthlines)):
					if truthlines[k]==newlines[k]:
						trueNodes=trueNodes+1
					else:
						print("incorrect pair: true then test")
						print(truthlines[k])
						print(newlines[k])
				trueNodeList.append(trueNodes)
		# else:
			# print('WRONG')
			# print(truth)
			# print(bestRun)
		avgs=[output[1][k]['min'] for k in range(0,len(output[1]))]
		#print(avgs)
		hofs.append(hof)
		temp=[]
		for hofind in range(0,len(hof)):
			tempVal=0
			for bit in range(0,len(hof[hofind])):
				if not hof[hofind][bit]==individual[bit]:
					tempVal=tempVal+1
			temp.append(tempVal)
		hofScores.append(temp)
		truthIndividuals.append(individual)
		truthAvgs.append(boolValues)
		gc.collect()
	# prefix='fatimaTest_noise_'+str(nodeNoise)+'_networkNoise_'+str(networkNoise)
	# outfile = open(prefix+'_hof.pkl', 'wb')
	# pickle.dump(hof, outfile)
	# outfile.close()
	f = open('hof_differences_'+str(nodeNoise)+str(networkNoise)+'.txt', 'w')
	g = open('hof_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	h = open('truth_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	f.write(str(hofScores))
	f.write('\n')
	for hofy in range(0,len(hofs)):
		g.write(str(hofs[hofy]))
		g.write('\n')
	h.write(str(truthIndividuals))
	h.write('\n')
	f.close()
	g.close()
	h.close()
	print('# true of # trials')
	for k in range(0,len(truthCounter)):
		print(truthCounter[k])
	print(trials)
	print("for the incorrects, by node")
	print(trueNodeList)
	print(len(truthlines))


def runFuzzyModel( individual, async, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValues, steps,sigmaNode, sigmaNetwork,earlyEvalNodes):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out.  
	#list of triples which gives the total list of nodes. 
	counter=0;
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))
	seq=range(0,len(individualParse))
	earlyEvalNodes=[]
	for node in range(0,len(earlyEvalNodes)):
		newValue[earlyEvalNodes[node]]=propUpdateNet(earlyEvalNodes[node],newValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
	for step in range(0,steps):
		oldValue=list(newValue)
		if async:
			shuffle(seq)
		for i in range(0,len(individualParse)):	
			if async:
				temp=propUpdateNet(seq[i],newValue,individual, individualParse[seq[i]], nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
			else:
				temp=propUpdateNet(seq[i],oldValue,individual, individualParse[seq[i]], nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)

			if sigmaNode==0:
				newValue[seq[i]]=temp
			else:
				newValue[seq[i]]=temp+sigmaNode*(1-temp/(.1+temp)-temp/10)*random()
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
	#for step in range(len(simData)-10,len(simData)):
		#for element in range(0,len(avg)):
			#stdev[element]=stdev[element]+(simData[step][element]-avg[element])^2
	#return (avg, stdev)
	return avg

def iterateFuzzyModel(async, iters, individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValues, steps ,sigmaNode, sigmaNetwork,earlyEvalNodes):
	sum=[0 for x in range(0,len(initValues))]
	for i in range(0,iters):
		avg,stdev=runFuzzyModel( individual, async, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValues, steps, sigmaNode, sigmaNetwork, earlyEvalNodes)
		for j in range(0,len(sum)):
			sum[j]=sum[j]+avg[j]
	# print(sum)
	avgs=list(sum)
	for i in range(0,len(sum)):
		avgs[i]=sum[i]/float(iters)
	return avgs
def propBoolModel(individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList,graph,sss1):
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps=returnSimParams()
	bestList=[]
	initValueList=[]
	for j in range(0,len(sss1)):
		initValueList.append([])
	for i in range(0,len(nodeList)):
		for j in range(0,len(sss1)):
			ss=sss1[j]
			if  nodeList[i] in sss1[0].keys():
				initValueList[j].append(ss[nodeList[i]])
			else:
				initValueList[j].append(0.5)
	for i in range(0,len(nodeList)):
		print(nodeList[i])
		print(individualParse[i])
		currentDev=10*len(sss1)
		best=[]
		if possibilityNumList[i]>0:
			for j in range(0,int(2**int(individualParse[i][2]-individualParse[i][0]+1))):
				#if possibilityNumList[i]==1:
					# print('single possibility')
					# print(j)
				bits=[]
				bits=bitList(j)
				while len(bits)<(individualParse[i][2]-individualParse[i][0]+1):
					bits.insert(0,0)
				deviation=0
				for steadyStateNum in range(0,len(sss1)):
					derivedVal=propUpdateNode(i,initValueList[steadyStateNum],bits, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList, )
					deviation=deviation+(derivedVal-sss1[steadyStateNum][nodeList[i]])**2
				print(bits)
				print(deviation)
				if(deviation<currentDev):
					print("best")
					best=bits
					currentDev=deviation	
		# print(i)
		# print(nodeList[i])
		bestList.append(best)
		# print(best)
	return [item for sublist in bestList for item in sublist]

def propUpdateNet(currentNode,oldValue,individual, triple, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList,modelType,corrMat):
	inputOrder=inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
	inputOrderInvert=inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
		inputOrder=list(inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]) # determine what order of inputs is
		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
		if len(inputOrder)==0:
			value=oldValue[currentNode] #if no inputs, maintain value
			return value
		elif len(inputOrder)==1:
			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
			if individual[triple[0]]==1:
				value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
			else:
				value=oldValue[currentNode]
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
					if(modelType=='prop'):
						tempVal=propAndE1(upstreamVals[counter],upstreamVals[counter+1],inputOrder[counter], inputOrder[counter+1],corrMat)
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
				tempVal=propOrE1(upstreamVals.pop(0),upstreamVals.pop(0))
				upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
			return upstreamVals[0]
	else:
		#returns savme value if now inputs
		return oldValue[currentNode]						

def propUpdateNode(currentNode,oldValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	
	oldTriple=individualParse[currentNode]
	triple=[]
	triple.append(0)
	triple.append(oldTriple[1]-oldTriple[0])
	triple.append(oldTriple[2]-oldTriple[0])

	retinue = propUpdateNet(currentNode,oldValue,individual, triple, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
	return retinue
					

def propOrE1(num1, num2, index1, index2,corrMat):
	return (num1+num2)*(2-corrMat[index1][index2])/2

def propAndE1(num1,num2,index1,index2,corrMat):
	return (num1+num2)*corrMat[index1][index2]/2

def setupGAparams(graph, sss):
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
			individualParse.append([int(counter),int(counter+logNum-1),int(counter+logNum+len(preds)-2)]) #set up limits: lognum-1 for list then len(preds)-1 for the activity flags
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
	return counter-1, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes
	#from scipy.stats import logistic
#import re
#from scoop import futures
#import pyximport; pyximport.install()
#import evaluator
#import cPickle as pickle
import networkx as nx

def hill(x,h,p,hillOn): #applies hill function if called for
	if hillOn:
		return ((1+h**p)*x**p)/(h**p+x**p)		
	else:
		return x
	
def logRegPrepNet(graph, ss):
	newNodes = [x for x in graph.nodes() if  (not (x in ss.keys()))]
	for node in newNodes:
		befores=graph.predecessors(node)
		afters=graph.successors(node)
		for before in befores:
			for after in afters:
				edge1=graph.get_edge_data(before,node)['signal']
				edge2=graph.get_edge_data(node,after)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(before,after,signal='i')
				else:
					graph.add_edge(before,after,signal='a')
		graph.remove_node(node)
	print(len(graph.nodes()))



def buildFatimaNetwork(): #build a network from Fatimas data and from IL1 pathways in KEGG
	dataFileName='C:/Users/Rohith/Desktop/Rohith_data/ruleMaker_GA/Data/fatima_fpkms.csv'
	#dataFileName='/home/rohith/Documents/fatima_fpkms.csv'
	#two dicts for the models
	nodeUpdateDict={}
	data=loadFpkms(dataFileName)
	sss=sortFpkms(data)
	gc.collect()
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
	
	graph=simplifyNetwork(graph,sss[0])
	return graph, sss

			

def randomizeInputs(sss,samples): # randomzes extant inputs for steady state simulations
	newSS=[]
	for i in range(0,samples):
		newSS.append({})
	for key in sss[0].keys():
		for i in range(0,samples):
			newSS[i][key]=random()
	return newSS

#above are old functions and imports
	
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

def testSimConvergence(async, iters, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList, nodeNoise,networkNoise,earlyEvalNodes, genSteps):
	newSSS=[]
	trueNodeList=[]
	for i in range(0,trials):
		individual=genRandBits(individualLength)
		newSSS=[]
		for k in range(0,len(sss)):
			if async:
				boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[k], genSteps ,nodeNoise,networkNoise,earlyEvalNodes)
			else:
				boolValues=runFuzzyModel(individual, async,individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[k], genSteps,nodeNoise,networkNoise,earlyEvalNodes)
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(evaluateNodes)):
				newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
			newSSS.append(newSS)
		#update eval nodes
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		for values in newInitValueList:
			print(values)

# def fuzzyUpdate(currentNode,oldValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
# 	triple=individualParse[currentNode] #find locations of bits we need in individual
# 	inputOrder=inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
# 	inputOrderInvert=inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
# 	if possibilityNumList[currentNode]>0:
# 		logicOperatorFlags=individual[triple[1]:triple[2]] # find demarcations on individual for bits we need
# 		inputOrder=inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # determine what order of inputs is
# 		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
# 		if len(inputOrder)==0:
# 			value=oldValue[currentNode] #if no inputs, maintain value
# 		elif len(inputOrder)==1:
# 			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
# 			if individual[triple[0]]==1:
# 				value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
# 			else:
# 				value=oldValue[currentNode]
# 		else:
# 			#update nodes with more than one input
# 			counter =0
# 			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
# 			if logicOperatorFlags[0]==0:
# 				value=fuzzyAnd(Inv(oldValue[inputOrder[0]],inputOrderInvert[0]),Inv(oldValue[inputOrder[1]],inputOrderInvert[1]))
# 			else:
# 				value=fuzzyOr(Inv(oldValue[inputOrder[0]],inputOrderInvert[0]),Inv(oldValue[inputOrder[1]],inputOrderInvert[1]))
# 			for i in range(2,len(inputOrder)):
# 				#combines subsequent inputs
# 				if logicOperatorFlags[i-1]==0:
# 					value=fuzzyAnd(value,Inv(oldValue[inputOrder[i]],inputOrderInvert[i]))
# 				else:
# 					value=fuzzyOr(value,Inv(oldValue[inputOrder[i]],inputOrderInvert[i]))
# 		return value
# 	else:
# 		#returns savme value if now inputs
# 		return oldValue[currentNode]						
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
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb=returnSimParams()
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))

	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	
	
	sigmaNode=0
	sigmaNetwork=0
	individual=[0,1,0,0,0,0]
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	sigmaNode=0
	sigmaNetwork=0
	individual=[1,0,0,0,0,0]
	print(individual)
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	sigmaNode=.01
	sigmaNetwork=.01
	avg, stdev=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList, steps, h,p,hillOn,sigmaNode, sigmaNetwork)
	print(avg)
	print(stdev)
	
	
	
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='i')
	graph.add_edge('one','two', signal='i')	
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	ss['zero']=1
	ss['one']=1
	ss['two']=1
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList= setupGAparams(graph, ss)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[0,1,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))
	
	individual=[1,0,0,0,0,0]
	print(individual)
	print(runFuzzyModel(individual, async, individualParse,nodeList, initValues))


def testFuzzySimTimeCourse():

	csvfile=open('output2.csv', 'w')
	writer = csv.writer(csvfile)
    
	simCount=0	
	#two dicts for the models
	nodeUpdateDict={}
	ss={}
	sss=[]
	for i in range(1,5):
		ss[str(i)]=0
	ss['zero']=1
	ss['one']=0
	ss['two']=0
	ss['three']=0
	print(ss.keys())
	sss.append(ss)
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=constructor.parseKEGGdict('ko00001.keg', aliasDict, dict)
	graph.add_edge('zero','one', signal='a')
	graph.add_edge('zero','two', signal='a')
	graph.add_edge('one','two', signal='i')	
	graph.add_edge('two','three', signal='a')
	graph.add_edge('three','one', signal='a')	

	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes= setupGAparams(graph, sss)
	async, iters, simSteps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb=returnSimParams()
	print(individualLength)
	print(nodeList)
	print(graph.edges())
	individual=[0,0,0,0,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async,individualParse, nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

	print("with_noise .01")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)





	
	sigmaNode=0
	sigmaNetwork=0
	individual=[0,1,0,0,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)


	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,1,0,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,0,1,0,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	


	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,0,0,1,0,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	




	sigmaNode=0
	sigmaNetwork=0
	individual=[0,0,0,0,0,1,0,0,0,0]
	print(individual)
	for currentNode in range(0,len(nodeList)):
		print(writeFuzzyNode(currentNode,individual, individualParse, nodeList,nodeOrderList, nodeOrderInvertList, possibilityNumList, h,p,hillOn))
	
	print("noiseless")
	steps=simSteps
	sigmaNode=0
	sigmaNetwork=0
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	

	print("with_noise .01")
	sigmaNode=.001
	sigmaNetwork=.001
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)
	
	print("with_noise .001")
	sigmaNode=.01
	sigmaNetwork=.01
	sims=runFuzzyModel(individual, async, individualParse,nodeList, nodeOrderList, nodeOrderInvertList, possibilityNumList, initValueList[0], steps, h,p,hillOn,sigmaNode, sigmaNetwork, earlyEvalNodes)
	writer.writerows(sims)

#above are old testing functions, below are current testing functions

def runNoiseSimTest(graph, sss, networkNoises, nodeNoises, trials):
	async, iters, steps, generations, popSize, h, p, hillOn, bitFlipProb, crossoverProb,mutationProb,mu,lambd, complexPenalty, genSteps= returnSimParams()
	individualLength, individualParse, nodeList, nodeOrderList, initValueList, evaluateNodes, possibilityNumList, nodeOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss)
	#open file to save, use noises in name
	#save graph seperately
	toolbox, stats=buildToolbox(individualLength, bitFlipProb, len(sss))
	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss,individualLength,evaluateNodes, individualParse,nodeList,nodeOrderList,nodeOrderInvertList,possibilityNumList,initValueList, steps, h,p,hillOn,nodeNoise,networkNoise,popSize,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd, complexPenalty, genSteps)




def runFatimaGA():
	#very simple function to just run GA from fatima data with no frills. just see output from a single call to eaMuCommaLambda
	graph, ss = buildFatimaNetwork()
	iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb = returnSimParams()
	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList= setupGAparams(graph, ss)
	toolbox, hof, stats=buildToolbox(iters,individualLength, bitFlipProb, ss, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList,simSteps, )
	toolbox.register("evaluate", evaluate, iters=iters,ss=ss,  evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValueList=initValues, steps=simSteps, h=h,p=p,hillOn=hillOn)
	population=toolbox.population(n=popSize)
	algo.eaMuCommaLambda(population, toolbox, mu=100, lambda_=200, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)

def bootstrapRules(graph,sss,reSampleNum,sampleSize, counter, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList):
	ruleList=[]
	for i in range(0,reSampleNum):
		sss2=numpy.random.choice(sss, size=sampleSize, replace=True, p=None)
		ruleList.append(propBoolModel( individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList,graph,sss2))
	ruleDict=[]

	# for i in range(0,len(nodeList)):
	# 	ruleDict.append({})
	# for trialNum in range(0,len(ruleList)):
	# 	for nodeNum in range(0,len(ruleList[trialNum])):
	# 		ruleDict[nodeNum][ruleList[trialNum][nodeNum]]=0
	# for trialNum in range(0,len(ruleList)):
	# 	for nodeNum in range(0,len(ruleList[trialNum])):
	# 		ruleDict[nodeNum][ruleList[trialNum][nodeNum]]=ruleDict[nodeNum][ruleList[trialNum][nodeNum]]+1
	return ruleList

def drawGraph3(representative, filename, KEGGdict):
	rep=representative.copy()
	dictionary={}
	names=nx.get_node_attributes(representative, 'name')
	for n in rep.nodes():
		if len(names[n].split())==1:
			if names[n] in KEGGdict.keys():
				dictionary[n]=KEGGdict[names[n]]
				#print(rep.node[n]['name'])
		else :
			translated=''
			for word in names[n].split():
				word1=word.lstrip('ko:')
				word1=word1.lstrip('gl:')
				if word1 in KEGGdict.keys():
					translated=translated+KEGGdict[word1]+'-'
				else:
					translated=translated+word1+'-'
			dictionary[n]=translated
	repar= nx.relabel_nodes(rep,dictionary)
	#print(repar.nodes())
	#print(repar.edges())
	B=agraph.to_agraph(repar)        # convert to a graphviz graph\
	B.layout()            # neato layout
	B.draw(filename)       # write postscript in k5.ps with neato layout

def drawGraph2(representative, filename):
	B=agraph.to_agraph(representative)        # convert to a graphviz graph\
	B.layout()            # neato layout
	B.draw(filename)       # write postscript in k5.ps with neato layout

def readKEGG(lines, graph, KEGGdict):
	#the network contained in a KEGG file to the graph. 
	grouplist=[] #sometimes the network stores a group of things all activated by a single signal as a "group". We need to rewire the arrows to go to each individual component of the group so we save a list of these groups and a dictionary between id numbers and lists of elements that are part of that group
	groups={}
	i=0 #i controls our movement through the lines of code. 
	dict={} #identifies internal id numbers with names (much like the code earlier does for groups to id numbers of elements of the group)
	while i< len(lines):
		line=lines[i]
		# print line
		#add nodes
		if "<entry" in line:		
			nameline=line.split('name="')
			namestring=nameline[1]
			nameline=namestring.split('"')
			name=str.lstrip(nameline[0],'ko:')
			name=str.lstrip(name,'gl:')
			name=name.replace(' ','-')
			name=name.replace('ko:','')
			name=name.replace(':','-')
			
			# print(name)
			
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			#build a dictionary internal to this file
			idline=line.split('id="')
			idstring=idline[1]
			idline=idstring.split('"')
			id=idline[0]
			
			
			namelist=[]
			if '-' in name: 
				namelines=name.split('-')
				for x in namelines:
					if x in KEGGdict.keys():
						namelist.append(KEGGdict[x])
					else:
						namelist.append(x)
				name=namelist[0]+'-'
				for x in range(1,len(namelist)):
					name=name+namelist[x]+'-'
				name=name[:-1:]
			if name in KEGGdict.keys():
				dict[id]=KEGGdict[name]
			else:
				dict[id]=name
			#print(dict[id])
			#if the node is a group, we create the data structure to deconvolute the group later
			if type=='group':
				i=i+1
				newgroup=[]
				grouplist.append(id)
				while(not '</entry>' in lines[i]):
					if 'component' in lines[i]:
						newcompline=lines[i].split('id="')
						newcompstring=newcompline[1]
						newcompline=newcompstring.split('"')
						newcomp=newcompline[0]
						newgroup.append(newcomp)
					i=i+1
				groups[id]=newgroup
			else:
				graph.add_node(dict[id],{'name':dict[id],'type':type})
				#print(name)
				#print(id)
				i=i+1
		#add edges from KEGG file
		elif "<relation" in line:
			color='black'
			signal='a'
			subtypes=[]
			
			entryline=line.split('entry1="')
			entrystring=entryline[1]
			entryline=entrystring.split('"')
			entry1=entryline[0]
			
			entryline=line.split('entry2="')
			entrystring=entryline[1]
			entryline=entrystring.split('"')
			entry2=entryline[0]
			

			
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			while(not '</relation>' in lines[i]):
				if 'subtype' in lines[i]:
					nameline1=lines[i].split('name="')
					namestring1=nameline1[1]
					nameline1=namestring1.split('"')
					subtypes.append(nameline1[0])
				i=i+1
			
			
			#color and activation assignament based on the type of interaction in KEGG
			if ('activation' in subtypes) or ('expression' in subtypes):
				color='green'
				signal='a'
			elif 'inhibition' in subtypes:
				color='red'
				signal='i'
			elif ('binding/association' in subtypes) or('compound' in subtypes):
				color='purple'
				signal='a'
			elif 'phosphorylation' in subtypes:
				color='orange'
				signal='a'
			elif 'dephosphorylation' in subtypes:
				color='pink'
				signal='i'
			elif 'indirect effect' in subtypes:
				color='cyan'
				signal='a'
			elif 'dissociation' in subtypes:
				color='yellow'
				signal='i'
			elif 'ubiquitination' in subtypes:
				color='cyan'
				signal='i'
			else:
				print('color not detected. Signal assigned to activation arbitrarily')
				print(subtypes)
				signal='a'
			
			#this code deconvolutes the "groups" that KEGG creates 
			#each element of the group gets an edge to or from it instead of the group
			if entry1 in grouplist:
				for current1 in groups[entry1]:
					node1=dict[current1]
					if entry2 in grouplist:
						for current2 in groups[entry2]:
							node2=dict[current2]
							graph.add_edge(node1, node2, color=color, subtype=name, type=type, signal=signal )
					else:
						node2=dict[entry2]
						graph.add_edge(node1, node2, color=color,  subtype=name, type=type, signal=signal )
			elif entry2 in grouplist:
				node1=dict[entry1]
				for current in groups[entry2]:
					if current in grouplist:
						for current1 in groups[current]:
							if current in dict.keys():
								node2=dict[current]
								graph.add_edge(node1,node2, color=color,  subtype=name, type=type, signal=signal )
							else: 
								print('groups on groups on groups')
						#print('groups on groups')
					elif current in dict.keys():
						node2=dict[current]
						graph.add_edge(node1,node2, color=color, subtype=name, type=type, signal=signal )
					else:
						print('not in the keys')
						print(current)
			else:
				node1=dict[entry1]
				node2=dict[entry2]
				
				graph.add_edge(node1,node2, color=color, subtype=name, type=type, signal=signal )
		i=i+1
def read_biopax(lines, graph):
	soup = BeautifulSoup(''.join(lines), 'xml')

	##old code, potentially useless
	# pathways = soup.find_all('pathway')
	# for pathway in pathways:
	# 	name = pathway.find('NAME').string
	# 	synonyms = [s.string for s in pathway.find_all('SYNONYMS')]
	# 	component_ids = [c.get('rdf:resource')[1:] for c in pathway.find_all('PATHWAY-COMPONENTS')]
	# 	print name
	# 	for c_id in component_ids:
	# 		component = soup.find(attrs={'rdf:ID': c_id})
	# 		print component.prettify()
	##end old code


	#treating edge classes as nodes for now, will simplify later
=======
	#constuct a BeautifulSoup4 parser object with the libxml backend
	#input in lines is given as an array, the join function converts to 1 long string
	#On windows this fails, try 'lxml' instead of 'lxml-xml'
	soup = BeautifulSoup(''.join(lines), 'xml')

	#in BioPAX, both physical objects (proteins, genes, etc.) and interactions (activation, 
	#inhibition, etc.) are represented as NODES.  The edges in the system represent membership 
	#or participation in these interactions.  Accordingly, we parse all objects into nodes below.

	#for every biopax class name (defined on lines 29-31)

## Rohith's OLD version, kept for history
# def parseKEGGdict(filename):
# 	#makes a dictionary to convert ko numbers from KEGG into real gene names
# 	#this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two.
# 	dict={}
# 	inputfile = open(filename, 'r')
# 	lines = inputfile.readlines()
# 	for line in lines:
# 		if line[0]=='D':
# 			kline=line.split('      ')
# 			kstring=kline[1]
# 			kline=kstring.split('  ')
# 			k=kline[0]
# 			nameline=line.replace('D      ', 'D')
# 			nameline=nameline.split('  ')
# 			namestring=nameline[1]
# 			nameline=namestring.split(';')
# 			name=nameline[0]
# 			dict[k]=name
# 	return dict

def runFatimaSim(networkNoises, nodeNoises, trials):
	#simulate data on Fatima network. Output how well it goes. 
	samples=20
	graph, sss=buildFatimaNetwork()
	sss1=synthesizeInputs(graph,samples)
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb ,mutationProb, mu,lambd= returnSimParams()
	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes=setupGAparams(graph, sss1)
	print("individualLength")
	print(individualLength)
	#open file to save, use noises in name
	#save graph seperately
	
	toolbox, stats=buildToolbox(individualLength, bitFlipProb, samples)

	for networkNoise in networkNoises:
		print('networkNoise')
		print(networkNoise)
		for nodeNoise in nodeNoises:
			print('nodeNoise')
			print(nodeNoise)
			fuzzySimulator(async, iters, toolbox, stats, graph, trials, sss1,individualLength,evaluateNodes, individualParse,nodeList,inputOrderList,inputOrderInvertList,possibilityNumList,initValueList,simSteps ,nodeNoise,networkNoise,popSize ,crossoverProb,mutationProb,generations,earlyEvalNodes, mu,lambd)
		
		

if __name__ == '__main__':
	# ss={}
	#data=loadFpkms('C:/Users/Rohith/Desktop/ruleMaker_GA/Data/fatima_fpkms.csv')
	#sortFpkms(data,ss)
	# print(ss)
	#teseter.LiuNetwork2SimTest()
	#runFatimaSim([.01,.05],[.01,.02],1)
	#runFatimaSim([0],[0.001],10)
	#tester.LiuNetwork2SimTest(20,1)
	# inputNum=50
	# counters=[]
	# graphStart=tester.LiuNetwork2Builder()
	# sss=synthesizeInputs(graphStart,inputNum)
	# graph=simplifyNetwork(graphStart,sss[0])
	# nx.write_graphml(graph,'Liu2.graphml')
	# for i in range(0,1):
	# 	sss=synthesizeInputs(graph,inputNum)
	# 	counters.append(testBootstrap(5,10,graph,sss))
	# print(counters)

def setupBootstrapParams(graph):
	repeat=True
	while(repeat):
		repeat=False
		for node in graph.nodes():
			if node in graph.successors(node):
				graph.remove_edge(node,node)
				repeat=True
	
	evaluateNodes=[] 
	individualParse=[] 
	inputOrderList=[] 
	inputOrderInvertList=[]
	possibilityNumList=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 100
	counter=0
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
		inputOrderList.append(possibilities)
		inputOrderInvertList.append(activities)
		possibilityNumList.append(len(possibilities))
		if len(possibilities)==0:
			logNum=0
		else:
			logNum=ceil(log(len(possibilities))/log(2))
		individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(possibilities)-1)])
		counter=counter+logNum+len(possibilities)
		# print(counter)
		# print(individualParse)
	return counter, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList

def testBootstrap(reSampleNum,sampleSize, graph,sss):

	# graph=nx.DiGraph()
	# graph.add_edge('1','3',signal='a')
	# graph.add_edge('2','3',signal='a')

	individualLength, individualParse, nodeList, inputOrderList, initValueList, evaluateNodes, possibilityNumList, inputOrderInvertList, earlyEvalNodes= setupGAparams(graph, sss)
	nodeNoise=0
	networkNoise=0
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps=returnSimParams()
	#for i in range(0,trials):
	individual=genRandBits(individualLength)
	# print(individualLength)
	# print(individual)
	
	truths=[]
	#print(writeFuzzyModel(individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))
	for i in range(0,len(nodeList)):
		truths.append(writeFuzzyNode(i,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))		
	newSSS=[]
	for k in range(0,len(sss)):
		if async:
			boolValues=iterateFuzzyModel(async, iters,individual, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps,nodeNoise,networkNoise,earlyEvalNodes)
		else:
			boolValues=runFuzzyModel(individual, async,individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList[k], genSteps ,nodeNoise,networkNoise,earlyEvalNodes)
		newSS=copy.deepcopy(sss[k])
		for j in range(0,len(evaluateNodes)):
			newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
		newSSS.append(newSS)
		#update eval nodes
		#print(evaluate(individual, async, iters, newSSS, evaluateNodes, individualParse, nodeList, inputOrderList, inputOrderInvertList, possibilityNumList, initValueList, steps, ,earlyEvalNodes,complexPenalty))
	
	initValueList=[]
	for j in range(0,len(newSSS)):
		initValueList.append([])
	for i in range(0,len(nodeList)):
		for j in range(0,len(newSSS)):
			ss=newSSS[j]
			if  nodeList[i] in newSSS[0].keys():
				initValueList[j].append(ss[nodeList[i]])
			else:
				initValueList[j].append(0.5)
				
	# print(individual)
	for i in range(0,len(nodeList)):
		# print('new Node')
		# print(individual)
		# print(individualParse[i])
		deviation=0
		# for steadyStateNum in range(0,len(newSSS)):
		# 	derivedVal=propUpdate(i,initValueList[steadyStateNum],individual[individualParse[i][0]:individualParse[i][2]], individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList, )
		# 	deviation=deviation+(derivedVal-newSSS[steadyStateNum][nodeList[i]])**2
		# print(deviation)

	ruleList=bootstrapRules(graph,newSSS,reSampleNum,sampleSize, individualLength, individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList)
	wrongCount=0
	correctCount=0
	for item in ruleList:
		for i in range(0,len(nodeList)):
			if possibilityNumList[i]>1:
				if truths[i]==writeFuzzyNode(i,item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
					correctCount=correctCount+1
				else:
					wrongCount=wrongCount+1
					print(truths[i])
					print(writeFuzzyNode(i,item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))

	return wrongCount, correctCount
		# print(item)
		# print(writeFuzzyModel(item, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList))


def hillInv(x, inverter,h,p,hillOn): #inverts if necessary then applies hill fun
	if inverter:
		return hill(1-x, h,p,hillOn)
	else:
		return hill(x, h,p,hillOn)
		
def returnSimParams():
	simSteps=100 # number of steps each individual is run when evaluating
	generations=50 # generations to run
	popSize=100 #size of population
	mu= 100#individuals selected
	lambd= 200#children produced
	bitFlipProb=.1 # prob of flipping bits inside mutation
	crossoverProb=.2 # prob of crossing over a particular parent
	mutationProb=.2 # prob of mutating a particular parent
	async=False # run in asynchronous mode
	iters=100 #number of simulations to try in asynchronous mode
	complexPenalty=False #penalize models with increased complexity
	genSteps=10000 # steps to find steady state with fake data
	return async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps


def loadFatimaData(filename,tempDict): #reads data from differential expression file
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]!='#':
			kline=line.split(' ')
			tempDict[str.upper(kline[0][3:])]=logistic.cdf(float(kline[1]))


def autoSimulator(trials, toolbox, stats,  sss,params, model, simulator):
	hofSize=10
	avgSSEs=[]
	hofs=[]
	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	truthCounter=[0 for number in xrange(hofSize)]
	trueNodeList=[]
	for i in range(0,trials):
		individual=genRandBits(individualLength)
		newSSS=[]
		for k in range(0,len(sss)):
			if async:
				individual, params, model, simulator, initValues, iters
				boolValues=averageResultModelSim(individual, params, model, simulator, initValueList[k], iters)
			else:
				boolValues=runFuzzyModel(individual, params, model, simulator, initValueList[k], 0)
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(evaluateNodes)):
				newSS[nodeList[evaluateNodes[j]]]=boolValues[evaluateNodes[j]]
			newSSS.append(newSS)
		#update eval nodes
		print(evaluate(individual, params, model, simulator, sss))
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		#print(len(sss))
		toolbox.register("evaluate", evaluate, params=params, model=model, simulator=simulator, sss=newSSS)
		#toolbox.register("evaluate", evaluate, ss=ss, evaluateNodes=evaluateNodes,individualParse= individualParse, nodeList=nodeList, inputOrderList=inputOrderList, inputOrderInvertList=inputOrderInvertList, possibilityNumList=possibilityNumList, initValues=boolValues, steps=steps, h=h,p=p,hillOn=hillOn)
		population=toolbox.population(n=popSize)
		hof = tools.HallOfFame(hofSize, similar=numpy.array_equal)
		output = tools.Logbook()
		output=algo.eaMuCommaLambda(population, toolbox, mu=mu, lambda_=lambd, stats=stats, cxpb=crossoverProb, mutpb=mutationProb, ngen=generations, verbose=True, halloffame=hof)
		truth=writeModel(individual, model)
		for j in range(0,10):
			bestRun=(writeModel(hof[j], model))
			if truth==bestRun:
				truthCounter[j]+=1
				break
			
			elif j==0:
				truthlines=truth.split('\n')
				newlines=bestRun.split('\n')
				trueNodes=0
				for k in range(0,len(truthlines)):
					if truthlines[k]==newlines[k]:
						trueNodes=trueNodes+1
					else:
						print("incorrect pair: true then test")
						print(truthlines[k])
						print(newlines[k])
				trueNodeList.append(trueNodes)
		avgs=[output[1][k]['min'] for k in range(0,len(output[1]))]
		#print(avgs)
		hofs.append(hof)
		temp=[]
		for hofind in range(0,len(hof)):
			tempVal=0
			for bit in range(0,len(hof[hofind])):
				if not hof[hofind][bit]==individual[bit]:
					tempVal=tempVal+1
			temp.append(tempVal)
		hofScores.append(temp)
		truthIndividuals.append(individual)
		truthAvgs.append(boolValues)
		gc.collect()
	f = open('hof_differences_'+str(nodeNoise)+str(networkNoise)+'.txt', 'w')
	g = open('hof_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	h = open('truth_individuals'+str(nodeNoise)+str(networkNoise)+'.txt',  'w')
	f.write(str(hofScores))
	f.write('\n')
	for hofy in range(0,len(hofs)):
		g.write(str(hofs[hofy]))
		g.write('\n')
	h.write(str(truthIndividuals))
	h.write('\n')
	f.close()
	g.close()
	h.close()
	print('# true of # trials')
	for k in range(0,len(truthCounter)):
		print(truthCounter[k])
	print(trials)
	print("for the incorrects, by node")
	print(trueNodeList)
	print(len(truthlines))