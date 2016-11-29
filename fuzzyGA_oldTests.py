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
		
def hillInv(x, inverter,h,p,hillOn): #inverts if necessary then applies hill fun
	if inverter:
		return hill(1-x, h,p,hillOn)
	else:
		return hill(x, h,p,hillOn)

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
