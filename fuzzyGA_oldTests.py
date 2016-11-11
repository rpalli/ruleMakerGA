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
