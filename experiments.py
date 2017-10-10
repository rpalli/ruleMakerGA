def deltaTester(ultimate, i, model,  newSSS, params, minny):
	print(i)
	modifyFlag=False
	copied=list(ultimate)
	copied[i]=1-copied[i]
	for node in range(0,len(model.nodeList)):
		better=True
		#get start and end indices for node in individual
		if node==(len(model.nodeList)-1):
			end=len(ultimate)-1
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		if model.andLenList[node]>0:
			if sum(copied[start:end])>0:
				newtot=numpy.sum(evaluateByNode(copied, params.cells, model,  newSSS, params, KOlist, KIlist))
				if newtot<minny:
					modifyFlag=True
	return modifyFlag

def simTester(graph, name):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	# set up params, sss, model, samples


	params=sim.paramClass()
	# rewire graph if necessary
	if params.rewire:
		graph=nc.rewireNetwork(graph)

	# # remove two node bunches
	# i=0
	# while i < range(len(graph.nodes())):
	# 	nodes=graph.nodes()
	# 	node=nodes[i]
	# 	if len(graph.predecessors(node))==0:
	# 		successorless=True
	# 		for successor in graph.successors(node):
	# 			if graph.successors(successor)>0 or graph.predecessors(successor)>1:
	# 				successorless=False
	# 		if successorless:
	# 			graph.remove_node(node)
	# 			i=i-1
	# 	i+=1

	# # remove orphan nodes
	# for i in range(len(graph.nodes())):
	# 	nodes=graph.nodes()
	# 	node=nodes[i]
	# 	if len(graph.predecessors(node))==0 and len(graph.successors(node))==0:
	# 		graph.remove_node(node)



	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	samples=params.samples

	# set up necessary lists for output
	truthList= [] 
	testList =[]
	devList=[]


	knockoutLists=[]
	knockinLists=[]

	for q in range(samples):
		temRand=rand.randint(0,len(model.nodeList))
		knockoutLists.append([])
		knockinLists.append([])


	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,params.trials):
		#generate random set of logic rules to start with
		individual=genBits(model)
		for node in range(0,len(model.nodeList)):
			if model.andLenList[node]>0:
				if node==len(model.nodeList)-1:
					end=len(model.nodeList)
				else:
					end=model.individualParse[node+1]
				if numpy.sum(individual[model.individualParse[node]:end])==0:
					individual[model.individualParse[node]]=1

		# generate Boolean model for this trial
		output=runProbabilityBooleanSims(individual, model, samples, params.cells, params, knockoutLists, knockinLists)
		# copy new output into newSSS and initial values
		newSSS=[]
		for k in range(0,samples):
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(model.nodeList)):
				newSS[model.nodeList[j]]=output[k][j]
			newSSS.append(newSS)
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(model.nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  model.nodeList[j] in sss[0].keys():
					newInitValueList[k].append(ss[model.nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		model.initValueList=newInitValueList
		#perform GA run
		dev, bruteOut=GAsearchModel(model, newSSS, params, knockoutLists, knockinLists)
		# add output of this trial to lists
		truthList.append(individual)
		testList.append(list(bruteOut))
		devList.append(dev)
	# set up output and save as a pickle
	outputList=[truthList,testList, devList,model.size, model.nodeList, model.individualParse,model.andNodeList ,model.andNodeInvertList, model.andLenList,	model.nodeList, model.nodeDict, model.initValueList]
	pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )

def partitionTester(graph,name):
	# first trim the graph to just what is in the dataset and the dataset to just what is in the graph

	if len(graph.nodes())>1:	
		# dev1, bruteOut1, model1= partitionTest(graph, False, False) # complicated method with max included
		# dev2, bruteOut2, model2= partitionTest(graph, True, False) # complicated method after throwing out max
		dev3, bruteOut3, model3= partitionTest(graph, False, True) # divide by max
		dev4, bruteOut4, model4= partitionTest(graph, True, True) # divide by second highest
		# set up output and save as a pickle
		# outputList=[[dev1,dev2,dev3,dev4],[bruteOut1,bruteOut2,bruteOut3, bruteOut4],[model1,model2,model3,model4]]
		outputList=[[dev3,dev4],[bruteOut3,bruteOut4],[model3,model4]]

		pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )
	else:
		print('not enough overlap')

def partitionTest(graph, maxRem, divider):
	
	# first construct the input -omics dataset
	if maxRem:
		if divider:
			fileName='data_T_T.pickle'
		else:
			fileName='data_T_F.pickle'
	else:
		if divider:
			fileName='data_F_T.pickle'
		else:
			fileName='data_F_F.pickle'
	ssDict=pickle.Unpickler(open( fileName, "rb" )).load()
	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, maxRem, divider)
	params=sim.paramClass()
	
	graph=nc.simplifyNetwork(graph, ssDict) # simplify the graph
	graph=mc.cutMotifs(graph)
	print(graph.nodes())
	keyList=ssDict.keys()
	# generate SSS which is a list of SSs where SS is a dict pointing to some values
	sss=[{} for i in range(len(ssDict[keyList[0]]))]
	newInitValueList=[[] for i in range(len(ssDict[keyList[0]]))]
	for i in range(len(sss)):
		for key in keyList:
			if key in graph.nodes():
				sss[i][key]=ssDict[key][i]
	
	model=sim.modelClass(graph,sss)
	# generate empty knockout and knockin lists for support
	knockoutLists=[]
	knockinLists=[]
	for q in range(len(sss)):
		knockoutLists.append([])
		knockinLists.append([])

	# generate a set of initial values which are just the values from the SS in order. 
	newInitValueList=[]
	for j in range(0,len(sss)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(sss)):
			ss=sss[k]
			if  model.nodeList[j] in sss[0].keys():
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	print('setup successful')
	#perform GA run
	dev, bruteOut= GAsearchModel(model, sss, params, knockoutLists, knockinLists)
	return dev, bruteOut, model

def indDistFind(name):
	ssDict=pickle.Unpickler(open( 'data_F_T.pickle', "rb" )).load()
	[[dev1,dev3],[bruteOut1,bruteOut3],[model1,model3]]=pickle.Unpickler(open( 'pickles/'+name+"_output.pickle", "rb" )).load()
	# generate SSS which is a list of SSs where SS is a dict pointing to some values
	keyList=ssDict.keys()
	sss=[{} for i in range(len(ssDict[keyList[0]]))]
	newInitValueList=[[] for i in range(len(ssDict[keyList[0]]))]
	for i in range(len(sss)):
		for key in keyList:
			if key in graph.nodes():
				sss[i][key]=ssDict[key][i]
	params=sim.paramClass()

	simulator=sim.simulatorClass('bool')
	importanceScores=sim.calcImportance(bruteOut1,params,model1,simulator, sss)

	pickle.dump( importanceScores, open( name+"_scores.pickle", "wb" ) )

def findClustDiff(lister):
	mixer= mix.GaussianMixture(n_components=2, covariance_type='full').fit([[a] for a in lister])
	return abs(mixer.means_[0]-mixer.means_[1])

def testDiscretizationSetup():

	findPathways()

	geneDict=mr.readData() # get -omics data

	ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, True, True)
	pickle.dump( ssDict, open( "data_T_T.pickle", "wb" ) )


	ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, False, True)
	pickle.dump( ssDict, open( "data_F_T.pickle", "wb" ) )

	ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, True, False)
	pickle.dump( ssDict, open( "data_T_F.pickle", "wb" ) )


	ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, False, False)
	pickle.dump( ssDict, open( "data_F_F.pickle", "wb" ) )

if __name__ == '__main__':
	import time
	start_time = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	results = parser.parse_args()
	graphName=results.graph
	name=graphName[:-8]+'_1'

	graph = nx.read_gpickle(graphName)
	indDistFind(name)


	# testDiscretizationSetup()


	# for numgraph in ['04066','04350','04380','04612','04630','04650','04657','04659','04660']:
	# valueList=[]
	# for key in geneDict.keys():
	# 	print(key)
	# 	valueList.append(findClustDiff(geneDict[key]))
	
	# pickle.dump( valueList, open( "distance_data.pickle", "wb" ) )
	# valueList=pickle.Unpickler(open( "distance_data.pickle", "rb" )).load()
	# print(len(valueList))
	# valueList=[valueList[i] for i in range(0,len(valueList),50)]
	# print(len(valueList))
	# bins = numpy.linspace(0,6000, 500)
	# n, bins, patches = plt.hist(valueList, 20, range=(0,20),  normed=0, facecolor='green')

	# plt.xlabel('distance between clusters')
	# plt.ylabel('Probability Density')
	# plt.title('Histogram of cluster distances')
	# # plt.axis([40, 160, 0, 0.03])
	# plt.grid(True)

	# plt.savefig('histogram_cluster_dist_1.png', bbox_inches='tight')




	# for numgraph in ['04657','04659','04660']:
	# 	graphName='hsa'+numgraph+'_unsimplified.gpickle'
	# 	name=graphName[:-8]+'_1'
	# 	graph = nx.read_gpickle(graphName)
	# 	partitionTester(graph, name, geneDict)
	# 	print("--- %s seconds ---" % (time.time() - start_time))