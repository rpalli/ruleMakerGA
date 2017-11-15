#import python packages
import networkx as nx
import operator
from sets import Set
import scipy.stats as stat
import matplotlib.pyplot as plt
import requests

# import other pieces of our software
import networkConstructor as nc
import murphyReader as mr

def read_gmt(filename):
	gmt_dict={}
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		newline=line.split('\t')
		gmt_dict[newline[0]]=Set(newline[2:])
	return gmt_dict

def find_overlaps(filename,geneDict):
	overlapsets=[]
	genes=Set(geneDict.keys())
	keggDict=read_gmt(filename)
	for key in keggDict.keys():
		if len(genes.intersection(keggDict[key]))>4:
			overlapsets.append(key)
			print(key)
			print(len(genes.intersection(keggDict[key])))
	return overlapsets

def retrieveGraph(name,aliasDict,dict1,dict2, geneDict):
	print(name)
	namelist=name.split('_')
	namelist.pop(0)
	requester='http://rest.kegg.jp/find/pathway/'+namelist.pop(0)
	for item in namelist:
		requester=requester+'+'+item

	# print(requester)

	r=requests.get(requester)
	genes=Set(geneDict.keys())

	lines=r.text
	print(lines)
	if len(lines.split('\n')[0].split(':'))>1:

		code=lines.split('\n')[0].split(':')[1][3:8] # KEGG number of overlapped pathway
		graph=nx.DiGraph()
		coder=str('ko'+code) 
		nc.uploadKEGGcodes([coder], graph, dict2)
		coder=str('hsa'+code)
		nc.uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		if len(list(nx.connected_component_subgraphs(graph.to_undirected() )))>0:
			newgraph = max(nx.connected_component_subgraphs(graph.to_undirected()), key=len)
			newOverlap=genes.intersection(set(newgraph.nodes()))
			print(len(newOverlap))
			graphLen= nx.average_shortest_path_length(newgraph, weight=None)
			graph=nc.simplifyNetwork(graph, geneDict)
			print('nodes, edges')
			print(len(graph.nodes()))
			print(len(graph.edges()))
			if len(newOverlap)>4:
				nx.write_graphml(graph,coder+'.graphml')
				nx.write_gpickle(graph,coder+'.gpickle')
				return coder+'.gpickle'
		else:
			graphLen=0
		print('graphlen:')
		print(graphLen)
		#check is average path length is sufficiently long


		
	else:
		print('not found:')
		print(requester)
		print(lines)
	return ''

def findPathways():
	geneDict=mr.readData() # get -omics data
	
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	namelist=find_overlaps('filtered.c2.cp.kegg.v3.0.symbols.gmt',geneDict)
	print('num of overlap nodes')
	print(len(namelist))
	for name in namelist:
		retrieveGraph(name,aliasDict,dict1,dict2, geneDict)
	# # read in list of codes then load them into network
	# #inputfile = open('ID_filtered.c2.cp.kegg.v3.0.symbols.txt', 'r')
	# for code in lines:
		
	# 	# if(len(graph.edges())>1):
	# 	# 	graph=simplifyNetwork(graph, data)
		
	# 	#graph = utils.simpleNetBuild()
	# 	#coder='unrolled'

	# 	#graph=utils.LiuNetwork1Builder()
	# 	# coder='liu'

	# 	# if the code has an interesting logic rule to find, run the algorithm... 
	# 	checker=False 
	# 	for x in graph.in_degree():
	# 		if x>1:
	# 			checker=True
	# 	if(checker):
	# 		print(coder)
	# 		codelist.append(coder)			
	# 		nx.write_graphml(graph,coder+'_unsimplified.graphml')
	# 		nx.write_gpickle(graph,coder+'_unsimplified.gpickle')

def simplifyNetworkpathwayAnalysis(graph, ss):
	#network simplification algorithm. 
	# # 1. remove nodes with no input data
	# # 2. remove edges to nodes from complexes they are a part of 
	# # 3. save graph for future use with influence score calculation
	# # 4. remove nodes with only 1 input
	# # 5. remove nodes with only 1 output
	# # 6. remove self edges

	# 1. remove nodes with no input data
	removeNodeList= [x for x in graph.nodes() if  not x  in ss.keys()]
	for rm in removeNodeList:
		for start in graph.predecessors(rm):
			for finish in graph.successors(rm):
				edge1=graph.get_edge_data(start,rm)['signal']
				edge2=graph.get_edge_data(rm,finish)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(start,finish,signal='i')
				else:
					graph.add_edge(start,finish,signal='a')
		graph.remove_node(rm)
	#print(graph.nodes())
	flag=True

	# 2. remove dependence of nodes on complexes that include that node
	for node in graph.nodes():
		predlist=graph.predecessors(node)
		for pred in predlist:
			if '-' in pred:
				# # print(pred)
				genes=pred.split('-')
				flag=True
				for gene in genes:
					if not gene in predlist:
						flag=False
				if flag:
					graph.remove_edge(pred,node)
	flag=True
		
	# 3. save graph here for use in constructing the model to measure influence scores


	# 4. rewire nodes that have only one upstream node
	#print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1) ]
	for rm in removeNodeList:
		if ss[rm] < ss[graph.predecessors(rm)[0]]:
			for after in graph.successors(rm):
				before=graph.predecessors(rm)[0]
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
		else:
			for before in graph.predecessors(graph.predecessors[rm][0]):
				newrm=graph.predecessors[rm][0]
				edge1=graph.get_edge_data(before,newrm)['signal']
				edge2=graph.get_edge_data(newrm,rm)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(before,rm,signal='i')
				else:
					graph.add_edge(before,rm,signal='a')
			graph.remove_node(newrm)
	flag=True

	# 5. rewire nodes that have only one downstream node
	removeNodeList= [x for x in graph.nodes() if (len(graph.successors(x))==1) ]
	for rm in removeNodeList:
		if len(graph.successors(x))==1:
			for start in graph.predecessors(rm):
				finish=graph.successors(rm)[0]
				edge1=graph.get_edge_data(start,rm)['signal']
				edge2=graph.get_edge_data(rm,finish)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(start,finish,signal='i')
				else:
					graph.add_edge(start,finish,signal='a')
			graph.remove_node(rm)
	flag=True

	# 6. remove self edges
	for edge in graph.edges():
		if edge[0]==edge[1]:
			graph.remove_edge(edge[0],edge[1])
	return graph

def findClustDiff(lister):
	mixer= mix.GaussianMixture(n_components=2, covariance_type='full').fit([[a] for a in lister])
	return abs(mixer.means_[0]-mixer.means_[1])

def indDistFind(name, ssDict):
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

#calculate average induced change by node
def calcImportance(individual,params,model,simulator, sss):
	importanceScores=[]
	knockoutLists=[]
	knockinLists=[]
	for q in range(len(sss)):
		knockoutLists.append([])
		knockinLists.append([])
	for node in range(len(model.nodeList)):
		tempImp=0
		SSEs=[]
		for j in range(0,len(sss)):
			ss=sss[j]
			initValues=model.initValueList[j]
			if initValues[node]>.1:
				initValues[node]=initValues[node]-.1
			else:
				initValues[node]=initValues[node]+.1
			SSE1=0
			boolValues =runModel(individual, model,simulator, model.initValueList[j], params, knockoutLists, knockinLists)
			for i in range(0, len(model.nodeList)):
				SSE1+=(boolValues[i]-ss[model.nodeList[i]])**2
				initValues=model.initValueList[j]
			if initValues[node]>.9:
				SSE2=SSE1
			else:
				initValues[node]=initValues[node]+.1
				SSE2=0
				boolValues=runModel(individual, model,simulator, model.initValueList[j], params, knockoutLists, knockinLists)
				for i in range(0, len(model.nodeList)):
					SSE2+=(boolValues[i]-ss[model.nodeList[i]])**2
					initValues=model.initValueList[j]
			SSEs.append(SSE1+SSE2)
		importanceScores.append(sum(SSEs))
	return importanceScores

def calcZscore(RAvals, importanceScores, RAnodeVals, bootstraps):
	randomsets=len(RAvals)
	resultset=[]
	for j in range(bootstraps):
		resultset.append(numpy.sum([score*RAvals[randomsets*random()] for score in importanceScores]))
	stdv= numpy.stdev(resultset)
	avg = numpy.mean(resultset)
	summer=0.
	for i in range(len(RAnodeVals)):
		summer+=RAnodeVals[i]*importanceScores[i]
	score= (summer-avg)/stdv
	return zscore

def calcPathwayScore():

	scoreDict={}
	calcImportance(individual,params,model,simulator, sss)
	for group in groupDict:
		scoreDict[group]=calcZscore(groupDict[group][0], importanceScores, groupDict[group][1], params.bootstraps)



