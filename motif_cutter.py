# import external files
import networkx as nx
import math
import networkx.algorithms.isomorphism as iso
import networkx.readwrite
from networkx.algorithms import isomorphism
import pickle
import os
# import ruleMakerGA files
# import networkConstructor as nc
# import utilities as utils

# print out the motif and cluster numbers in a neat format
def printMotifsClusters(motifs, clusters):
	for i in clusters.keys():
		print('cluster '+str(i)+' has '+str(len(clusters[i]))+' elements')
	for i in motifs.keys():
		print('motif '+str(i)+' has '+str(len(motifs[i]))+' elements')

# generate a hash table to use to search 3 node motifs
def fillHashTable(num, dict):
	picklefile=open('motifs/clusters.pkl','rb')
	clusters=pickle.load(picklefile)
	picklefile.close()
	clustering=[]
	motifs=[]
	for i in range(0,num):
		newnum=[]
		newmotif=[]
		clustering.append(newnum)
		motifs.append(newmotif)
		graph = nx.read_gpickle('motifs/pickles/model'+str(i)+'.pkl')
		
		#code here adds signal property to each edge. this property is either a for activation or i for inhibition
		#we continue this notation when we load from the KEGG files and then use that information to compare the edges
		for currentEdge in graph.edges():
			coloring= graph.get_edge_data(currentEdge[0],currentEdge[1])['color']
			if coloring=='green':
				graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
			elif coloring=='blue':
				graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
			elif coloring=='signal' or coloring=='Signal':
				graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
			else:
				graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'i'
		hashvalue= hashGraph(graph)
		#print(graph.edges())
		dict[tuple(hashvalue)].append([graph,clusters[i],i])
	# print(clusters)
	return clustering, motifs

# generate a hash table to use to search 3 node motifs
def fillHashTable3node(dict):
	print('filling hash table')
	folderlist=['3_node']

	classification={}
	motifs={}
	for folder in folderlist:
		classificationDict=[8,8, 7, 9, 8, 2, 1, 1, 7, 5, 2, 4, 1, 1, 6, 9, 2, 3, 1, 1, 7, 9 ,3, 3]
		filelist= [ f for f in os.listdir('motif_library/'+folder+'/pickles') ]
		for i in range(0,len(filelist)):
			graph = nx.read_gpickle('motif_library/'+folder+'/pickles/model'+str(i)+'.pkl')
			#code here adds signal property to each edge. this property is either a for activation or i for inhibition
			#we continue this notation when we load from the KEGG files and then use that information to compare the edges

			graph.remove_node('Signal')
			for currentEdge in graph.edges():
				coloring= graph.get_edge_data(currentEdge[0],currentEdge[1])['color']
				if coloring=='green':
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
				elif coloring=='blue':
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
				elif coloring=='signal' or coloring=='Signal':
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
				else:
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'i'
			hashvalue= hashGraph(graph)
			#print(graph.edges())
			#print(i)
			if not tuple(hashvalue) in dict.keys():
				dict[tuple(hashvalue)]=[]
			# if classificationDict[i]=='short':
				# drawGraph3(graph,'short-'+str(i)+'.bmp',KEGGdict)
			dict[tuple(hashvalue)].append([graph,classificationDict[i],folder,i])
			newnum=[]
			classification[classificationDict[i]]=newnum
			newmotif=[]
			motifs[folder,i]=newmotif
		# print(folder)
	# for element in dict.keys():
		# for graph1 in dict[element]:
			# print(graph1[0].edges())
		# print(classification)
	return classification, motifs

# generate a hash table to use to search 4 node motifs
def fillHashTable5node(dict):
	print('filling hash table')
	folderlist=['6','7_degree4','7_maxdegree3','8']

	classification={}
	motifs={}
	for folder in folderlist:
		picklefile=open('motif_library/'+folder+'/classificationDict.pkl','rb')
		classificationDict=pickle.load(picklefile)
		picklefile.close()
		filelist= [ f for f in os.listdir('motif_library/'+folder+'/FF/pickles') ]
		for i in range(0,len(filelist)):
			graph = nx.read_gpickle('motif_library/'+folder+'/FF/pickles/model'+str(i)+'.pkl')
			#code here adds signal property to each edge. this property is either a for activation or i for inhibition
			#we continue this notation when we load from the KEGG files and then use that information to compare the edges

			graph.remove_node('Signal')
			for currentEdge in graph.edges():
				coloring= graph.get_edge_data(currentEdge[0],currentEdge[1])['color']
				if coloring=='green':
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
				elif coloring=='blue':
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
				elif coloring=='signal' or coloring=='Signal':
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'a'
				else:
					graph.edge[currentEdge[0]][currentEdge[1]]['signal'] = 'i'
			hashvalue= hashGraph(graph)
			#print(graph.edges())
			#print(i)
			if not tuple(hashvalue) in dict.keys():
				dict[tuple(hashvalue)]=[]
			# if classificationDict[i]=='short':
				# drawGraph3(graph,'short-'+str(i)+'.bmp',KEGGdict)
			dict[tuple(hashvalue)].append([graph,classificationDict[i],folder,i])
			newnum=[]
			classification[classificationDict[i]]=newnum
			newmotif=[]
			motifs[folder,i]=newmotif
		print(folder)
	# for element in dict.keys():
		# for graph1 in dict[element]:
			# print(graph1[0].edges())
		print(classification)
	return classification, motifs
# build an empty hash table to be filled in
def buildEmptyHashTable():
	dict={}
	# print('done building table')
	return dict
# hash function to speed search
def hashGraph(graph):
	hashvalue=[]
	nodelist=graph.nodes()
	ins=graph.in_degree()
	outs=graph.out_degree()
	#print(ins)
	for i in range(0,len(ins)):
		hashvalue.append(0)
		hashvalue.append(0)
	for node in nodelist:	
		hashvalue[ins[node]]=hashvalue[ins[node]]+1
		if(outs[node]>=len(ins)):
			print('hash error')
			print(node)
			print(outs[node])
			print(graph.successors(node))
		hashvalue[outs[node]+len(ins)]=hashvalue[outs[node]+len(ins)]+1
	return hashvalue	
	
# classify 3 node motifs
def classifyMotif(graph, nodes, dict, clusters, motifs, edgeDict, classList):
	mini=graph.subgraph(nodes)
	#print(mini.nodes())
	if tuple(hashGraph(mini)) in dict.keys():
		for candidate in dict[tuple(hashGraph(mini))]:
			if nx.is_isomorphic(candidate[0],mini,edge_match=iso.categorical_edge_match('signal', 'a')):
				clusters[candidate[1]].append([mini,candidate[2]])
				motifs[candidate[2],candidate[3]].append(mini)
				if candidate[3] in classList:
					for edge in mini.edges():
						if not edge in edgeDict.keys():
							edgeDict[edge]=0
						edgeDict[edge]=edgeDict[edge]+1

# classify 4 node motifs
def classifyMotif5node(graph, nodes, dict, classifications, motifs, edgeDict, classList):
	mini=graph.subgraph(nodes)
	#print(mini.nodes())
	if tuple(hashGraph(mini)) in dict.keys():
		# print('new candidate')
		for candidate in dict[tuple(hashGraph(mini))]:
			if nx.is_isomorphic(candidate[0],mini,edge_match=iso.categorical_edge_match('signal', 'a')):
				classifications[candidate[1]].append([mini,candidate[2],candidate[3]])
				motifs[candidate[2],candidate[3]].append(mini)
				if candidate[1] in classList:
					for edge in mini.edges():
						if not edge in edgeDict.keys():
							edgeDict[edge]=0
						edgeDict[edge]=edgeDict[edge]+1

				# drawGraph3(mini,'tested.bmp',KEGGdict)
				# drawGraph3(candidate[0],'dict.bmp',KEGGdict)
				# print(mini.edges())
				# print(candidate[0].edges())
				# os.system("pause")

# search 4 node motifs systematically			
def searchMotifs5(givenGraph, dict, clusters,motifs):
	#algorithm taken from Itzhack, Mogilevski, and Louzoun 2007, Physica A (381; 482-490), "An optimal algorithm for countin network motifs" 
	#we iterate through each node, count all motifs involved with that node, remove the node, and move on. 
	i=0
	j=0
	k=0
	L=0
	graph=givenGraph.copy()
	nodelist=graph.nodes()
	# print('searching motifs')
	for v in range(0,len(nodelist)):
		vneighbors=[]
		vneighbors.extend(graph.successors(nodelist[v]))
		vneighbors.extend(graph.predecessors(nodelist[v]))
		while nodelist[v] in vneighbors:
			vneighbors.remove(nodelist[v])
		for i in range(0,len(vneighbors)):
			ineighbors=[]
			ineighbors.extend(graph.successors(vneighbors[i]))
			ineighbors.extend(graph.predecessors(vneighbors[i]))
			while vneighbors[i] in ineighbors:
				ineighbors.remove(vneighbors[i])
			for j in range(i+1,len(vneighbors)):
				jneighbors=[]
				jneighbors.extend(graph.successors(vneighbors[j]))
				jneighbors.extend(graph.predecessors(vneighbors[j]))
				while vneighbors[j] in jneighbors:
					jneighbors.remove(vneighbors[j])
				for k in range(j+1, len(vneighbors)):
					kneighbors=[]
					kneighbors.extend(graph.successors(vneighbors[k]))
					kneighbors.extend(graph.predecessors(vneighbors[k]))
					while vneighbors[k] in kneighbors:
						kneighbors.remove(vneighbors[k])
					if not vneighbors[i]==vneighbors[k] and not vneighbors[i]==vneighbors[j] and not vneighbors[j]==vneighbors[k]:
						#1/4
						for L in range(k+1, len(vneighbors)):
							if(not vneighbors[L]==vneighbors[i] and not vneighbors[L]==vneighbors[j] and not vneighbors[L]== vneighbors[k] ):
								classifyMotif5node(graph,[nodelist[v],vneighbors[i],vneighbors[j],vneighbors[k],vneighbors[L]], dict, clusters,motifs)
				
						union=list(set(kneighbors).union(set(ineighbors).union(set(jneighbors))))
						for L in union:
							if not L in nodelist[v] and not L in vneighbors:
								classifyMotif5node(graph,[nodelist[v],vneighbors[i],vneighbors[j],vneighbors[k],L], dict, clusters,motifs)
				union=list(set(ineighbors).union(set(jneighbors)))
				#the below magically covers 3/8 and 4/8 by considering everything in the next generation as equivalent
				
				for k in range(0,len(union)):
					if not union[k]==nodelist[v] and not union[k]==vneighbors[j] and not union[k]==vneighbors[i] and not vneighbors[i]==vneighbors[j]:
						for L in range(k,len(union)):
							if not union[L]==nodelist[v] and not union[L] in vneighbors:
								classifyMotif5node(graph,[nodelist[v],vneighbors[i],vneighbors[j],union[k],union[L]], dict, clusters,motifs)
						
			for j in range(0,len(ineighbors)):
				if (not ineighbors[j] in vneighbors) and (not ineighbors[j]==nodelist[v]):
					jneighbors=[]
					jneighbors.extend(graph.successors(ineighbors[j]))
					jneighbors.extend(graph.predecessors(ineighbors[j]))
					while ineighbors[j] in jneighbors:
						jneighbors.remove(ineighbors[j])
					for k in range(j+1,len(ineighbors)):
						kneighbors=[]
						kneighbors.extend(graph.successors(ineighbors[k]))
						kneighbors.extend(graph.predecessors(ineighbors[k]))
						while ineighbors[k] in kneighbors:
							kneighbors.remove(ineighbors[k])
						if (not ineighbors[k]in vneighbors) and (not ineighbors[k]==nodelist[v]) and  (not ineighbors[k]==ineighbors[j]):
							#5/8
							for L in range(k+1,len(ineighbors)):
								if not ineighbors[L] in vneighbors and not ineighbors[L]==nodelist[v] and not ineighbors[L]==ineighbors[j] and not ineighbors[L]==ineighbors[k]:
									classifyMotif5node(graph,[nodelist[v],vneighbors[i],ineighbors[j],ineighbors[k], ineighbors[L]], dict, clusters,motifs)
							union=list(set(kneighbors).union(set(jneighbors)))
							#6/8
							for L in range(0,len(union)):
								if not union[L] in ineighbors and not union[L] in vneighbors and not union[L]==nodelist[v] and not union[L]==ineighbors[k] and not union[L]==ineighbors[j]:
									classifyMotif5node(graph,[nodelist[v],vneighbors[i],ineighbors[j],ineighbors[k], union[L]], dict, clusters,motifs)
					for k in range(0,len(jneighbors)):
						kneighbors=[]
						kneighbors.extend(graph.successors(jneighbors[k]))
						kneighbors.extend(graph.predecessors(jneighbors[k]))
						while jneighbors[k] in kneighbors:
							kneighbors.remove(jneighbors[k])
						if (not jneighbors[k] in ineighbors) and (not jneighbors[k] in vneighbors)and not (jneighbors[k]==vneighbors[i]) and not jneighbors[k]==nodelist[v]:
							#7/8
							for L in range(k,len(jneighbors)):
								if not jneighbors[L] in ineighbors and not jneighbors[L]==jneighbors[k] and not jneighbors[L] in vneighbors and not jneighbors[L]==nodelist[v]:
									classifyMotif5node(graph,[nodelist[v],vneighbors[i],ineighbors[j],jneighbors[k],jneighbors[L]], dict, clusters,motifs)
							#8/8
							for L in range(0,len(kneighbors)):
								if not kneighbors[L] in jneighbors and not kneighbors[L] in ineighbors and not kneighbors[L] in vneighbors and not kneighbors[L]==nodelist[v]:
									classifyMotif5node(graph,[nodelist[v],vneighbors[i],ineighbors[j],jneighbors[k],kneighbors[L]], dict, clusters,motifs)
	
		graph.remove_node(nodelist[v])
		# print('removed node '+str(v))
# search 3 node motifs systematically
def searchMotifs4(givenGraph, dict, clusters,motifs):
	#algorithm taken from Itzhack, Mogilevski, and Louzoun 2007, Physica A (381; 482-490), "An optimal algorithm for countin network motifs" 
	#we iterate through each node, count all motifs involved with that node, remove the node, and move on. 
	i=0
	j=0
	k=0
	graph=givenGraph.copy()
	nodelist=graph.nodes()
	for v in range(0,len(nodelist)):
		vneighbors=[]
		vneighbors.extend(graph.successors(nodelist[v]))
		vneighbors.extend(graph.predecessors(nodelist[v]))
		while nodelist[v] in vneighbors:
			vneighbors.remove(nodelist[v])
		
		
		
		for i in range(0,len(vneighbors)):
			ineighbors=[]
			ineighbors.extend(graph.successors(vneighbors[i]))
			ineighbors.extend(graph.predecessors(vneighbors[i]))
			for j in range(i+1,len(vneighbors)):
				jneighbors=[]
				jneighbors.extend(graph.successors(vneighbors[j]))
				jneighbors.extend(graph.predecessors(vneighbors[j]))
				for k in range(j+1, len(vneighbors)):
					if not vneighbors[i]==vneighbors[k] and not vneighbors[i]==vneighbors[j] and not vneighbors[j]==vneighbors[k]:
						classifyMotif(graph,[nodelist[v],vneighbors[i],vneighbors[j],vneighbors[k]], dict, clusters,motifs)
				union=list(set(ineighbors).union(set(jneighbors)))
				for k in union:
					if not k==nodelist[v] and not k==vneighbors[j] and not k==vneighbors[i] and not vneighbors[i]==vneighbors[j]:
						classifyMotif(graph,[nodelist[v],vneighbors[i],vneighbors[j],k], dict, clusters,motifs)
			for j in range(0,len(ineighbors)):
				if (not ineighbors[j] in vneighbors) and (not ineighbors[j]==nodelist[v]):
					jneighbors=[]
					jneighbors.extend(graph.successors(ineighbors[j]))
					jneighbors.extend(graph.predecessors(ineighbors[j]))
					for k in range(j+1,len(ineighbors)):
						if (not ineighbors[k]in vneighbors) and (not ineighbors[k]==nodelist[v]) and  (not ineighbors[k]==ineighbors[j]):
							classifyMotif(graph,[nodelist[v],vneighbors[i],ineighbors[j],ineighbors[k]], dict, clusters,motifs)
					for k in range(0,len(jneighbors)):
						if (not jneighbors[k] in ineighbors) and (not jneighbors[k] in vneighbors)and not (jneighbors[k]==vneighbors[i]) and not jneighbors[k]==nodelist[v]:
							classifyMotif(graph,[nodelist[v],vneighbors[i],ineighbors[j],jneighbors[k]], dict, clusters,motifs)
	graph.remove_node(nodelist[v])

# classify motifs without a signal direction
def classifySingalLessMotif(graph, nodes, hashtable, classifications, motifs, edgeDict, classList):
	
	for node in nodes:
		tempnodes=list(nodes)
		tempnodes.remove(node)
		flag=False
		predecessorList=graph.predecessors(node)
		for tempNode in tempnodes:
			if tempNode in predecessorList:
				predecessorList.remove(tempNode)
		for checkNode in predecessorList:
			flag = True
			
			for tempNode in tempnodes:
				#print(checkNode)
				#print(graph.predecessors(tempNode))
				if checkNode in graph.predecessors(tempNode) or checkNode in graph.successors(tempNode):
					flag=False
					#print('same')
			if flag:
				passNode=checkNode
				break
		if flag:
			tempnodes.append(passNode)
			tempnodes.append(node)
			#print('tempNodes')
			# print('passed value')
			# print(tempnodes) edgeDict, classList
			classifyMotif5node(graph, tempnodes, hashtable, classifications, motifs, edgeDict, classList)
		
# better way of searching 4 node motifs
def searchMotifs4plus(givenGraph, hashtable, clusters,motifs):
	#algorithm taken from Itzhack, Mogilevski, and Louzoun 2007, Physica A (381; 482-490), "An optimal algorithm for countin network motifs" 
	#we iterate through each node, count all motifs involved with that node, remove the node, and move on. 
	classList=['delay2','delay1', 'short']
	edgeDict={}
	i=0
	j=0
	k=0
	graph=givenGraph.copy()
	nodelist=graph.nodes()
	for v in range(0,len(nodelist)):
		vneighbors=[]
		vneighbors.extend(graph.successors(nodelist[v]))
		vneighbors.extend(graph.predecessors(nodelist[v]))
		while nodelist[v] in vneighbors:
			vneighbors.remove(nodelist[v])
		for i in range(0,len(vneighbors)):
			ineighbors=[]
			ineighbors.extend(graph.successors(vneighbors[i]))
			ineighbors.extend(graph.predecessors(vneighbors[i]))
			for j in range(i+1,len(vneighbors)):
				jneighbors=[]
				jneighbors.extend(graph.successors(vneighbors[j]))
				jneighbors.extend(graph.predecessors(vneighbors[j]))
				for k in range(j+1, len(vneighbors)):
					if not vneighbors[i]==vneighbors[k] and not vneighbors[i]==vneighbors[j] and not vneighbors[j]==vneighbors[k]:
						classifyMotif5node(graph,[nodelist[v],vneighbors[i],vneighbors[j],vneighbors[k]], hashtable, clusters,motifs, edgeDict, classList)
				union=list(set(ineighbors).union(set(jneighbors)))
				for k in union:
					if not k==nodelist[v] and not k==vneighbors[j] and not k==vneighbors[i] and not vneighbors[i]==vneighbors[j]:
						classifyMotif5node(graph,[nodelist[v],vneighbors[i],vneighbors[j],k], hashtable, clusters,motifs, edgeDict, classList)
			for j in range(0,len(ineighbors)):
				if (not ineighbors[j] in vneighbors) and (not ineighbors[j]==nodelist[v]):
					jneighbors=[]
					jneighbors.extend(graph.successors(ineighbors[j]))
					jneighbors.extend(graph.predecessors(ineighbors[j]))
					for k in range(j+1,len(ineighbors)):
						if (not ineighbors[k]in vneighbors) and (not ineighbors[k]==nodelist[v]) and  (not ineighbors[k]==ineighbors[j]):
							classifyMotif5node(graph,[nodelist[v],vneighbors[i],ineighbors[j],ineighbors[k]], hashtable, clusters,motifs, edgeDict, classList)
					for k in range(0,len(jneighbors)):
						if (not jneighbors[k] in ineighbors) and (not jneighbors[k] in vneighbors)and not (jneighbors[k]==vneighbors[i]) and not jneighbors[k]==nodelist[v]:
							classifyMotif5node(graph,[nodelist[v],vneighbors[i],ineighbors[j],jneighbors[k]], hashtable, clusters,motifs, edgeDict, classList)
		graph.remove_node(nodelist[v])
		# print('removed node '+str(v))
	return edgeDict

# better way of searching 3 node motifs
def searchMotifs3(givenGraph, dicty, clusters,motifs, classList):
	#algorithm taken from Itzhack, Mogilevski, and Louzoun 2007, Physica A (381; 482-490), "An optimal algorithm for countin network motifs" 
	#we iterate through each node, count all motifs involved with that node, remove the node, and move on. 
	edgeDict={}
	i=0
	j=0
	graph=givenGraph.copy()
	nodelist=list(graph.nodes())
	for v in range(0,len(nodelist)):
		vneighbors=[]
		vneighbors.extend(graph.successors(nodelist[v]))
		vneighbors.extend(graph.predecessors(nodelist[v]))
		for i in range(0,len(vneighbors)):
			ineighbors=[]
			ineighbors.extend(graph.successors(vneighbors[i]))
			ineighbors.extend(graph.predecessors(vneighbors[i]))
			for j in range(i+1,len(vneighbors)):
				if not vneighbors[i]==vneighbors[j]:
					classifyMotif(graph,[nodelist[v],vneighbors[i],vneighbors[j]], dicty, clusters,motifs, edgeDict, classList)
			for j in range(0,len(ineighbors)):
				if (not ineighbors[j] in vneighbors) and (not ineighbors[j]==nodelist[v]) and (not ineighbors[j]==vneighbors[i]):
					classifyMotif(graph,[nodelist[v],vneighbors[i],ineighbors[j]], dicty, clusters,motifs, edgeDict, classList)			
		graph.remove_node(nodelist[v])
	return edgeDict
# search across the 3 node motifs
def motifs3(network):
	classList=[0,1,4,8,15,18]
	hashtable=buildEmptyHashTable()
	clusters, motifs = fillHashTable3node(hashtable)
	edgeDict=searchMotifs3(network,hashtable, clusters, motifs, classList)
	return clusters, motifs, edgeDict
	
def motifs4(num, network):
	hashtable=buildEmptyHashTable(4)
	clusters, motifs = fillHashTable(num,hashtable)
	#print(clusters)
	searchMotifs4(network,hashtable, clusters, motifs)
	return clusters, motifs
	
def motifs5(network):
	print('motifs5')
	hashtable=buildEmptyHashTable()
	classifications, motifs = fillHashTable5node(hashtable)
	print('built hash table')
	#print(clusters)
	edgeDict=searchMotifs4plus(network,hashtable, classifications, motifs)
	return classifications, motifs, edgeDict

def classifyMotifs(motifs):
	counts={}
	picklefile=open('classificationDict.pkl','rb')
	classificationDict=pickle.load(picklefile)
	picklefile.close()
	for value in classificationDict.values():
		counts[value]=0
	for i in range(0,len(motifs)):
		counts[classificationDict[i]]=counts[classificationDict[i]]+len(motifs[i])
	return counts
	
def cut_3_motifs(graphName):
	graph = nx.read_gpickle(graphName+'.gpickle')
	for node in graph.nodes():
		if node in graph.successors(node):
			graph.remove_edge(node,node)

	print(graphName)
	classifications, motifs, edgeDict = motifs3(graph)
	print('start edges')
	print(len(graph.edges()))
	nx.write_graphml(graph, 'before.graphml')
	print('end edges')
	motifList=[]
	for i in [0,1,4,8,15,18]:
		motifList.extend(motifs[('3_node',i)])
	while len(motifList)>0:
		maxnum=0
		for edge in edgeDict.keys():
			if edgeDict[edge]>maxnum:
				maxnum=edgeDict[edge]
				removal=edge
		edgeDict[removal]=0
		graph.remove_edge(removal[0], removal[1])
		tempmotiflist=list(motifList)
		for motif in tempmotiflist:
			if removal in motif.edges():
				motifList.remove(motif) 
				for edge in motif.edges():
					edgeDict[edge]=edgeDict[edge]-1
	nx.write_graphml(graph, graphName+'_3cut.graphml')
	nx.write_gpickle(graph, graphName+'_3cut.gpickle')

	print(len(graph.edges()))

	classifications, motifs, edgeDict = motifs5(graph)
	motifList=classifications['delay2']
	motifList.extend(classifications['delay1'])
	motifList.extend(classifications['short'])
	while len(motifList)>0:
		maxnum=0
		for edge in edgeDict.keys():
			if edgeDict[edge]>maxnum:
				maxnum=edgeDict[edge]
				removal=edge
		edgeDict[removal]=0
		graph.remove_edge(removal[0], removal[1])
		tempmotiflist=list(motifList)
		for motif in tempmotiflist:
			if removal in motif[0].edges():
				motifList.remove(motif) 
				for edge in motif[0].edges():
					edgeDict[edge]=edgeDict[edge]-1
	nx.write_graphml(graph, graphName+'_4cut.graphml')
	nx.write_gpickle(graph, graphName+'_4cut.gpickle')
	print(len(graph.edges()))

def cutMotifs(graph):
	for node in graph.nodes():
		if node in graph.successors(node):
			graph.remove_edge(node,node)

	classifications, motifs, edgeDict = motifs3(graph)
	
	motifList=[]
	for i in [0,1,4,8,15,18]:
		motifList.extend(motifs[('3_node',i)])
	while len(motifList)>0:
		maxnum=0
		for edge in edgeDict.keys():
			if edgeDict[edge]>maxnum:
				maxnum=edgeDict[edge]
				removal=edge
		edgeDict[removal]=0
		graph.remove_edge(removal[0], removal[1])
		tempmotiflist=list(motifList)
		for motif in tempmotiflist:
			if removal in motif.edges():
				motifList.remove(motif) 
				for edge in motif.edges():
					edgeDict[edge]=edgeDict[edge]-1

	classifications, motifs, edgeDict = motifs5(graph)
	motifList=classifications['delay2']
	motifList.extend(classifications['delay1'])
	motifList.extend(classifications['short'])
	while len(motifList)>0:
		maxnum=0
		for edge in edgeDict.keys():
			if edgeDict[edge]>maxnum:
				maxnum=edgeDict[edge]
				removal=edge
		edgeDict[removal]=0
		graph.remove_edge(removal[0], removal[1])
		tempmotiflist=list(motifList)
		for motif in tempmotiflist:
			if removal in motif[0].edges():
				motifList.remove(motif) 
				for edge in motif[0].edges():
					edgeDict[edge]=edgeDict[edge]-1
	return graph

if __name__ == '__main__':
	graphNames=['hsa04066','hsa04350','hsa04380','hsa04612','hsa04630','hsa04650','hsa04658','hsa04659','hsa04660','Liu_1','Liu_2','Liu_3']
	for name in graphNames:
		cut_3_motifs(name)