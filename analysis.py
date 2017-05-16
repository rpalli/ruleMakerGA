import matplotlib.pyplot as plt
import numpy as numpy
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns

def compareIndividualsNodeWise(truthList, testList, model):
	nodesensitivity=[]
	nodespecificity=[]
	netOnes=[]
	netZeros=[]
	netNegOnes=[]
	for node in range(len(model.nodeList)):
		ones=[]
		zeros=[]
		negones=[]
		# get indices to extract individual for this node
		if node==len(model.nodeList)-1:
			end=len(model.nodeList)
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		sumindividual=[]
		#loop over individuals provided and calculate relevant values
		for i in range(len(truthList)):
			truth= truthList[i][start:end]
			test= testList[i][start:end]
			sumindividual.append(numpy.sum(truth))
			newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
			ones.append(newindividual.count(1))
			zeros.append(newindividual.count(0))
			negones.append(newindividual.count(-1))
		# append node-wise breakdowns to list of breakdowns for the model as a whole
		netOnes.append(numpy.sum(ones))
		netZeros.append(numpy.sum(zeros))
		netNegOnes.append(numpy.sum(negones))
		
		# calculate sensitivity and specificity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(1.*numpy.sum(temp)/len(temp))
		temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			specificity=100
		else:
			specificity=(1.*numpy.sum(temp)/len(temp))
		# add to list of sensitivity and specificity by node
		nodesensitivity.append(sensitivity)
		nodespecificity.append(specificity)
	#calculate sensitivity and specificity on the network as a whole
	ones=[]
	zeros=[]
	negones=[]
	sumindividual=[]

	for i in range(len(truthList)):
		truth= truthList[i]
		test= testList[i]
		sumindividual.append(1.*numpy.sum(truth))
		newindividual=[a_i - b_i for a_i, b_i in zip(truth, test)]
		ones.append(newindividual.count(1))
		zeros.append(newindividual.count(0))
		negones.append(newindividual.count(-1))
	temp=[100 if (sumindividual[i])==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		sensitivity=100
	else:
		sensitivity=(1.*numpy.sum(temp)/len(temp))
	temp=[100 if (len(newindividual)-sumindividual[i])==0 else (1.*len(newindividual)-sumindividual[i]-negones[i])/(len(newindividual)-sumindividual[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		specificity=100
	else:
		specificity=(1.*numpy.sum(temp)/len(temp))
	return sensitivity, specificity, nodesensitivity, nodespecificity

class modelHolder:
	def __init__(self,valueList):
		[ evaluateNodes,individualParse, andNodeList, andNodeInvertList, andLenList,earlyEvalNodes,nodeList, nodeDict, initValueList]=valueList
		self.evaluateNodes=evaluateNodes
		self.individualParse=individualParse
		self.andNodeList=andNodeList
		self.andNodeInvertList=andNodeInvertList
		self.andLenList=andLenList
		self.earlyEvalNodes=earlyEvalNodes
		self.nodeList=nodeList
		self.nodeDict=nodeDict
		self.initValueList=initValueList

def analyzeRewire(fileName):
	# upload data for this replicate
	outputList=pickle.Unpickler(open( fileName, "rb" )).load()
	[truthList,testList, devList,size, evaluateNodes,individualParse, andNodeList, andNodeInvertList, andLenList,earlyEvalNodes,nodeList, nodeDict, initValueList]=outputList
	model=modelHolder([evaluateNodes,individualParse, andNodeList, andNodeInvertList, andLenList,earlyEvalNodes,nodeList, nodeDict, initValueList])

	truthIndividuals=[]
	truthAvgs=[]
	hofScores=[]
	newSSS=[]
	trueNodeList=[]
	zeros=[]
	ones=[]
	negones=[]
	truthCounter=0
	sumindividual=[]

	print('devList')
	print(devList)

	tuple1=compareIndividualsNodeWise(truthList, testList, model)
	sensitivities=[]
	specificities=[]
	inEdgeNums=[]
	overlaps=[]
	# run correction that simplifies truth rules
	TPsum=0
	TNsum=0
	FNsum=0
	FPsum=0
	newtruths=[]
	for k in range(len(truthList)):
		newtruths.append([])
	newtests=[]
	for k in range(len(testList)):
		newtests.append([])
	for node in range(0,len(model.nodeList)):
		FP=0
		TP=0
		TN=0
		FN=0
		if node==(len(model.nodeList)-1):
			end=len(truthList[0])-1
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		andNodeList=model.andNodeList[node]
		inEdges=[]
		for lister in andNodeList:
			inEdges.append(set(lister))
		for k in range(len(truthList)):
			truth=truthList[k][start:end]
			test= testList[k][start:end]
			for i in range(len(truth)):
				if truth[i]==1:
					for j in range(len(truth)):
						if truth[j]==1 and not i==j:
							if inEdges[i].issubset(inEdges[j]):
								truth[j]=0
			for i in range(len(test)):
				if test[i]==1:
					for j in range(len(test)):
						if test[j]==1 and not i==j:
							if inEdges[i].issubset(inEdges[j]):
								test[j]=0
			newtests[k].extend(test)
			newtruths[k].extend(truth)
		for k in range(len(truthList)):
			truth=newtruths[k][start:end]
			test= newtests[k][start:end]
			truthSet=Set([])
			testSet=Set([])
			baseSet=Set([])
			for i in range(0,len(truth)):
				if truth[i]==1:
					for nodeToAdd in andNodeList[i]:
						truthSet.add(nodeToAdd)
			for i in range(0,len(test)):
				if test[i]==1:
					for nodeToAdd in andNodeList[i]:
						testSet.add(nodeToAdd)
				for nodeToAdd in andNodeList[i]:
					baseSet.add(nodeToAdd)
			FP+=1.*len(testSet.difference(truthSet))
			TP+=1.*len(testSet.intersection(truthSet))
			TN+=1.*(len((baseSet.difference(truthSet)).difference(testSet)))
			FN+=1.*len(truthSet.difference(testSet))
		if (TP+FN)>0:
			sensitivity=1.*TP/(TP+FN)
		else:
			sensitivity=100
		if TN+FP>0:
			specificity=1.*TN/(TN+FP)
		else:
			specificity=100
		sensitivities.append(sensitivity)
		specificities.append(specificity)
		TPsum+=TP
		TNsum+=TN
		FNsum+=FN
		FPsum+=FP
		inEdgeNums.append(len(baseSet))
		overlaper=-1
		if len(baseSet)==2:
			overlapSet=Set([])
			lappersets=[]
			inneredge=[]
			for upstream in baseSet:
				lapperset=Set([])
				if upstream==(len(model.nodeList)-1):
					end=len(truthList[0])-1
				else:
					end=model.individualParse[upstream+1]
				start=model.individualParse[upstream]
				andNodeList=model.andNodeList[upstream]
				for lister in andNodeList:
					lapperset.add(Set(lister))
				lappersets.append(lapperset)
			overlaper=len(lappersets[0].intersection(lappersets[1]))
			if overlaper==0:
				for upstream in baseSet:
					if upstream in lappersets[0] or upstream in lappersets[1]:
						print("feed-forward")
						overlaper=1
		overlaps.append(overlaper)
	tuple2=compareIndividualsNodeWise(newtruths, testList, model)
	if (TPsum+FNsum)>0:
		sensitivity=1.*TPsum/(TPsum+FNsum)
	else:
		sensitivity=100
	if (FPsum+TNsum)>0:
		specificity= 1.*TNsum/(FPsum+TNsum)
	else:
		specificity=100
	tuple3= (sensitivity, specificity, sensitivities, specificities)

	tuple4=compareIndividualsNodeWise(newtruths, newtests, model)
	for i in range(len(newtruths)):
		print("truth, found")
		print(utils.writeModel(newtruths[i], model))
		print(utils.writeModel(newtests[i], model))

	return tuple1, tuple2, tuple3, tuple4, devList, inEdgeNums, overlaps, [testList, truthList, newtruths, newtests]


def analyzeExperiment():
	sensitivities=[]
	specificities=[]
	specificityStds=[]
	sensitivityStds=[]
	nodesensitivities=[[],[],[],[]]
	nodespecificities=[[],[],[],[]]
	truthholder=[]
	edgeDegree=[]
	overlapNodes=[]

	#FOR .gpickle!
	for code in lines:
		# graph=nx.DiGraph()
		# coder=str('ko'+code[:-1])
		# nc.uploadKEGGcodes([coder], graph, dict2)
		# coder=str('hsa'+code[:-1])
		# nc.uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		# if(len(graph.edges())>1):
		# 	graph=nc.simplifyNetwork(graph, data)
		#graph = simpleNetBuild()
		graph=liu.LiuNetwork1Builder()
		# coder='liu'
	
		print(coder)
		print(len(graph.edges()))
		
		nx.write_graphml(graph,coder+'.graphml')
		tempsensitivities=[[],[],[],[]]
		tempspecificities=[[],[],[],[]]
		truthlists=[[],[],[],[]]
		devLists=[]

		#FOR EACH RUN WITH THAT GPICKLE
		for i in range(bioReplicates):
			tuple1, tuple2, tuple3, tuple4, devList, inEdgeNums, overlaps, truthlisttemp=analyzeRewire(name1, name2,name3, replicates)
			sensitivity1, specificity1, nodesensitivity1, nodespecificity1 = tuple1
			sensitivity2, specificity2, nodesensitivity2, nodespecificity2 = tuple2
			sensitivity3, specificity3, nodesensitivity3, nodespecificity3 = tuple3
			sensitivity4, specificity4, nodesensitivity4, nodespecificity4 = tuple4
			tempsensitivities[0].append(sensitivity1)
			tempsensitivities[1].append(sensitivity2)
			tempsensitivities[2].append(sensitivity3)
			tempsensitivities[3].append(sensitivity4)
			tempspecificities[0].append(specificity1)
			tempspecificities[1].append(specificity2)
			tempspecificities[2].append(specificity3)
			tempspecificities[3].append(specificity4)
			nodesensitivities[0].extend(nodesensitivity1)
			nodesensitivities[1].extend(nodesensitivity2)
			nodesensitivities[2].extend(nodesensitivity3)
			nodesensitivities[3].extend(nodesensitivity4)
			nodespecificities[0].extend(nodespecificity1)
			nodespecificities[1].extend(nodespecificity2)
			nodespecificities[2].extend(nodespecificity3)
			nodespecificities[3].extend(nodespecificity4)
			devLists.append(devList)
			overlapNodes.extend(overlaps)
			edgeDegree.extend(inEdgeNums)
			truthlists[0].extend(truthlisttemp[0])
			truthlists[1].extend(truthlisttemp[1])
			truthlists[2].extend(truthlisttemp[2])
			truthlists[3].extend(truthlisttemp[3])
		for i in range(len(tempsensitivities)):
			tempsensitivities[i]=filter(lambda a: a != 100, tempsensitivities[i])
			if len(tempsensitivities[i])==0:
				tempsensitivities[i].append(100)
		for tempholder in tempspecificities:
			tempspecificities[i]=filter(lambda a: a != 100, tempspecificities[i])
			if len(tempspecificities[i])==0:
				tempspecificities[i].append(0.)
		sensitivity=[1.*numpy.sum(tempsensitivities[0])/len(tempsensitivities[0]),1.*numpy.sum(tempsensitivities[1])/len(tempsensitivities[1]),1.*numpy.sum(tempsensitivities[2])/len(tempsensitivities[2]),1.*numpy.sum(tempsensitivities[3])/len(tempsensitivities[3])]
		sensitivityStd=[1.*numpy.std(tempsensitivities[0])/len(tempsensitivities[0]),1.*numpy.std(tempsensitivities[1])/len(tempsensitivities[1]),1.*numpy.std(tempsensitivities[2])/len(tempsensitivities[2]),1.*numpy.std(tempsensitivities[3])/len(tempsensitivities[3])]
		specificity=[1.*numpy.sum(tempspecificities[0])/len(tempspecificities[0]),1.*numpy.sum(tempspecificities[1])/len(tempspecificities[1]),1.*numpy.sum(tempspecificities[2])/len(tempspecificities[2]),1.*numpy.sum(tempspecificities[3])/len(tempspecificities[3])]
		specificityStd=[1.*numpy.std(tempspecificities[0])/len(tempspecificities[0]),1.*numpy.std(tempspecificities[1])/len(tempspecificities[1]),1.*numpy.std(tempspecificities[2])/len(tempspecificities[2]),1.*numpy.std(tempspecificities[3])/len(tempspecificities[3])]
		truthholder.append(truthlists)
		specificities.append(specificity)
		sensitivities.append(sensitivity)
		specificityStds.append(specificityStd)
		sensitivityStds.append(sensitivityStd)
		print(sensitivity)
		print(sensitivityStd)
		print(specificity)
		print(specificityStd)
		print(nodesensitivities)
		print(nodespecificities)
		print(devLists)
	nodeLookup={}
	for number in edgeDegree:
		nodeLookup[number]=[[],[],[],[],[],[],[],[]]
	for i in range(len(edgeDegree)):
		nodeLookup[edgeDegree[i]][0].append(nodesensitivities[0][i])
		nodeLookup[edgeDegree[i]][1].append(nodesensitivities[1][i])
		nodeLookup[edgeDegree[i]][2].append(nodesensitivities[2][i])
		nodeLookup[edgeDegree[i]][3].append(nodesensitivities[3][i])
		nodeLookup[edgeDegree[i]][4].append(nodespecificities[0][i])
		nodeLookup[edgeDegree[i]][5].append(nodespecificities[1][i])
		nodeLookup[edgeDegree[i]][6].append(nodespecificities[2][i])
		nodeLookup[edgeDegree[i]][7].append(nodespecificities[3][i])
	overlapLookup={}
	for overlap in overlapNodes:
		overlapLookup[overlap]=[[],[],[],[],[],[],[],[]]
	for i in range(len(overlapNodes)):
		overlapLookup[overlapNodes[i]][0].append(nodesensitivities[0][i])
		overlapLookup[overlapNodes[i]][1].append(nodesensitivities[1][i])
		overlapLookup[overlapNodes[i]][2].append(nodesensitivities[2][i])
		overlapLookup[overlapNodes[i]][3].append(nodesensitivities[3][i])
		overlapLookup[overlapNodes[i]][4].append(nodespecificities[0][i])
		overlapLookup[overlapNodes[i]][5].append(nodespecificities[1][i])
		overlapLookup[overlapNodes[i]][6].append(nodespecificities[2][i])
		overlapLookup[overlapNodes[i]][7].append(nodespecificities[3][i])
	finalNodeData=[]
	finalExtendData=[]
	for key in nodeLookup.keys():
		templisting=[]
		tempExtended=[]
		for lister in nodeLookup[key]:
			newlist=filter(lambda a: a != 100, lister)
			tempExtended.append(newlist)
			if len(newlist)==0:
				templisting.append(100)
			else:
				templisting.append(1.* numpy.sum(newlist)/len(newlist))
		finalNodeData.append(templisting)
		finalExtendData.append(tempExtended)
	finalOverlapExtendData=[]
	for key in overlapLookup.keys():
		templisting=[]
		tempOverlap=[]
		for lister in overlapLookup[key]:
			newlist=filter(lambda a: a != 100, lister)
			tempOverlap.append(newlist)
		finalOverlapExtendData.append(tempOverlap)
	print(nodeLookup.keys())
	print(overlapLookup.keys())
	print(finalNodeData)
	print(finalExtendData) 

def plotResults():
	networkNodeNums=[18,32,7,23,68,78,72,28]

	finalNodeData=pickle.Unpickler(open( "node_by_node_data.pickle", "rb" )).load()
	extendedNodeData=pickle.Unpickler(open( "extended_node_by_node_data.pickle", "rb" )).load()
	[sensitivities,sensStd, specificities, specsStd]=pickle.Unpickler(open( "network_by_network_data.pickle", "rb" )).load()



	sensitivities=[sens[3] for sens in sensitivities]
	specificities=[spec[3] for spec in specificities]
	senStd=[sens[3] for sens in sensStd]
	specStd=[sens[3] for sens in specsStd]
	specStd.pop(0)
	senStd.pop(0)
	sensitivities.pop(0)
	specificities.pop(0)
	print(sensitivities)
	sns.set_style("ticks")
	sns.set(font_scale=2) 
	sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
	threshold=.8

	y_pos = np.arange(len(networkNodeNums))
	plt.bar(y_pos, sensitivities, align='center', alpha=0.5,yerr=senStd, capsize=20)
	plt.xticks(y_pos, networkNodeNums)
	plt.ylabel('Sensitivity')
	plt.xlabel('Node number')
	plt.title('Sensitivity by Nodes in Pathway')
	plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
	plt.tight_layout()
	sns.despine()
	plt.savefig('IFNG_network_sensitivities.png', bbox_inches='tight')

	plt.clf()

	plt.bar(y_pos, specificities, align='center', alpha=0.5,yerr=specStd)
	plt.xticks(y_pos, networkNodeNums)

	plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
	plt.ylabel('Specificity')
	plt.xlabel('Node number')
	plt.title('Specificity by Nodes in Pathway')
	sns.despine()
	plt.tight_layout()
	plt.savefig('IFNG_network_specificities.png', bbox_inches='tight')







	plt.clf()

	# nodeNums=[2, 3, 4, 5, 6, 7, 9]
	nodeNums=[0,1,2,3,4,5,6]

	nodeArray=[]
	extendedNodeData.pop(0)
	for tempExtended in extendedNodeData:
		nodeArray.append(tempExtended[3])
	for array in nodeArray:
		if len(array)==0:
			array.append(0)
	# nodeArray.pop(0)
	nodeArray.pop(0)
	nodeNums.pop(0)
	# nodeArray.pop(0)

	dfdictlist=[]
	for i in range(len(nodeArray)):
		for element in nodeArray[i]:
			dfdictlist.append({'In-degree':nodeNums[i],'Sensitivity':element})
	df=pd.DataFrame(dfdictlist)

	# plt.title('Sensitivity by in-degree')

	#print(df)
	sns.set(font_scale=2) 
	sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})

	ax=sns.swarmplot(data=df,x='In-degree', y='Sensitivity', color='orange')

	plt.ylabel('Sensitivity')
	plt.xlabel('In-degree')
	plt.title('Sensitivity by in-degree', color='white')
	plt.ylim([0,1])
	plt.tight_layout()
	plt.ylim([0,1])
	plt.savefig('IFNG_clear_sens.png', bbox_inches='tight', transparent=True)
	plt.clf()


	nodeNums=[0,1,2,3,4,5,6]

	nodeArray=[]
	for tempExtended in extendedNodeData:
		nodeArray.append(tempExtended[7])
	for array in nodeArray:
		if len(array)==0:
			array.append(0)
	nodeArray.pop(0)
	nodeNums.pop(0)
	nodeArray.pop(0)
	nodeNums.pop(0)
	# nodeArray.pop(0)


	dfdictlist=[]
	for i in range(len(nodeArray)):
		for element in nodeArray[i]:
			dfdictlist.append({'In-degree':nodeNums[i],'Specificity':element})
	df=pd.DataFrame(dfdictlist)


	sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})


	plt.xlabel('In-degree', color='white')
	# plt.title('Sensitivity by in-degree')

	sns.set(font_scale=2) 


	ax=sns.swarmplot(data=df,x='In-degree', y='Specificity', color='orange')
	# ax.set_xticklabels(nodeNums, color='white')
	# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')

	plt.ylabel('Specificity')
	plt.xlabel('In-degree')
	plt.title('Specificity by in-degree', color='white')
	plt.ylim([0,1])
	plt.tight_layout()
	# sns.despine()

	plt.savefig('IFNG_clear_spec.png', bbox_inches='tight', transparent=True)


	plt.clf()


print(analyzeRewire('hsa04066_6_output.pickle'))