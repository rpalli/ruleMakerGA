# import python packages
import matplotlib.pyplot as plt
import numpy as numpy
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from sets import Set
import os as os
# import other parts of our code
import utils as utils
def compareIndividualsNodeWise(truthList, testList, model):
	nodesensitivity=[]
	nodePPV=[]
	netOnes=[]
	netZeros=[]
	netNegOnes=[]
	ruleTruths=[]
	SDs=[]
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
		

		# calculate % rules correct and structural distance
		temp=[1 if (ones[i]+negones[i])==0 else 0 for i in range(0,len(ones))]
		ruleTruths.append(1.*numpy.sum(temp)/len(temp))
		temp=[ones[i]+negones[i] for i in range(0,len(ones))]
		SDs.append(1.*numpy.sum(temp)/len(temp))

		# calculate sensitivity and specificity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(1.*numpy.sum(temp)/len(temp))

		temp=[100 if (sumindividual[i]-ones[i]+negones[i])==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]-ones[i]+negones[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			PPV=100
		else:
			PPV=(1.*numpy.sum(temp)/len(temp))
		# add to list of sensitivity and specificity by node
		nodesensitivity.append(sensitivity)
		nodePPV.append(PPV)
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
	temp=[100 if (len(newindividual)-sumindividual[i])==0 else 1.*(sumindividual[i]-ones[i])/((sumindividual[i]-ones[i])+negones[i]) for i in range(0,len(ones))]
	temp=filter(lambda a: a != 100, temp)
	if len(temp)==0:
		PPV=100
	else:
		PPV=(1.*numpy.sum(temp)/len(temp))
	return sensitivity, PPV, nodesensitivity, nodePPV, ruleTruths, SDs

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
	tuple1=compareIndividualsNodeWise(truthList, testList, model)
	sensitivities=[]
	specificities=[]
	inEdgeNums=[]
	# run correction that simplifies truth rules
	TPsum=0
	TNsum=0
	FNsum=0
	FPsum=0
	PPVs=[]
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
		if TP+FP>0:
			PPV=1.*TP/(TP+FP)
		else:
			PPV=100
		sensitivities.append(sensitivity)
		PPVs.append(PPV)
		TPsum+=TP
		TNsum+=TN
		FNsum+=FN
		FPsum+=FP
		inEdgeNums.append(len(baseSet))
	if (TPsum+FNsum)>0:
		sensitivity=1.*TPsum/(TPsum+FNsum)
	else:
		sensitivity=100
	if (FPsum+TPsum)>0:
		PPV= 1.*TPsum/(FPsum+TPsum)
	else:
		PPV=100
	tuple2=compareIndividualsNodeWise(newtruths, testList, model)
	tuple3= (sensitivity, PPV, sensitivities, PPVs)
	tuple4=compareIndividualsNodeWise(newtruths, newtests, model)
	# for i in range(len(newtruths)):
	# 	print("truth, found")
	# 	print(utils.writeModel(newtruths[i], model))
	# 	print(utils.writeModel(newtests[i], model))
	return tuple1, tuple2, tuple3, tuple4, devList, inEdgeNums, [testList, truthList, newtruths, newtests], len(model.andNodeList), model

def analyzeExperiment(codes):
	sensitivities=[]
	specificities=[]
	specificityStds=[]
	sensitivityStds=[]
	nodesensitivities=[[],[],[],[]]
	nodespecificities=[[],[],[],[]]
	truthholder=[]
	degree2dictlist=[]

	combinedInDegrees={}
	ruleTruths=[]
	SDs=[]
	nodeNumbers=[]
	# for each gpickle
	for code in codes:
		print(code)		
		nodeLookup={}
		tempsensitivities=[[],[],[],[]]
		tempspecificities=[[],[],[],[]]
		truthlists=[[],[],[],[]]
		devLists=[]
		ruletruthtemp=[]
		SDtemp=[]
		#FOR EACH RUN WITH THAT GPICKLE
		for i in range(1,11):
			edgeDegree=[]
			nodesensitivities=[[],[],[],[]]
			nodespecificities=[[],[],[],[]]


			tuple1, tuple2, tuple3, tuple4, devList, inEdgeNums,  truthlisttemp, nodeNum, model=analyzeRewire('pickles/'+code+'_'+str(i)+'_output.pickle')
			if i==1:
				nodeNumbers.append(nodeNum)
				print(nodeNum)
			sensitivity1, specificity1, nodesensitivity1, nodespecificity1, ruleTruths1, SDs1 = tuple1
			sensitivity2, specificity2, nodesensitivity2, nodespecificity2, ruleTruths2, SDs2 = tuple2
			sensitivity3, specificity3, nodesensitivity3, nodespecificity3 = tuple3
			sensitivity4, specificity4, nodesensitivity4, nodespecificity4, ruleTruth, SD = tuple4
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
			ruletruthtemp.extend(ruleTruth)
			SDtemp.extend(SD)
			devLists.append(devList)
			edgeDegree.extend(inEdgeNums)
			truthlists[0].extend(truthlisttemp[0])
			truthlists[1].extend(truthlisttemp[1])
			truthlists[2].extend(truthlisttemp[2])
			truthlists[3].extend(truthlisttemp[3])
			for number in edgeDegree:
				if not number in nodeLookup.keys():
					nodeLookup[number]=[[],[],[],[],[],[],[],[]]
				if not number in combinedInDegrees.keys():
					combinedInDegrees[number]=[[],[],[],[],[],[],[],[]]
			for i in range(len(edgeDegree)):
				for j in range(0,4):
					nodeLookup[edgeDegree[i]][j].append(nodesensitivities[j][i])
				for j in range(0,4):
					nodeLookup[edgeDegree[i]][j+4].append(nodespecificities[j][i])
			for number in nodeLookup.keys():
				for i in range(0,8):
					tempcomb=filter(lambda a: a != 100, nodeLookup[number][i])
					combinedInDegrees[number][i].append(numpy.mean(tempcomb))
			if 2 in nodeLookup.keys():
				tempcomb=filter(lambda a: a != 100, nodeLookup[2][3])
				tempcomb1=filter(lambda a: a != 100, nodeLookup[2][7])
				degree2dictlist.append({'Node_Num':nodeNum,'Sensitivity':numpy.mean(tempcomb), 'PPV': numpy.mean(tempcomb1)})
			if code=='hsa04630':
				print(model.andNodeList)
		SDs.append(numpy.mean(SDtemp))
		ruleTruths.append(numpy.mean(ruletruthtemp))
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
		# print(sensitivity)
		# print(sensitivityStd)
		# print(specificity)
		# print(specificityStd)
		# print(nodesensitivities)
		# print(nodespecificities)
		# print(devLists)
	print(SDs)
	print(ruleTruths)

		
	finalNodeData=[]
	finalExtendData=[]
	for key in combinedInDegrees.keys():
		templisting=[]
		tempExtended=[]
		for lister in combinedInDegrees[key]:
			newlist=filter(lambda a: a != 100, lister)
			tempExtended.append(newlist)
			if len(newlist)==0:
				templisting.append(100)
			else:
				templisting.append(1.* numpy.sum(newlist)/len(newlist))
		finalNodeData.append(templisting)
		finalExtendData.append(tempExtended)


	sensitivities=[sens[3] for sens in sensitivities]
	specificities=[spec[3] for spec in specificities]
	senStd=[sens[3] for sens in sensitivityStds]
	specStd=[sens[3] for sens in specificityStds]
	# specStd.pop(0)
	# senStd.pop(0)
	# sensitivities.pop(0)
	# specificities.pop(0)
	print(sensitivities)
	sns.set_style("ticks")
	sns.set(font_scale=2) 
	sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
	threshold=.8
	# nodeNumbers.pop(0)
	y_pos = np.arange(len(nodeNumbers))
	plt.bar(y_pos, sensitivities, align='center', alpha=0.5,yerr=senStd, capsize=20)
	plt.xticks(y_pos, nodeNumbers)
	plt.ylabel('Sensitivity')
	plt.xlabel('Node number')
	plt.title('Sensitivity by Nodes in Pathway')
	plt.plot( [-.5,len(nodeNumbers)-.5],[threshold, threshold], "k--")
	plt.tight_layout()
	sns.despine()
	plt.savefig('IFNG_network_sensitivities.png', bbox_inches='tight')

	plt.clf()

	plt.bar(y_pos, specificities, align='center', alpha=0.5,yerr=specStd)
	plt.xticks(y_pos, nodeNumbers)

	plt.plot( [-.5,len(nodeNumbers)-.5],[threshold, threshold], "k--")
	plt.ylabel('PPV')
	plt.xlabel('Node number')
	plt.title('PPV by Nodes in Pathway')
	sns.despine()
	plt.tight_layout()
	plt.savefig('IFNG_network_PPV.png', bbox_inches='tight')


	plt.clf()


	plt.bar(y_pos, ruleTruths, align='center', alpha=0.5)
	plt.xticks(y_pos, nodeNumbers)

	plt.plot( [-.5,len(nodeNumbers)-.5],[threshold, threshold], "k--")
	plt.ylabel('Specificity')
	plt.xlabel('Node number')
	plt.title('% Rules True by Nodes in Pathway')
	sns.despine()
	plt.tight_layout()
	plt.savefig('IFNG_network_Percent_T_Rules.png', bbox_inches='tight')


	plt.clf()




	# nodeNums=[2, 3, 4, 5, 6, 7, 9]
	nodeNums=[0,1,2,3,4,5,6]

	nodeArray=[]
	for tempExtended in finalExtendData:
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
			if nodeNums[i]!=1:
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
	for tempExtended in finalExtendData:
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

	df=pd.DataFrame(degree2dictlist)
	sns.set(font_scale=2) 
	sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})

	ax=sns.swarmplot(data=df,x='Node_Num', y='Sensitivity', color='orange')

	plt.ylabel('Sensitivity')
	plt.xlabel('Node Number')
	plt.title('Sensitivity at in-degree 2 by Node Number', color='white')
	plt.ylim([0,1])
	plt.tight_layout()
	plt.ylim([0,1])
	plt.savefig('IFNG_sens_node_num.png', bbox_inches='tight', transparent=True)
	plt.clf()

	sns.set(font_scale=2) 
	sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})

	ax=sns.swarmplot(data=df,x='Node_Num', y='PPV', color='orange')

	plt.ylabel('PPV')
	plt.xlabel('Node Number')
	plt.title('PPV at in-degree 2 by Node Number', color='white')
	plt.ylim([0,1])
	plt.tight_layout()
	plt.ylim([0,1])
	plt.savefig('IFNG_PPV_node_num.png', bbox_inches='tight', transparent=True)
	plt.clf()




if __name__ == '__main__':
	codes=[]
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			codes.append(file[:-8])
	print(codes)
	analyzeExperiment(codes)
