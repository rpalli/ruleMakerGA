
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
from scipy.stats import variation
import networkx as nx
# import other parts of our code
import utils as utils

def styleSetter():
	sns.set(font_scale=2) 
	sns.set_style("whitegrid")
	# sns.set_context("talk")
	# sns.set_style("dark")
	# sns.axes_style(style='white', rc={'xtick.color': 'white','axes.labelcolor': 'white','text.color': 'white'})
	# plt.ylabel( color='white')
	# plt.xlabel( color='white')
	# plt.xticks( color='white')
	# plt.yticks( color='white')
	# sns.set_context("talk", rc={'xtick.color': '.1','axes.labelcolor': 'white','text.color': 'white'})
	# sns.set_style("dark")
def plot_df01swarm(dictlist, xval, yval, title, xlabeler, ylabeler, plotname):
	styleSetter()
	df=pd.DataFrame(dictlist)
	ax=sns.swarmplot(data=df,x=xval, y=yval, color='black')
	# ax.set_xticklabels(nodeNums, color='white')
	# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	# plt.title(title, color='black')
	# plt.ylim([0,1])
	plt.tight_layout()
	# sns.despine()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()
def plot_scatter(dictlist, xval, yval, title, xlabeler, ylabeler, plotname):

	df=pd.DataFrame(dictlist)
	styleSetter()
	ax=sns.lmplot(data=df,x=xval, y=yval, fit_reg=True, )
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	plt.tight_layout()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()
def plot_scatter2(dictlist, xval, yval, title, xlabeler, ylabeler, hue, plotname):

	df=pd.DataFrame(dictlist)
	styleSetter()
	ax=sns.lmplot(data=df,x=xval, y=yval,hue=hue, fit_reg=False)
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()
def jointplotter(dictlist, xval, yval, title, xlabeler, ylabeler, plotname,xlims):
	df=pd.DataFrame(dictlist)
	styleSetter()
	ax=sns.jointplot(data=df,x=xval, y=yval, kind='kde')
	plt.tight_layout()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()
def plot_df01box(dictlist, xval, yval, title, xlabeler, ylabeler, plotname, hueBool, huer):
	df=pd.DataFrame(dictlist)
	sns.set_palette(sns.color_palette("Greys", 6))
	styleSetter()
	if not hueBool:
		ax=sns.barplot(data=df,x=xval, y=yval, color=sns.color_palette("Greys", 4)[0], ci=68.26, capsize=.07)
	else:
		ax=sns.barplot(data=df,x=xval, y=yval, color='Gray', hue=huer, ci=68.26, capsize=.07)
		plt.legend(bbox_to_anchor=(1.1, .8), bbox_transform=plt.gcf().transFigure)
	sns.despine()
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()
def pointplot(dictlist, xval, yval, title, xlabeler, ylabeler, plotname, hueBool, huer):
	df=pd.DataFrame(dictlist)
	#sns.set_style("dark")
	sns.set_palette(sns.color_palette("Greys", 6))
	sns.set(font_scale=2) 
	sns.set_style("ticks")
	if not hueBool:
		ax=sns.regplot(data=df,x=xval, y=yval, color='black', fit_reg=False)
	else:
		ax=sns.regplot(data=df,x=xval, y=yval, fit_reg=False)
		plt.legend(bbox_to_anchor=(1.1, .8), bbox_transform=plt.gcf().transFigure)
	sns.despine()
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	# plt.title(title, color='black')
	# plt.ylim([0,1])
	# plt.tight_layout()
	# sns.despine()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()

def findEnds(model, node, indiv):
	if node==len(model.nodeList)-1:
		end1=len(indiv)
	else:
		end1=model.individualParse[node+1]
	start1=model.individualParse[node]
	return start1, end1

def findUpstreamChars(stem, nodePPV, nodeSens, nodeRuleTruth, maxRT, model):
	graph = nx.read_gpickle('gpickles/'+stem[8:-1]+'.gpickle')
	sources=[]
	for element in graph.nodes():
		if graph.in_degree()[element]==0:
			sources.append(element)
	upstreamPPV, upstreamSens, upstreamRT, upstreamMaxRT=[],[],[], []
	for node in model.nodeList:
		tempPPV=[]
		tempSens=[]
		tempRT=[]
		tempMaxRT=[]
		templist= graph.predecessors(node)
		for tempNode in templist:
			if (nodePPV[model.nodeDict[tempNode]] <100):
				tempPPV.append(nodePPV[model.nodeDict[tempNode]])
				tempSens.append(nodeSens[model.nodeDict[tempNode]])
				tempRT.append(nodeRuleTruth[model.nodeDict[tempNode]])
				tempMaxRT.append(maxRT[model.nodeDict[tempNode]])
		if len(tempPPV)>0:
			upstreamPPV.append(numpy.mean(tempPPV))
			upstreamSens.append(numpy.mean(tempSens))
			upstreamRT.append(numpy.mean(tempRT))
			upstreamMaxRT.append(numpy.mean(tempMaxRT))
		else:
			upstreamPPV.append(1.)
			upstreamSens.append(1.)
			upstreamRT.append(1.)
			upstreamMaxRT.append(1.)

	return upstreamPPV, upstreamSens, upstreamRT, upstreamMaxRT
def compareIndividualsNodeWise(truthList, testList, model1s, model2s,covs, equivs):

	modeler=model1s[0]
	SDs=[0. for q in truthList]
	nodeSDs=[]
	edgeSens, inDegrees, edgePPVs=[], [], []
	inCoV=[]
	TPsum, TNsum, FPsum, FNsum=0, 0,0,0
	for node in range(0,len(modeler.nodeList)):
		tempSD=0.
		FP, TP, TN, FN=0, 0,0,0
		# simplify rules at the node and find the edge-wise PPV, sens, and SDs
		inCovTemper=[]
		for k in range(len(truthList)):
			inCovtemp=[]
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[k], node, truthList[k])
			start2,end2 = findEnds(model2s[k], node, testList[k])			
			
			# find the shadow and nodes for each model
			truthInEdges=[]
			testInEdges=[]
			for lister in model1s[k].andNodeList[node]:
				tempTup=tuple(lister)
				truthInEdges.append(set(tempTup))
			for lister in model2s[k].andNodeList[node]:
				tempTup=tuple(lister)
				testInEdges.append(set(tempTup))

			# find the values for just this node
			truth= truthList[k][start1:end1]
			test= testList[k][start2:end2]
			# simplify ground truth rule
			for i in range(len(truth)):
				if truth[i]==1:
					for j in range(len(truth)):
						if truth[j]==1 and not i==j:
							if truthInEdges[i].issubset(truthInEdges[j]):
								truth[j]=0
			# simplify GA identified rule
			for i in range(len(test)):
				if test[i]==1:
					for j in range(len(test)):
						if test[j]==1 and not i==j:
							if testInEdges[i].issubset(testInEdges[j]):
								test[j]=0
			testList[k][start2:end2] = test
			truthList[k][start1:end1] = truth
		
			# find SD, PPV, etc....
			truthSet=Set([])
			testSet=Set([])
			baseSet=Set([])
			for i in range(0,len(truth)):
				if truth[i]==1:
					for nodeToAdd in model1s[k].andNodeList[node][i]:
						truthSet.add(nodeToAdd)
						inCovtemp.append(covs[k][node])
				for nodeToAdd in model1s[k].andNodeList[node][i]:
					baseSet.add(nodeToAdd)
			for i in range(0,len(test)):
				if test[i]==1:
					for nodeToAdd in model2s[k].andNodeList[node][i]:
						testSet.add(nodeToAdd)
			SDs[k]=SDs[k]+len(truthSet.difference(testSet))+len(testSet.difference(truthSet))
			tempSD=tempSD+len(truthSet.difference(testSet))+len(testSet.difference(truthSet))
			FP+=1.*len(testSet.difference(truthSet))
			TP+=1.*len(testSet.intersection(truthSet))
			TN+=1.*(len((baseSet.difference(truthSet)).difference(testSet)))
			FN+=1.*len(truthSet.difference(testSet))
			inCovTemper.append(numpy.mean(inCovtemp))
		if (TP+FN)>0:
			sensitivity=1.*TP/(TP+FN)
		else:
			sensitivity=100
		if TP+FP>0:
			PPV=1.*TP/(TP+FP)
		else:
			PPV=100
		nodeSDs.append(tempSD/len(truthList))
		edgeSens.append(sensitivity)
		edgePPVs.append(PPV)
		TPsum+=TP
		TNsum+=TN
		FNsum+=FN
		FPsum+=FP
		inDegrees.append(len(baseSet))
		inCoV.append(numpy.mean(inCovTemper))
	if (TPsum+FNsum)>0:
		edgeSens=1.*TPsum/(TPsum+FNsum)
	else:
		edgeSens=100
	if (FPsum+TPsum)>0:
		edgePPV= 1.*TPsum/(FPsum+TPsum)
	else:
		edgePPV=100



	nodeSens=[] # sensitivity by node
	nodePPV=[] # PPV by node
	nodeRTs=[] # rules true by node
	nodePsens=[]
	nodepPPV=[]
	nodelister=model1s[0].nodeList # gives a node List (should be the same across all trials in a network...)
	sampleRTs=[[] for item in truthList] # Rules True for each trial
	samplePPVs=[[] for item in truthList] # PPV for each trial
	sampleSenss=[[] for item in truthList] # Sens for each trial
	equivRTsens=[[] for item in truthList] # RT sensitivity of equivalents for each trial
	equivSens=[[] for item in truthList] # sensitivity for equivalents for each trial
	equivNodeRTsens=[]
	equivNodeSens=[]

	# iterate over all nodes in the network
	for node in range(len(nodelister)):
		rtTemp=[] # stores rules true for this node across all networks
		ones=[] # stores the number of false negative and rules
		zeros=[] # stores the number of  correct and rules
		negones=[] # stores the number of false positive and rules 
		pOnes=[]
		pZeros=[]
		pNegOnes=[]
		sumindividual=[] # total number true positive and rules
		equivRTsensNode=[]
		equivSensNode=[]
		
		#loop over individuals provided and calculate sens, PPV, rules true
		for i in range(len(truthList)):
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[i], node, truthList[i])
			start2,end2 = findEnds(model2s[i], node, testList[i])			
			# find the values for just this node
			truth= truthList[i][start1:end1]
			test= testList[i][start2:end2]




			# set up empty lists for ands, edges, and the shadow and nodes associated with this node in each model
			
			truthAnds=[]
			testAnds=[]
			truthEdges=set([])
			testEdges=set([])
			inEdges1=[]
			inEdges2=[]
			# add all possible incoming shadow and nodes to a list for each model
			for lister in model1s[i].andNodeList[node]:
				inEdges1.append(set(lister))
			for lister in model2s[i].andNodeList[node]:
				inEdges2.append(set(lister))
			# print(truth)
			# print(model1s[i].andNodeList[node])
			# get the set of all shadow and nodes that are actually used in each rule
			for j in range(len(model1s[i].andNodeList[node])):
				if truth[j]>0:
					truthAnds.append(tuple(model1s[i].andNodeList[node][j]))
					truthEdges.update(set(model1s[i].andNodeList[node][j]))
			for j in range(len(model2s[i].andNodeList[node])):
				if test[j]>0:
					testAnds.append(tuple(model2s[i].andNodeList[node][j]))
					testEdges.update(set(model2s[i].andNodeList[node][j]))
			truthAnd=tuple(truthAnds)
			truthAnd=set(truthAnd)
			testAnd=set(tuple(testAnds))
			
			equivAnds=[]
			# print(equivs[i])
			for test1 in equivs[i][node]:
				truthAnds1=[]
				testAnds1=[]
				truthEdges1=set([])
				testEdges1=set([])
				for j in range(len(model2s[i].andNodeList[node])):
					if test1[j]>0:
						testAnds1.append(tuple(model2s[i].andNodeList[node][j]))
						testEdges1.update(set(model2s[i].andNodeList[node][j]))
				testAnd1=set(tuple(testAnds1))
				equivAnds.append(testAnd1)
			RTequiv=0.
			possibilityOnes=[]
			possibilityZeros=[]
			possibilityZNetones=[]
			for testAnder1 in equivAnds:	
				if (truthAnd==testAnder1):
					RTequiv=1.
				possibilityOnes.append(len(truthAnd.difference(testAnd)))
				possibilityZeros.append(len(truthAnd.intersection(testAnd)))
				possibilityZNetones.append(len(testAnd.difference(truthAnd)))
			maxpossibilityZeros=max(possibilityZeros)
			minpossiblityOnes=min(possibilityOnes)
			minpossibilityNegOnes=min(possibilityZNetones)
			pOnes.append(minpossiblityOnes)
			pZeros.append(maxpossibilityZeros)
			pNegOnes.append(minpossibilityNegOnes)


			equivRTsensNode.append(RTequiv)
			equivRTsens[i].append(RTequiv)
			# calculate true positives, false positives, false negatives, and total slots for this node, trial and save
			onetemp=len(truthAnd.difference(testAnd))
			zerotemp=len(truthAnd.intersection(testAnd))
			negonetemp=len(testAnd.difference(truthAnd))
			sumindtemp=len(truthAnd)
			ones.append(onetemp)
			zeros.append(zerotemp)
			negones.append(negonetemp)
			sumindividual.append(sumindtemp)			
			# add Rules true first sample-wise then node-wise
			if len(inEdges1)>1:
				if (truthAnd==testAnd):
					sampleRTs[i].append(1.)
				else:
					sampleRTs[i].append(0.)
				if (sumindtemp-onetemp+negonetemp)>0:
					samplePPVs[i].append(1.*(sumindtemp-onetemp)/(sumindtemp-onetemp+negonetemp))
				else:
					samplePPVs[i].append(100)
				if (sumindividual[i])>0:
					sampleSenss[i].append(1.*(sumindtemp-onetemp)/(sumindtemp))
				else:
					sampleSenss[i].append(100)
			if (truthAnd==testAnd):
				rtTemp.append(1.)
			else:
				rtTemp.append(0.)

		nodeRTs.append(numpy.mean(rtTemp)) # node-wise Rules true added
		equivNodeRTsens.append(numpy.mean(equivRTsensNode))

		# calculate sensitivity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(1.*numpy.sum(temp)/len(temp))


		# calculate max sensitivity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-pOnes[i])/(sumindividual[i]) for i in range(0,len(pOnes))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			psensitivity=100
		else:
			psensitivity=(1.*numpy.sum(temp)/len(temp))
		nodePsens.append(psensitivity)

		# calculate PPV for the node
		temp=[100 if (sumindividual[i]-ones[i]+negones[i])==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]-ones[i]+negones[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			PPV=100
		else:
			PPV=(1.*numpy.sum(temp)/len(temp))



		# calculate PPV for the node
		temp=[100 if (sumindividual[i]-pOnes[i]+pNegOnes[i])==0 else 1.*(sumindividual[i]-pOnes[i])/(sumindividual[i]-pOnes[i]+pNegOnes[i]) for i in range(0,len(pOnes))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			pPPV=100
		else:
			pPPV=(1.*numpy.sum(temp)/len(temp))
		nodepPPV.append(pPPV)

		# add to list of sensitivity and PPV by 
		nodeSens.append(sensitivity)
		nodePPV.append(PPV)
	sampleEquivRT=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in equivRTsens] # Rules True for each trial
	sampleRT=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in sampleRTs] # Rules True for each trial
	samplePPV=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in samplePPVs] # PPV for each trial
	sampleSens=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in sampleSenss] # Sens for each trial
	return sampleEquivRT, equivNodeRTsens, nodePsens,nodepPPV,sampleSens, samplePPV, nodeSens, nodePPV, sampleRT, nodeRTs, edgeSens, edgePPV, SDs, nodeSDs, len(modeler.nodeList), inDegrees, inCoV

class modelHolder:
	def __init__(self,valueList):
		[ size, nodeList, individualParse, andNodeList, andNodeInvertList, andLenList,nodeList, nodeDict, initValueList]=valueList
		self.size=size
		self.individualParse=individualParse
		self.andNodeList=andNodeList
		self.andNodeInvertList=andNodeInvertList
		self.andLenList=andLenList
		self.nodeList=nodeList
		self.nodeDict=nodeDict
		self.initValueList=initValueList

def analyzeGraph(stem, replicates):
	# upload data for this replicate
	model1s=[]
	model2s=[]
	truthList=[]
	testList=[]
	varLists=[]
	equivalents=[]
	devs=[]
	for i in range(1,replicates+1):
		outputList=pickle.Unpickler(open( stem+str(i)+'_output.pickle', "rb" )).load()
		[truth, test, truthmodel, testmodel, equivs, dev]=outputList
		equivalents.append(equivs)
		model1s.append(modelHolder(truthmodel))
		model2s.append(modelHolder(testmodel))
		truthList.append(truth)
		testList.append(test)
		inputList=pickle.Unpickler(open( stem+str(i)+'_input.pickle', "rb" )).load()
		varList=[]
		for j in range(len(inputList[0])):
			varList.append(numpy.std([inputList[k][j] for k in range(len(inputList))]))
		varLists.append(varList)
	print(len(truthList))
	print(len(model1s))
	return model1s[0],compareIndividualsNodeWise(truthList, testList, model1s,model2s, varLists, equivalents), varLists

def analyzeExperiment(codes, end):
	degree2dictlist, degree3dictlist, SDs, nodeSDS, nodeNumbers,nodeNumlistcut, nodeNumlistuncut, dfdictlist, sensspeclist=[], [], [], [],[],[],[], [], []
	nodewiseList=[]
	nodewiseList3=[]

	nodewiseList2=[]
	# for each gpickle
	for code in codes:
		print(code)		
		nodeLookup={}
		tempsensitivities=[[],[],[],[]]
		tempspecificities=[[],[],[],[]]

		ruleTruthtots=[]
		#FOR EACH RUN WITH THAT GPICKLE
		for noiseNum in range(1,6):
			if noiseNum==1:
				noise=2
			elif noiseNum==2:
				noise=3
			elif noiseNum==3:
				noise=4
			elif noiseNum==4:
				noise=5
			elif noiseNum==5:
				noise=10	
			truthmodel,result, varList=analyzeGraph('pickles/'+code+'_'+str(noiseNum)+'_',10)
			sampleEquivRT, equivNodeRTsens,nodePsens,nodepPPV, tempSens, tempPPV, nodeSens, nodePPV, ruleTruth, nodeRuleTruth, edgeSens, edgePPV, SD, nodeSDs, nodeNum, edgeDegree, inCoV = result
			upstreamPPV, upstreamSens, upstreamRT, upstreamMaxRT=findUpstreamChars('pickles/'+code+'_', nodePPV, nodeSens, nodeRuleTruth, equivNodeRTsens, truthmodel)
			for i in range(len(nodePPV)):
				cov=numpy.mean([varList[k][i] for k in range(10)])
				if not nodePPV[i]==100:
					nodewiseList.append({'equiv Sens': nodePsens[i],'equiv PPV': nodepPPV[i],'equiv RT':equivNodeRTsens[i], 'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'CoV': cov,'in_CoV': inCoV[i], 'SD':nodeSDs[i], 'in_degree':edgeDegree[i], 'upstream Max RT': upstreamMaxRT[i], 'upstream RT': upstreamRT[i], 'noise':noise})
					if edgeDegree[i]==3:
						nodewiseList3.append({'equiv Sens': nodePsens[i],'equiv PPV': nodepPPV[i],'equiv RT':equivNodeRTsens[i],'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'CoV': cov,'in_CoV': inCoV[i], 'SD':nodeSDs[i], 'in_degree':edgeDegree[i], 'upstream Max RT': upstreamMaxRT[i], 'upstream RT': upstreamRT[i], 'noise':noise})
					if edgeDegree[i]==2:
						nodewiseList2.append({'equiv Sens': nodePsens[i],'equiv PPV': nodepPPV[i],'equiv RT':equivNodeRTsens[i],'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'CoV': cov,'in_CoV': inCoV[i], 'SD':nodeSDs[i], 'in_degree':edgeDegree[i], 'upstream Max RT': upstreamMaxRT[i], 'upstream RT': upstreamRT[i], 'noise':noise})
			nodeNumbers.append(nodeNum)
			if len(ruleTruth)>0:
				ruleTruthtots.append(numpy.mean(ruleTruth))
			SDs.append(SD)
			nodeSensDict={}
			nodePPVDict={}
			nodeRTdict={}
			# assign what the var of interest is
			for i in range(len(tempSens)):
				sensspeclist.append({'Node_Num':nodeNum,'Sensitivity':tempSens[i],'Equiv RT Sens':sampleEquivRT[i],'PPV':tempPPV[i], 'Proportion Rules True': ruleTruth[i],'GA': 'no Adapt','Noise':noise})
		print(nodeNumbers[-1])	
	plot_df01box(sensspeclist,'Node_Num', 'Equiv RT Sens','Rules True by Node Number', 'Node Number', 'Maximum Proportion Rules True','Adapt_GA_Equiv_RT_sens'+end+'.png', True, 'Noise')
	plot_df01box(sensspeclist,'Node_Num', 'Proportion Rules True','Rules True by Node Number', 'Node Number', 'Proportion Rules True','Adapt_GA_Rules_True'+end+'.png', True, 'Noise')
	plot_df01box(sensspeclist,'Noise', 'Equiv RT Sens','Rules True by Sample Number', 'Sample Number', 'Maximum Proportion Rules True','Noise_GA_Equiv_RT_sens'+end+'.png', False, 'Noise')
	plot_df01box(sensspeclist,'Noise', 'Proportion Rules True','Rules True by Sample Number', 'Sample Number', 'Proportion Rules True','Noise_GA_Rules_True'+end+'.png', False, 'Noise')

	df=pd.DataFrame(sensspeclist)
	groupedvalues=df.groupby('Node_Num').mean().reset_index()
	print(groupedvalues)
	plot_df01box(sensspeclist,'Node_Num', 'Sensitivity','Sensitivity by Node Number', 'Node Number', 'Sensitivity','Adapt_GA_Sensitivity'+end+'.png', False, 'Noise')
	plot_df01box(sensspeclist,'Node_Num', 'PPV','PPV by Node Number', 'Node Number', 'PPV','Adapt_GA_PPV'+end+'.png', False, 'Noise')

def compareIndividualsNodeWise2(truthList, testList, model1s, model2s):

	modeler=model1s[0]
	SDs=[0. for q in truthList]
	nodeSDs=[]
	edgeSens, inDegrees, edgePPVs=[], [], []
	TPsum, TNsum, FPsum, FNsum=0, 0,0,0
	for node in range(0,len(modeler.nodeList)):
		tempSD=0.
		FP, TP, TN, FN=0, 0,0,0
		# simplify rules at the node and find the edge-wise PPV, sens, and SDs
		for k in range(len(truthList)):
			
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[k], node, truthList[k])
			start2,end2 = findEnds(model2s[k], node, testList[k])			
			
			# find the shadow and nodes for each model
			truthInEdges=[]
			testInEdges=[]
			for lister in model1s[k].andNodeList[node]:
				tempTup=tuple(lister)
				truthInEdges.append(set(tempTup))
			for lister in model2s[k].andNodeList[node]:
				tempTup=tuple(lister)
				testInEdges.append(set(tempTup))

			# find the values for just this node
			truth= truthList[k][start1:end1]
			test= testList[k][start2:end2]
			# simplify ground truth rule
			for i in range(len(truth)):
				if truth[i]==1:
					for j in range(len(truth)):
						if truth[j]==1 and not i==j:
							if truthInEdges[i].issubset(truthInEdges[j]):
								truth[j]=0
			# simplify GA identified rule
			for i in range(len(test)):
				if test[i]==1:
					for j in range(len(test)):
						if test[j]==1 and not i==j:
							if testInEdges[i].issubset(testInEdges[j]):
								test[j]=0
			testList[k][start2:end2] = test
			truthList[k][start1:end1] = truth
		
			# find SD, PPV, etc....
			truthSet=Set([])
			testSet=Set([])
			baseSet=Set([])
			for i in range(0,len(truth)):
				if truth[i]==1:
					for nodeToAdd in model1s[k].andNodeList[node][i]:
						truthSet.add(nodeToAdd)
				for nodeToAdd in model1s[k].andNodeList[node][i]:
					baseSet.add(nodeToAdd)
			for i in range(0,len(test)):
				if test[i]==1:
					for nodeToAdd in model2s[k].andNodeList[node][i]:
						testSet.add(nodeToAdd)
			SDs[k]=SDs[k]+len(truthSet.difference(testSet))+len(testSet.difference(truthSet))
			tempSD=tempSD+len(truthSet.difference(testSet))+len(testSet.difference(truthSet))
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
		nodeSDs.append(tempSD/len(truthList))
		edgeSens.append(sensitivity)
		edgePPVs.append(PPV)
		TPsum+=TP
		TNsum+=TN
		FNsum+=FN
		FPsum+=FP
		inDegrees.append(len(baseSet))
	if (TPsum+FNsum)>0:
		edgeSens=1.*TPsum/(TPsum+FNsum)
	else:
		edgeSens=100
	if (FPsum+TPsum)>0:
		edgePPV= 1.*TPsum/(FPsum+TPsum)
	else:
		edgePPV=100



	nodeSens=[] # sensitivity by node
	nodePPV=[] # PPV by node
	nodeRTs=[] # rules true by node
	nodelister=model1s[0].nodeList # gives a node List (should be the same across all trials in a network...)
	sampleRTs=[[] for item in truthList] # Rules True for each trial
	samplePPVs=[[] for item in truthList] # PPV for each trial
	sampleSenss=[[] for item in truthList] # Sens for each trial
	# iterate over all nodes in the network
	for node in range(len(nodelister)):
		rtTemp=[] # stores rules true for this node across all networks
		ones=[] # stores the number of false negative and rules
		zeros=[] # stores the number of  correct and rules
		negones=[] # stores the number of false positive and rules 
		sumindividual=[] # total number true positive and rules
		
		#loop over individuals provided and calculate sens, PPV, rules true
		for i in range(len(truthList)):
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[k], node, truthList[k])
			start2,end2 = findEnds(model2s[k], node, testList[k])			
			# find the values for just this node
			truth= truthList[i][start1:end1]
			test= testList[i][start2:end2]

			# set up empty lists for ands, edges, and the shadow and nodes associated with this node in each model
			truthAnds=[]
			testAnds=[]
			truthEdges=set([])
			testEdges=set([])
			inEdges1=[]
			inEdges2=[]
			# add all possible incoming shadow and nodes to a list for each model
			for lister in model1s[i].andNodeList[node]:
				inEdges1.append(set(lister))
			for lister in model2s[i].andNodeList[node]:
				inEdges2.append(set(lister))
			# print(truth)
			# print(model1s[i].andNodeList[node])
			# get the set of all shadow and nodes that are actually used in each rule
			for j in range(len(model1s[i].andNodeList[node])):
				if truth[j]>0:
					truthAnds.append(tuple(model1s[i].andNodeList[node][j]))
					truthEdges.update(set(model1s[i].andNodeList[node][j]))
			for j in range(len(model2s[i].andNodeList[node])):
				if test[j]>0:
					testAnds.append(tuple(model1s[i].andNodeList[node][j]))
					testEdges.update(set(model1s[i].andNodeList[node][j]))
			truthAnd=tuple(truthAnds)
			truthAnd=set(truthAnd)
			testAnd=set(tuple(testAnds))
			
			# calculate true positives, false positives, false negatives, and total slots for this node, trial and save
			onetemp=len(truthAnd.difference(testAnd))
			zerotemp=len(truthAnd.intersection(testAnd))
			negonetemp=len(testAnd.difference(truthAnd))
			sumindtemp=len(truthAnd)
			ones.append(onetemp)
			zeros.append(zerotemp)
			negones.append(negonetemp)
			sumindividual.append(sumindtemp)			

			# add Rules true first sample-wise then node-wise
			if len(inEdges1)>1:
				if (truthAnd==testAnd):
					sampleRTs[i].append(1.)
				else:
					sampleRTs[i].append(0.)
				if (sumindtemp-onetemp+negonetemp)>0:
					samplePPVs[i].append(1.*(sumindtemp-onetemp)/(sumindtemp-onetemp+negonetemp))
				else:
					samplePPVs[i].append(100)
				if (sumindividual[i])>0:
					sampleSenss[i].append(1.*(sumindtemp-onetemp)/(sumindtemp))
				else:
					sampleSenss[i].append(100)
			if (truthAnd==testAnd):
				rtTemp.append(1.)
			else:
				rtTemp.append(0.)

		nodeRTs.append(rtTemp) # node-wise Rules true added


		# calculate sensitivity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(1.*numpy.sum(temp)/len(temp))

		# calculate PPV for the node
		temp=[100 if (sumindividual[i]-ones[i]+negones[i])==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]-ones[i]+negones[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			PPV=100
		else:
			PPV=(1.*numpy.sum(temp)/len(temp))
		# add to list of sensitivity and PPV by 
		nodeSens.append(sensitivity)
		nodePPV.append(PPV)

	sampleRT=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in sampleRTs] # Rules True for each trial
	samplePPV=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in samplePPVs] # PPV for each trial
	sampleSens=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in sampleSenss] # Sens for each trial
	return sampleSens, samplePPV, nodeSens, nodePPV, sampleRT, nodeRTs, edgeSens, edgePPV, SDs, nodeSDs, len(modeler.nodeList), inDegrees

if __name__ == '__main__':
	codes=[]
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			codes.append(file[:-8])
	print(codes)

	analyzeExperiment(codes,'')
	fitnesses=[]
	model1s=[]
	model2s=[]
	dictmaker=[]
	dictmaker2=[]
	dictmaker3, dictmaker4, dictmaker5, dictmaker6, dictmaker7, dictmaker8, dictmaker9, dictmaker10, dictmaker11, dictmaker12= [],[],[],[],[],[],[],[],[],[]
	
	dataminer=[]