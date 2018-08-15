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
# import other parts of our code
import utils as utils


def plot_df01swarm(df, xval, yval, title, xlabeler, ylabeler, plotname):

	sns.set(font_scale=2) 
	sns.set_style("whitegrid")
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
def plot_box(df, xval, yval, title, xlabeler, ylabeler, plotname, limits):

	sns.set(font_scale=2) 
	sns.set_style("whitegrid")
	ax=sns.boxplot(data=df,x=xval, y=yval, color='black')
	# ax.set_xticklabels(nodeNums, color='white')
	# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	# plt.title(title, color='black')
	if limits:
		plt.ylim([0,16])
	plt.tight_layout()
	# sns.despine()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()

def plot_scatter(dictlist, xval, yval, title, xlabeler, ylabeler, plotname):

	df=pd.DataFrame(dictlist)
	sns.set(font_scale=2) 
	sns.set_style("whitegrid")
	ax=sns.lmplot(data=df,x=xval, y=yval, fit_reg=True, )
	# ax.set_xticklabels(nodeNums, color='white')
	# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	# plt.title(title, color='black')
	plt.tight_layout()
	# sns.despine()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()

def plot_scatter2(dictlist, xval, yval, title, xlabeler, ylabeler, hue, plotname):

	df=pd.DataFrame(dictlist)
	sns.set(font_scale=2) 
	sns.set_style("whitegrid")
	ax=sns.lmplot(data=df,x=xval, y=yval,hue=hue, fit_reg=False)
	# ax.set_xticklabels(nodeNums, color='white')
	# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')
	plt.ylabel(ylabeler)
	plt.xlabel(xlabeler)
	# plt.title(title, color='black')
	# sns.despine()
	plt.savefig(plotname, bbox_inches='tight', transparent=True)
	plt.clf()

def plot_df01box(dictlist, xval, yval, title, xlabeler, ylabeler, plotname, hueBool, huer):
	df=pd.DataFrame(dictlist)
	#sns.set_style("dark")
	sns.set_palette(sns.color_palette("Greys", 6))
	sns.set(font_scale=2) 
	sns.set_style("ticks")
	if not hueBool:
		ax=sns.barplot(data=df,x=xval, y=yval, color='black', ci=68.26, capsize=.07,errcolor='black')
	else:
		ax=sns.barplot(data=df,x=xval, y=yval, color='gray', hue=huer, ci=68.26, capsize=.07)
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


def plot_df01box2(df, xval, yval, title, xlabeler, ylabeler, plotname, hueBool, huer):
	#sns.set_style("dark")
	sns.set_palette(sns.color_palette("Greys", 6))
	sns.set(font_scale=2) 
	sns.set_style("ticks")
	if not hueBool:
		ax=sns.barplot(data=df,x=xval, y=yval, color='black', ci=68.26, capsize=.07,errcolor='black')
	else:
		ax=sns.barplot(data=df,x=xval, y=yval, color='gray', hue=huer, ci=68.26, capsize=.07)
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
		endy=len(indiv)
	else:
		endy=model.individualParse[node+1]
	starty=model.individualParse[node]
	return starty, endy

def compareIndividualsNodeWise(truthList, testList, model1s, model2s):
	modeler=model1s[0]
	SDs=[0. for q in truthList]
	nodeSDs=[]
	edgeSens, inDegrees, edgePPVs=[], [], []
	TPsum, TNsum, FPsum, FNsum=0,0,0,0
	for node in range(0,len(modeler.nodeList)):
		tempSD=0.
		FP, TP, TN, FN=0,0,0,0
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
		# print('node'+str(node))
		#loop over individuals provided and calculate sens, PPV, rules true
		for i in range(len(truthList)):
			# print('i')
			# print(i)
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[i], node, truthList[i])
			start2,end2 = findEnds(model2s[i], node, testList[i])			
			# find the values for just this node
			truth= truthList[i][start1:end1]
			test= testList[i][start2:end2]
			# print(len(truth))
			# print('startend')
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
			# print(test)
			# print(model2s[i].andNodeList[node])
			# get the set of all shadow and nodes that are actually used in each rule
			for j in range(len(model1s[i].andNodeList[node])):
				if truth[j]>0:
					truthAnds.append(tuple(model1s[i].andNodeList[node][j]))
					truthEdges.update(set(model1s[i].andNodeList[node][j]))
			for j in range(len(model2s[i].andNodeList[node])):
				# print(test)
				# print(model2s[i].andNodeList[node])
				# print(j)
				if test[j]>0:
					testAnds.append(tuple(model2s[i].andNodeList[node][j]))
					testEdges.update(set(model2s[i].andNodeList[node][j]))
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

def analyzeGraph(stem, replicates,gen1):
	# upload data for this replicate
	model1s=[]
	model2s=[]
	truthList=[]
	testList=[]
	varLists=[]
	print(gen1)
	for i in range(1,replicates+1):
		outputList=pickle.Unpickler(open( stem+str(gen1)+"_"+str(i)+'_output.pickle', "rb" )).load()
		[truthListtemp, testListtemp, truthmodel, testmodel]=outputList
		model1s.append(modelHolder(truthmodel))
		model2s.append(modelHolder(testmodel))
		if (len(testListtemp[0])!= model2s[len(model2s)-1].individualParse[-1]):
			print(len(testListtemp[0]))
			print(model2s[len(model2s)-1].individualParse[-1])
		truthList.extend(truthListtemp)
		testList.extend(testListtemp)
	return compareIndividualsNodeWise(truthList, testList, model1s,model2s)

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
		for gen2 in [0,2,5,10,15,20,40]:
			result=analyzeGraph('pickles/'+code+'_',25, gen2)
			tempSens, tempPPV, nodeSens, nodePPV, ruleTruth, nodeRuleTruth, edgeSens, edgePPV, SD, nodeSDs, nodeNum, edgeDegree = result
			# for i in range(len(nodePPV)):
			# 	if not nodePPV[i]==100:
			# 		nodewiseList.append({'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i], 'SD':nodeSDs[i], 'in_degree':edgeDegree[i]})
			# 		if edgeDegree[i]==3:
			# 			nodewiseList3.append({'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i], 'SD':nodeSDs[i], 'in_degree':edgeDegree[i]})
			# 		if edgeDegree[i]==2:
			# 			nodewiseList2.append({'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'SD':nodeSDs[i], 'in_degree':edgeDegree[i]})
			nodeNumbers.append(nodeNum)
			if len(ruleTruth)>0:
				ruleTruthtots.append(numpy.mean(ruleTruth))
			SDs.append(SD)
			nodeSensDict={}
			nodePPVDict={}
			nodeRTdict={}
			# assign what the var of interest is
			for i in range(len(tempSens)):
				sensspeclist.append({'Gen 2':gen2/10.,'Node_Num':nodeNum,'Sensitivity':tempSens[i],'PPV':tempPPV[i], 'Proportion Rules True': ruleTruth[i], 'SD': SD[i]})
	
	df=pd.DataFrame(sensspeclist)

	print(df)
	plot_box(df, 'Gen 2', 'SD', 'SD by noise', 'noise', 'Structural Distance', 'SD_box1.png',True)
	plot_box(df, 'Gen 2', 'SD', 'SD by noise', 'noise', 'Structural Distance', 'SD_box2.png',False)


	# plot_df01box(sensspeclist,'Node_Num', 'Proportion Rules True','Rules True by Node Number', 'Node Number', 'Proportion Rules True','Adapt_GA_Rules_True'+end+'.png', True, 'Gen 2')
	# df=pd.DataFrame(sensspeclist)
	# groupedvalues=df.groupby('Node_Num').mean().reset_index()
	# gb = df.groupby(['Node_Num'])
	# for k, gp in gb:
	# 	plot_df01box2(gp,'Gen 1', 'Proportion Rules True','Rules True by Initial GA generations', 'Gen 1', 'Proportion Rules True','GA_Rules_True'+end+'_'+str(k)+'.png', False, 'Gen 2')
	# 	plot_df01box2(gp,'Gen 1', 'Sensitivity','Sensitivity by Initial GA generations', 'Gen 1', 'Sensitivity','GA_Sensitivity'+end+'_'+str(k)+'.png', False, 'Gen 2')
	# 	plot_df01box2(gp,'Gen 1', 'PPV','PPV by Initial GA generations', 'Gen 1', 'PPV','GA_PPV'+'_'+str(k)+end+'.png', False, 'Gen 2')


if __name__ == '__main__':
	codes=[]
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			codes.append(file[:-8])
	print(codes)

	analyzeExperiment(codes,'')