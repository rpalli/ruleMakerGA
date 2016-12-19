#import stuff
from random import random
from random import shuffle
from random import gauss
import numpy as numpy
from math import floor, ceil, log
import operator
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
import networkx as nx
import fuzzyGA_oldTests as tester
import fuzzyNetworkConstructor as constructor
import csv
import itertools as itertool
import copy as copy
import gc as gc
import re as re





def Inv(x, inverter): #inverts if necessary then applies hill fun
	if inverter:
		return (1-x)
	else:
		return (x)

def hillInv(x, inverter,h,p,hillOn): #inverts if necessary then applies hill fun
	if inverter:
		return hill(1-x, h,p,hillOn)
	else:
		return hill(x, h,p,hillOn)


def fuzzyAnd(num1, num2):
	#return num1*num2
	return min(num1,num2)
			
def fuzzyOr(num1, num2):
	#return (num1+num2-(num1*num2))
	return max(num1, num2)

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
		