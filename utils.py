# import necessary python packages
import numpy as numpy
from random import random
import csv as csv

def genRandBits(individualLength): #makes a random bitstring
	arr = numpy.random.randint(2, size=(int(individualLength),))
	return list(arr) 

def bitList(n, x):
	templist=[1 if digit=='1' else 0 for digit in bin(n)[2:]]
	while len(templist)<x:
		templist.insert(0,0)
	return templist

def loadFpkms(filename): #loads data from fpkms tab delimited csv file
	with open(filename) as csvfile:
		data=[]
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			data.append(row)
		return data

def sortFpkms(data): #puts fpkms data into appropriate lists of steady state dictionaries following normalization to largest value (as denominator)
	sss=[]
	for j in range(1,len(data[1])):
		sss.append({})
	for i in range(1,len(data)):
		maxdata=0
		for j in range(1,len(data[i])):
			if float(data[i][j])>maxdata:
				maxdata=float(data[i][j])
		if maxdata==0:
			maxdata=1
		for j in range(0,len(data[i])-1):
			sss[j][str.upper(data[i][0])]=float(data[i][1])/maxdata
	return sss

def synthesizeInputs(graph,samples): # generates synthetic completely random inputs for simulation to steady state
	sss=[]
	for i in range(0,samples):
		sss.append({})
	for node in graph.nodes():
		for i in range(0,samples):
			sss[i][node]=random()
	return sss

def writeModel(individual, model):
	#iterate over nodes to generate a BooleanNet representation for the entire model
	addString=''
	print(model.andLenList)
	print(model.individualParse)
	for i in range(0,len(model.nodeList)):
		addString=addString+writeNode(i,individual[model.individualParse[i]:model.individualParse[i+1]], model)
		addString=addString+'\n'
	return addString

def writeBruteNode(currentNode,individual,model):
	padindividual=[0 for x in range(0,model.individualParse[currentNode][0])]
	padindividual.extend(individual)
	return(writeNode(currentNode, padindividual,model))

def writeNode(currentNode,nodeIndividual, model):
	#write out evaluation instructions in BooleanNet format. 
	# This follows the exact same code as updateNode (for switch=0), but writes a string instead of actually updating the values of the nodes
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	writenode=''+model.nodeList[currentNode]+'=' # set up the initial string to use to write node


	if model.andLenList[currentNode]==0 or sum(nodeIndividual)==0:
		return writenode + ' ' + model.nodeList[currentNode] #if no inputs, maintain value
	elif len(andNodes)==1: 
		#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
		# print(model.andLenList[currentNode])
		# print(nodeIndividual)
		value=''	
		#if only one input, then set to that number
		if andNodeInvertList[0][0]==0:
			value= value+ 'not ' + model.nodeList[andNodes[0][0]]
		else:
			value= value + model.nodeList[andNodes[0][0]]
		return writenode + value 
	else:
		#update nodes with more than one input

		# first deal with case of simple logic without need of linear regression
		orset=[]
		# go through list of possible shadow and nodes to see which ones actually contribute
		for andindex in range(len(nodeIndividual)):
			newval='('
			if nodeIndividual[andindex]==1:
				# if a shadow and contributes, compute its value using its upstream nodes
				if andNodeInvertList[andindex][0]:
					newval=newval+'not '
				newval=newval+model.nodeList[andNodes[andindex][0]]
				for addnode in range(1,len(andNodes[andindex])):
					newval= newval + ' and '
					if andNodeInvertList[andindex][addnode]:
						newval=newval+' not '
					newval=newval+model.nodeList[andNodes[andindex][addnode]]
				orset.append(newval +')')
			#combine the shadow and nodes with or operations
		print(nodeIndividual)
		writenode=writenode + orset.pop()
		for val in orset:
			writenode = writenode + ' or ' + val
		return writenode














	triple=model.individualParse[currentNode]
	inputOrder=model.inputOrderList[currentNode]
	inputOrderInvert=model.inputOrderInvertList[currentNode]
	
	# print(individual[triple[0]:triple[1]])
	# print(triple)
	if model.possibilityNumList[currentNode]>0:
		logicOperatorFlags=individual[triple[1]:triple[2]]
		inputOrder=inputOrder[bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]]
		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%model.possibilityNumList[currentNode]]
		if len(inputOrder)==1:
				if inputOrderInvert[0]:
					writenode=writenode+' not '
				writenode=writenode+model.nodeList[inputOrder[0]]+' '
		else:
			counter =0
			if inputOrderInvert[0]:
				writenode=writenode+' not '
			writenode=writenode+model.nodeList[inputOrder[0]]+' '
			if logicOperatorFlags[0]==0:
				writenode=writenode+' and '
			else:
				writenode=writenode+' or '
			if inputOrderInvert[1]:
				writenode=writenode+' not '
			writenode=writenode+model.nodeList[inputOrder[1]]+' '
			for i in range(2,len(inputOrder)):
				if logicOperatorFlags[i-1]==0:
					writenode=writenode+' and '
				else:
					writenode=writenode+' or '
				if inputOrderInvert[i]:
					writenode=writenode+' not '
				writenode=writenode+model.nodeList[inputOrder[i]]+' '
	else:
		writenode=writenode+' '+ model.nodeList[currentNode]
	return writenode	
