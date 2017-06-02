# import necessary python packages
import numpy as numpy
from random import random
import csv as csv

def genRandBits(individualLength): #makes a random bitstring
	arr = numpy.random.randint(2, size=(int(individualLength),))
	return list(arr) 


def bitList(n, x):
	templist=[1 if digit=='1' else 0 for digit in bin(n)[::-1]]
	while len(templist)<x:
		templist.append(0)
	while(len(templist))>x:
		templist.pop()
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
	for i in range(0,len(model.nodeList)):
		addString=addString+writeNode(i,individual[model.individualParse[i]:model.individualParse[i+1]], model)
		addString=addString+'\n'
	return addString

def writeBruteNode(currentNode,individual,model):
	padindividual=[0 for x in range(0,model.individualParse[currentNode][0])]
	padindividual.extend(individual)
	return (writeNode(currentNode, padindividual,model))

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
		writenode=writenode + orset.pop()
		for val in orset:
			writenode = writenode + ' or ' + val
		return writenode


def simpleNetBuild():
	graph = nx.DiGraph()
	graph.add_edge('a1','a', signal='a')	
	graph.add_edge('a2','a', signal='a')
	graph.add_edge('a3','a', signal='a')
	graph.add_edge('a4','a', signal='a')
	graph.add_edge('b1','b', signal='a')
	graph.add_edge('b2','b', signal='a')
	graph.add_edge('b3','b', signal='a')
	graph.add_edge('a','d', signal='a')	
	graph.add_edge('b','d', signal='a')
	graph.add_edge('d','e', signal='a')
	return graph



def LiuNetwork1Builder():
	graph = nx.DiGraph()
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='i')	
	graph.add_edge('f','k', signal='i')
	graph.add_edge('a','c', signal='a')
	graph.add_edge('b','d', signal='a')
	graph.add_edge('c','f', signal='a')	
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')

	return graph