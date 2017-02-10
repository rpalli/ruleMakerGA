import numpy as numpy
from random import random
def genRandBits(individualLength): #makes a random bitstring
	arr = numpy.random.randint(2, size=(int(individualLength),))
	return list(arr) 

def bitList(n):
    return [1 if digit=='1' else 0 for digit in bin(n)[2:]]

def bit2int(bitlist): #converts bitstring to integer
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out


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

def synthesizeInputs(graph,samples): # generates synthetic random inputs for steady state simulations
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
		addString=addString+writeNode(i,individual, model)
		addString=addString+'\n'
	return addString

def writeBruteNode(currentNode,individual,model):
	padindividual=[0 for x in range(0,model.individualParse[currentNode][0])]
	padindividual.extend(individual)
	return(writeNode(currentNode, padindividual,model))
def writeNode(currentNode,individual, model):
	#write out evaluation instructions in BooleanNet format. This follows the exact same code as fuzzyUpdate, but writes a string instead of actually updating the values of the nodes
	triple=model.individualParse[currentNode]
	inputOrder=model.inputOrderList[currentNode]
	inputOrderInvert=model.inputOrderInvertList[currentNode]
	writenode=''+model.nodeList[currentNode]+'='
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

	