
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

def loadFatimaData(filename,tempDict): #reads data from differential expression file
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]!='#':
			kline=line.split(' ')
			tempDict[str.upper(kline[0][3:])]=logistic.cdf(float(kline[1]))

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


class paramClass:
	def __init__(self):     
		self.simSteps=100 # number of steps each individual is run when evaluating
		self.generations=50 # generations to run
		self.popSize=100 #size of population
		self.mu= 100#individuals selected
		self.lambd= 200#children produced
		self.bitFlipProb=.1 # prob of flipping bits inside mutation
		self.crossoverProb=.2 # prob of crossing over a particular parent
		self.mutationProb=.2 # prob of mutating a particular parent
		self.async=False # run in asynchronous mode
		self.iters=100 #number of simulations to try in asynchronous mode
		self.complexPenalty=False #penalize models with increased complexity
		self.genSteps=10000 # steps to find steady state with fake data

def returnSimParams():
	simSteps=100 # number of steps each individual is run when evaluating
	generations=50 # generations to run
	popSize=100 #size of population
	mu= 100#individuals selected
	lambd= 200#children produced
	bitFlipProb=.1 # prob of flipping bits inside mutation
	crossoverProb=.2 # prob of crossing over a particular parent
	mutationProb=.2 # prob of mutating a particular parent
	async=False # run in asynchronous mode
	iters=100 #number of simulations to try in asynchronous mode
	complexPenalty=False #penalize models with increased complexity
	genSteps=10000 # steps to find steady state with fake data
	return async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps

def synthesizeInputs(graph,samples): # generates synthetic random inputs for steady state simulations
	sss=[]
	for i in range(0,samples):
		sss.append({})
	for node in graph.nodes():
		for i in range(0,samples):
			sss[i][node]=random()
	return sss

def writeFuzzyNode(currentNode,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	#write out evaluation instructions in BooleanNet format. This follows the exact same code as fuzzyUpdate, but writes a string instead of actually updating the values of the nodes
	triple=individualParse[currentNode]
	inputOrder=inputOrderList[currentNode]
	inputOrderInvert=inputOrderInvertList[currentNode]
	writenode=''+nodeList[currentNode]+'='
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=individual[triple[1]:triple[2]]
		inputOrder=inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		if len(inputOrder)==1:
			#print(inputOrder[0])
			if individual[triple[0]]==1:
				if inputOrderInvert[0]:
					writenode=writenode+' not '
				writenode=writenode+nodeList[inputOrder[0]]+' '
			else:
				writenode=writenode+' '+ nodeList[currentNode]
		else:
			counter =0
			if inputOrderInvert[0]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[inputOrder[0]]+' '
			if logicOperatorFlags[0]==0:
				writenode=writenode+' and '
			else:
				writenode=writenode+' or '
			if inputOrderInvert[1]:
				writenode=writenode+' not '
			writenode=writenode+nodeList[inputOrder[1]]+' '
			for i in range(2,len(inputOrder)):
				if logicOperatorFlags[i-1]==0:
					writenode=writenode+' and '
				else:
					writenode=writenode+' or '
				if inputOrderInvert[i]:
					writenode=writenode+' not '
				writenode=writenode+nodeList[inputOrder[i]]+' '
	else:
		writenode=writenode+' '+ nodeList[currentNode]
	return writenode					