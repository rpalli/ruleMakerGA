#import stuff
from random import random
import numpy as numpy
from math import floor, ceil, log
import operator
from deap import base
from deap import creator
from deap import gp
from deap import tools
from deap import algorithms as algo
import networkx as nx
from scipy.stats import logistic
import re
import urllib2 


def setupGAparams(graph):
	global individualParse #keep a triplet of how to parse out each individual
	global steps
	global nodeList
	global nodeOrderList
	global initValues
	global ss
	global evaluateNodes #list of indices in nodeList to be evaluated from nodelist
	global possibilityNumList
	global h
	global p
	p=.5
	h= 3
	evaluateNodes=[] 
	individualParse=[] 
	nodeOrderList=[] 
	possibilityNumList=[]
	nodeList=graph.nodes()
	nodeDict={}
	for i in range(0,len(nodeList)):
		nodeDict[nodeList[i]]=i
	steps = 10
	counter=0
	initValues=[]
	for i in range(0,len(nodeList)):
		preds=graph.predecessors(nodeList[i])
		for j in range(0,len(preds)):
			preds[j]=nodeDict[preds[j]]
		from itertools import product, repeat
		with_nones = zip(preds, repeat(None))
		possibilities=list(product(*with_nones))
		for j in range(0,len(possibilities)):
			possibilities[j]=list(possibilities[j])
			while None in possibilities[j]:
				possibilities[j].remove(None)
			while [] in possibilities[j]:
				possibilities[j].remove([])
			# print('possible')
			# print possibilities[j]
		while [] in possibilities:
			possibilities.remove([])
		nodeOrderList.append(possibilities)
		possibilityNumList.append(len(possibilities))
		if len(possibilities)==0:
			logNum=0
		else:
			logNum=ceil(log(len(possibilities))/log(2))
		individualParse.append([int(counter),int(counter+logNum),int(counter+logNum+len(possibilities)-1)])
		counter=counter+logNum+len(possibilities)

		if  nodeList[i] in ss.keys():
			initValues.append(ss[nodeList[i]])
			evaluateNodes.append(i)
		else:
			initValues.append(0.5)
	# print(initValues)
	# print(ss)
	# print(len(ss.keys()))
	return counter

def parseKEGGdict(filename):
	dict={}
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]=='D':
			kline=line.split('      ')
			kstring=kline[1]
			kline=kstring.split('  ')
			k=kline[0]
			nameline=line.replace('D      ', 'D')
			nameline=nameline.split('  ')
			namestring=nameline[1]
			nameline=namestring.split(';')
			name=nameline[0]
			dict[k]=name
	return dict

def readKEGG(lines, graph, KEGGdict):
	grouplist=[]
	groups={}
	construct=graph
	stillread=True
	i=0
	dict={}
	while i< len(lines):
		line=lines[i]
		#add nodes
		if "<entry" in line:		
			nameline=line.split('name="')
			namestring=nameline[1]
			nameline=namestring.split('"')
			name=str.lstrip(nameline[0],'ko:')
			name=str.lstrip(name,'gl:')
			name=name.replace(' ','-')
			name=name.replace('ko:','')
			name=name.replace(':','-')
			#print(name)
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			#build a dictionary internal to this file
			idline=line.split('id="')
			idstring=idline[1]
			idline=idstring.split('"')
			id=idline[0]
			if name in KEGGdict.keys():
				dict[id]=KEGGdict[name]
			else:
				dict[id]=name
			#print(dict[id])
			#if the node is a group, we create the data structure to deconvolute the group later
			if type=='group':
				i=i+1
				newgroup=[]
				grouplist.append(id)
				while(not '</entry>' in lines[i]):
					if 'component' in lines[i]:
						newcompline=lines[i].split('id="')
						newcompstring=newcompline[1]
						newcompline=newcompstring.split('"')
						newcomp=newcompline[0]
						newgroup.append(newcomp)
					i=i+1
				groups[id]=newgroup
			else:
				graph.add_node(dict[id],{'name':dict[id],'type':type})
				#print(name)
				#print(id)
				i=i+1
		#add edges from KEGG file
		elif "<relation" in line:
			color='black'
			signal='a'
			subtypes=[]
			
			entryline=line.split('entry1="')
			entrystring=entryline[1]
			entryline=entrystring.split('"')
			entry1=entryline[0]
			
			entryline=line.split('entry2="')
			entrystring=entryline[1]
			entryline=entrystring.split('"')
			entry2=entryline[0]
			

			
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			while(not '</relation>' in lines[i]):
				if 'subtype' in lines[i]:
					nameline1=lines[i].split('name="')
					namestring1=nameline1[1]
					nameline1=namestring1.split('"')
					subtypes.append(nameline1[0])
				i=i+1
			
			
			#color and activation assignament based on the type of interaction in KEGG
			if ('activation' in subtypes) or ('expression' in subtypes):
				color='green'
				signal='a'
			elif 'inhibition' in subtypes:
				color='red'
				signal='i'
			elif ('binding/association' in subtypes) or('compound' in subtypes):
				color='purple'
				signal='a'
			elif 'phosphorylation' in subtypes:
				color='orange'
				signal='a'
			elif 'dephosphorylation' in subtypes:
				color='pink'
				signal='i'
			elif 'indirect effect' in subtypes:
				color='cyan'
				signal='a'
			elif 'dissociation' in subtypes:
				color='yellow'
				signal='i'
			elif 'ubiquitination' in subtypes:
				color='cyan'
				signal='i'
			else:
				print('color not detected. Signal assigned to activation arbitrarily')
				print(subtypes)
			
			#this code deconvolutes the "groups" that KEGG creates 
			#each element of the group gets an edge to or from it instead of the group
			if entry1 in grouplist:
				for current1 in groups[entry1]:
					node1=dict[current1]
					if entry2 in grouplist:
						for current2 in groups[entry2]:
							node2=dict[current2]
							graph.add_edge(node1, node2, color=color, subtype=name, type=type, signal=signal )
					else:
						node2=dict[entry2]
						graph.add_edge(node1, node2, color=color,  subtype=name, type=type, signal=signal )
			elif entry2 in grouplist:
				node1=dict[entry1]
				for current in groups[entry2]:
					if current in grouplist:
						for current1 in groups[current]:
							if current in dict.keys():
								node2=dict[current]
								graph.add_edge(node1,node2, color=color,  subtype=name, type=type, signal=signal )
							else: 
								print('groups on groups on groups')
						#print('groups on groups')
					elif current in dict.keys():
						node2=dict[current]
						graph.add_edge(node1,node2, color=color, subtype=name, type=type, signal=signal )
					else:
						print('not in the keys')
						print(current)
			else:
				node1=dict[entry1]
				node2=dict[entry2]
				
				graph.add_edge(node1,node2, color=color, subtype=name, type=type, signal=signal )
		i=i+1
def parseKEGGdict(filename, aliasDict, dict):
	#reads KEGG dictionary of identifiers between orthology numbers and actual protein names and saves it to a python dictionary
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]=='D':
			kline=line.split('      ')
			kstring=kline[1]
			kline=kstring.split('  ')
			k=kline[0]
			nameline=line.replace('D      ', 'D')
			nameline=nameline.split('  ')
			namestring=nameline[1]
			nameline=namestring.split(';')
			name=nameline[0]
			#print(name)
			if ',' in name:
				nameline=name.split(',')
				name=nameline[0]
				for entry in range(1,len(nameline)):
					#print(nameline[entry])
					aliasDict[nameline[entry].strip()]=name
			dict[k]=name
	return dict


def uploadKEGGfiles(filelist, graph, foldername):
	for file in filelist:
		inputfile = open(foldername+'/'+file, 'r')
		lines = inputfile.readlines()
		readKEGG(lines, graph)

		
def uploadKEGGfolder(foldername, graph):
#uploads an entire folder of KEGG files given by foldername to the netx graph
	filelist= [ f for f in listdir(foldername+'/') if isfile(join(foldername,join('/',f))) ]
	uploadKEGGfiles(filelist, graph, foldername)
	
	
def uploadKEGGcodes(codelist, graph, KEGGdict):
#queries the KEGG for the pathways with the given codes then uploads to KEGG
	for code in codelist:
		url=urllib2.urlopen('http://rest.kegg.jp/get/'+code+'/kgml')
		text=url.readlines()
		readKEGG(text, graph, KEGGdict)
		#print(graph.nodes())

		

def genRandBits():
	global individualLength
	arr = numpy.random.randint(2, size=(individualLength,))
	print(arr)
	return list(arr)
	
def loadFatimaData(filename,tempDict):
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]!='#':
			kline=line.split(' ')
			tempDict[str.upper(kline[0][3:])]=logistic.cdf(float(kline[1]))

	
def hill(x):
	global h
	global p
	return ((1+h**p)*x**p)/(h**p+x**p)		
	
# def hill(x, h, p):
	# return ((1+h**p)*x**p)/(h**p+x**p)
			
def fuzzyAnd(num1, num2):
	return min(num1,num2)
			
def fuzzyOr(num1, num2):
	return max(num1, num2)

	
	
def bit2int(bitlist):
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out


def fuzzyUpdate(currentNode,oldValue,individual):
	global individualParse
	global nodeList
	global nodeOrderList
	global possibilityNumList

	triple=individualParse[currentNode]
	# print(triple[0])
	# print(triple[1])
	# print(individual[triple[0]:triple[1]])
	#print(nodeOrderList)
	nodeOrder=nodeOrderList[currentNode]
	#print(nodeOrderList)
	if possibilityNumList[currentNode]>0:
		#print(bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode])
		logicOperatorFlags=individual[triple[1]:triple[2]]
		#print(logicOperatorFlags)
		#print(nodeOrder)
		nodeOrder=nodeOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]
		if len(nodeOrder)==0:
			value=oldValue[currentNode]
		elif len(nodeOrder)==1:
			#print(nodeOrder[0])
			value=hill(oldValue[nodeOrder[0]])
		else:
			counter =0
			if logicOperatorFlags[0]==0:
				value=fuzzyAnd(hill(oldValue[nodeOrder[0]]),hill(oldValue[nodeOrder[1]]))
			else:
				value=fuzzyAnd(hill(oldValue[nodeOrder[0]]),hill(oldValue[nodeOrder[1]]))
			for i in range(2,len(nodeOrder)):
				if logicOperatorFlags[0]==0:
					value=fuzzyAnd(value,hill(oldValue[nodeOrder[i]]))
				else:
					value=fuzzyAnd(value,hill(oldValue[nodeOrder[i]]))
		return value
	else:
		return oldValue[currentNode]
						
				
def runFuzzySim(individual):
	#do fuzzy simulation. individual is just a long bitstring... need to seperate it out. 
	global invididualParse
	global nodeList
	global initValues
	#list of triples which gives the total list of nodes. 
	counter=0;
	oldValue=list(initValues)
	newValue=list(initValues)
	simData=[]
	simData.append(oldValue)
	global steps
	for step in range(0,steps):
		for i in range(0,len(individualParse)):
			newValue[i]=fuzzyUpdate(i,oldValue,individual)
		oldValue=list(newValue)
		simData.append(newValue)
	array= [0 for x in range(0,len(newValue))]
	for step in range(len(simData)-5,len(simData)):
		#print(simData[step])
		for element in range(0,len(array)):
			array[element]=array[element]+simData[step][element]
	for element in range(0,len(array)):
		array[element]=array[element]/5
	return array
	
	
			
#this is the evaluation function. Need to create a graph then recreate everything with it.
def evaluate(individual):
	global ss
	global evaluateNodes
	#print evaluateNodes
	RME=1.0
	#for i in range(0,len(individual)):
		# print(individual)
		# print(i)
		# print(individual[int(i)])
	boolValues=runFuzzySim(individual)	
	#print(len(boolValues))
	for i in range(0, len(evaluateNodes)):
		RME=RME+(boolValues[evaluateNodes[i]]-ss[nodeList[evaluateNodes[i]]])**2
		#print(RME)
	return RME,

	
	
def mutate(individual):
	global nodeList
	
	return individual,


if __name__ == '__main__':
	filename='inputDataFatima.txt'
	KEGGfileName='ko04060.xml'
		
	#two dicts for the models
	nodeUpdateDict={}
	global ss
	ss={}
	loadFatimaData(filename,ss)
	#print(ss.keys())
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=parseKEGGdict('ko00001.keg', aliasDict, dict)
	
	inputfile = open(KEGGfileName, 'r')
	lines = inputfile.readlines()
	readKEGG(lines, graph, KEGGdict)
	
	# KEGGfileName='ko04060.xml'
	# inputfile = open(KEGGfileName, 'r')
	# lines = inputfile.readlines()
	# readKEGG(lines, graph, KEGGdict)
	
	KEGGfileName='ko04062.xml'
	inputfile = open(KEGGfileName, 'r')
	lines = inputfile.readlines()
	readKEGG(lines, graph, KEGGdict)
	#print(graph.nodes())
	
	
	
	currentfile='IL1b_pathways.txt'
	inputfile = open(currentfile, 'r')
	line = inputfile.read()
	codelist=re.findall('ko\d\d\d\d\d',line)	
	print(codelist)
	uploadKEGGcodes(codelist, graph, KEGGdict)
	for node in graph.nodes():
		if node in graph.successors(node):
			graph.remove_edge(node,node)
	global individualLength
	individualLength=setupGAparams(graph)
	#graph input stuff




	#setup toolbox
	toolbox = base.Toolbox()

	pset = gp.PrimitiveSet("MAIN", arity=1)
	pset.addPrimitive(operator.add, 2)
	pset.addPrimitive(operator.sub, 2)
	pset.addPrimitive(operator.mul, 2)

	#make a fitness minimization function
	creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
	#create a class of individuals that are lists from networkx
	creator.create("Individual", list, fitness=creator.FitnessMin)


	#how to create aliases for your individuals... 
	#toolbox.register("attr_float", random.random)
	#need this alias to create new graphs... but we should just be using the first one.... 

	toolbox.register("genRandomBitString", genRandBits)
	toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.individual)
	
	#ind1=toolbox.individual()
	population=toolbox.population(n=1)

	
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)
	hof = tools.HallOfFame(1, similar=numpy.array_equal)
	
	#finish registering the toolbox functions
	toolbox.register("mate", tools.cxTwoPoint)
	toolbox.register("mutate", mutate)
	toolbox.register("select", tools.selBest)
	toolbox.register("evaluate", evaluate)
	toolbox.register("similar", numpy.array_equal)
	algo.eaSimple(population, toolbox, stats=stats, cxpb=.2, mutpb=.2, ngen=2, verbose=True)
	
	
	
	# graphy=nx.DiGraph()
	# graphy.add_edge(1,2,color='blue',type='typing')	

	# ind2=population[2]
	# child1, child2 = [toolbox.clone(ind) for ind in (ind1, ind2)]
	# tools.cxBlend(child1, child2, 0.5)
	# del child1.fitness.values
	# del child2.fitness.values