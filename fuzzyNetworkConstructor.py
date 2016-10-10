#import stuff
from random import random
from random import shuffle
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
import fuzzyNetworkConstructor as constructor
import csv 
import sys
from bs4 import BeautifulSoup
import itertools as it

## Rohith's OLD version, kept for history
# def parseKEGGdict(filename):
# 	#makes a dictionary to convert ko numbers from KEGG into real gene names
# 	#this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two.
# 	dict={}
# 	inputfile = open(filename, 'r')
# 	lines = inputfile.readlines()
# 	for line in lines:
# 		if line[0]=='D':
# 			kline=line.split('      ')
# 			kstring=kline[1]
# 			kline=kstring.split('  ')
# 			k=kline[0]
# 			nameline=line.replace('D      ', 'D')
# 			nameline=nameline.split('  ')
# 			namestring=nameline[1]
# 			nameline=namestring.split(';')
# 			name=nameline[0]
# 			dict[k]=name
# 	return dict

def parseKEGGdict(filename, aliasDict, dict):
	#reads KEGG dictionary of identifiers between orthology numbers and actual protein names and saves it to a python dictionary
	#extends above method to also keep track of gene names that are identical so we can recover those from input data as well
	#again, read in lines from the file, save ko number and gene name, identify them in the dictionary. 
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
			if ',' in name:
				nameline=name.split(',')
				name=nameline[0]
				for entry in range(1,len(nameline)):
					aliasDict[nameline[entry].strip()]=name
			dict[k]=name
	return dict

def readKEGG(lines, graph, KEGGdict):
	#the network contained in a KEGG file to the graph. 
	grouplist=[] #sometimes the network stores a group of things all activated by a single signal as a "group". We need to rewire the arrows to go to each individual component of the group so we save a list of these groups and a dictionary between id numbers and lists of elements that are part of that group
	groups={}
	i=0 #i controls our movement through the lines of code. 
	dict={} #identifies internal id numbers with names (much like the code earlier does for groups to id numbers of elements of the group)
	while i< len(lines):
		line=lines[i]
		# print line
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
			
			# print(name)
			
			typeline=line.split('type="')
			typestring=typeline[1]
			typeline=typestring.split('"')
			type=typeline[0]
			
			#build a dictionary internal to this file
			idline=line.split('id="')
			idstring=idline[1]
			idline=idstring.split('"')
			id=idline[0]
			
			
			namelist=[]
			if '-' in name: 
				namelines=name.split('-')
				for x in namelines:
					if x in KEGGdict.keys():
						namelist.append(KEGGdict[x])
					else:
						namelist.append(x)
				name=namelist[0]+'-'
				for x in range(1,len(namelist)):
					name=name+namelist[x]+'-'
				name=name[:-1:]
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
				signal='a'
			
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

def readKEGGnew(lines, graph, KEGGdict):
	#read all lines into a bs4 object using libXML parser
	soup = BeautifulSoup(''.join(lines), 'lxml-xml')
	groups = {} # store group IDs and list of sub-ids
	id_to_name = {} # map id numbers to names

	for entry in soup.find_all('entry'):
		entry_name = entry['name'].split(':')[1] if ':' in entry['name'] else entry['name']
		entry_name = KEGGdict[entry_name] if entry_name in KEGGdict.keys() else entry_name
		entry_type = entry['type']
		entry_id = entry['id']

		id_to_name[entry_id] = entry_name

		if entry_type == 'group':
			group_ids = []
			for component in entry.find_all('component'):
				group_ids.append(component['id'])
			groups[entry_id] = group_ids
		else:
			graph.add_node(entry_name, {'name': entry_name, 'type': entry_type})

	for relation in soup.find_all('relation'):
		(color, signal) = ('black', 'a')

		relation_entry1 = relation['entry1']
		relation_entry2 = relation['entry2']
		relation_type = relation['type']

		subtypes = []

		for subtype in relation.find_all('subtype'):
			subtypes.append(subtype['name'])

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
			signal='a'

		#given: a node ID that may be a group
		#returns: a list that contains all group IDs deconvoluted
		def expand_groups(node_id):
			node_list = []
			if node_id in groups.keys():
				for component_id in groups[node_id]:
					node_list.extend(expand_groups(component_id))
			else:
				node_list.extend([node_id])
			return node_list

		entry1_list = expand_groups(relation_entry1)
		entry2_list = expand_groups(relation_entry2)

		for (entry1, entry2) in it.product(entry1_list, entry2_list):
			node1 = id_to_name[entry1]
			node2 = id_to_name[entry2]
			graph.add_edge(node1,node2, color=color, subtype='/'.join(subtypes), type=relation_type, signal=signal)


def uploadKEGGfiles(filelist, graph, foldername, KEGGdict):
	#upload KEGG files from a particular folder.
	#just provide the folder, file names as a list, and the graph you want to put things in.
	#iteratively calls teh readKEGG method
	for file in filelist:
		inputfile = open(foldername+'/'+file, 'r')
		lines = inputfile.readlines()
		readKEGG(lines, graph, KEGGdict)

		
def uploadKEGGfolder(foldername, graph, KEGGdict):
#uploads an entire folder of KEGG files given by foldername to the netx graph
#just picks up names from the folder and passes them to the uploadKEGGfiles method
#be sure not to pass in a folder that has files other than the KEGG files in it
	filelist= [ f for f in listdir(foldername+'/') if isfile(join(foldername,join('/',f))) ]
	uploadKEGGfiles(filelist, graph, foldername, KEGGdict)
	
	
def uploadKEGGcodes(codelist, graph, KEGGdict):
#queries the KEGG for the pathways with the given codes then uploads to graph. Need to provide the KEGGdict so that we can name the nodes with gene names rather than KO numbers
	for code in codelist:
		url=urllib2.urlopen('http://rest.kegg.jp/get/'+code+'/kgml')
		text=url.readlines()
		readKEGG(text, graph, KEGGdict)
		print(code)
		#print(graph.nodes())

if __name__ == '__main__':
	graph = nx.DiGraph()
	dict={}
	aliasDict={}
	KEGGdict=parseKEGGdict('ko00001.keg', aliasDict, dict)
	
	if len(sys.argv) > 1: #if we have supplied a command-line argument
		KEGGfileName=sys.argv[1]
		inputfile = open(KEGGfileName, 'r')
		lines = inputfile.readlines()
		readKEGG(lines, graph, KEGGdict)
		print(graph.nodes())
	
	
	
	currentfile='IL1b_pathways.txt'
	inputfile = open(currentfile, 'r')
	line = inputfile.read()
	
	codelist=re.findall('ko\d\d\d\d\d',line)	
	uploadKEGGcodes(codelist, graph, KEGGdict)
	for node in graph.nodes():
		if node in graph.successors(node):
			graph.remove_edge(node,node)
	global individualLength

	#this line breaks the code
	#individualLength=setupGAparams(graph)

	#graph input stuff

