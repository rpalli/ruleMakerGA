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
import itertools as it
import sys
from bs4 import BeautifulSoup
# import pygraphviz
# import matplotlib.patches as mpatches
# from matplotlib.collections import LineCollection
# from matplotlib.colors import ListedColormap, BoundaryNorm
import networkx.drawing.nx_agraph as agraph

#definitions from BioPAX level 3 reference manual (http://www.biopax.org/mediawiki/index.php?title=Specification)
#these biopax classes are iteracted over in the biopax methods
edge_classes = ['Interaction', 'GeneticInteraction', 'MolecularInteraction', 'TemplateReaction', 'Control', 'Catalysis', 'TemplateReactionRegulation', 'Modulation', 'Conversion', 'ComplexAssembly', 'BiochemicalReaction', 'Degradation', 'Transport', 'TransportWithBiochemicalReaction']
node_classes = ['PhysicalEntity', 'Complex', 'Dna', 'DnaRegion', 'Protein', 'Rna', 'RnaRegion', 'SmallMolecule', 'Gene']
graph_classes = ['Pathway']


def drawGraph3(representative, filename, KEGGdict):
	rep=representative.copy()
	dictionary={}
	names=nx.get_node_attributes(representative, 'name')
	for n in rep.nodes():
		if len(names[n].split())==1:
			if names[n] in KEGGdict.keys():
				dictionary[n]=KEGGdict[names[n]]
				#print(rep.node[n]['name'])
		else :
			translated=''
			for word in names[n].split():
				word1=word.lstrip('ko:')
				word1=word1.lstrip('gl:')
				if word1 in KEGGdict.keys():
					translated=translated+KEGGdict[word1]+'-'
				else:
					translated=translated+word1+'-'
			dictionary[n]=translated
	repar= nx.relabel_nodes(rep,dictionary)
	#print(repar.nodes())
	#print(repar.edges())
	B=agraph.to_agraph(repar)        # convert to a graphviz graph\
	B.layout()            # neato layout
	B.draw(filename)       # write postscript in k5.ps with neato layout

def drawGraph2(representative, filename):
	B=agraph.to_agraph(representative)        # convert to a graphviz graph\
	B.layout()            # neato layout
	B.draw(filename)       # write postscript in k5.ps with neato layout

def parseKEGGdict(filename):
	#makes a dictionary to convert ko numbers from KEGG into real gene names
	#this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two.
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
	soup = BeautifulSoup(''.join(lines), 'xml')
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

#initial: lines is an array representing strings of each line of the biopax XML+RDF file
#         graph is an (empty) networkX graph 
#final: graph has been populated with the data from the file in lines 
def read_biopax(lines, graph):
	#constuct a BeautifulSoup4 parser object with the libxml backend
	#input in lines is given as an array, the join function converts to 1 long string
	#On windows this fails, try 'lxml' instead of 'lxml-xml'
	soup = BeautifulSoup(''.join(lines), 'xml')

	#in BioPAX, both physical objects (proteins, genes, etc.) and interactions (activation, 
	#inhibition, etc.) are represented as NODES.  The edges in the system represent membership 
	#or participation in these interactions.  Accordingly, we parse all objects into nodes below.

	#for every biopax class name (defined on lines 29-31)
	for biopax_class in it.chain(node_classes, edge_classes, graph_classes):
		#use XML parser to get all objects of the current class
		biopax_objects = soup.find_all(biopax_class)
		
		#for each instance of this class
		for biopax_object in biopax_objects:
			#get the ID of this object.  They should all be unique
			node_id = biopax_object.get('rdf:ID')

			#this is IMPROPER biopax, but some pathwaycommons files have no ID and specify 
			#only a URL to an external resource. in this case we use the URL as the ID
			if node_id == None:
				node_id = biopax_object.get('rdf:about')

			#every biopax object should have a display name that we can make use of 
			#but sometimes it is not present.  in that case, we re-use the ID as a name
			try:
				node_name = unicode(biopax_object.find('displayName').string)
			except Exception as e:
				node_name = node_id

			#add the node to the graph
			graph.add_node(node_id, {'name': node_name, 'type': biopax_class})

	#at this point, we have all the nodes in the system.  what's left is to connect these nodes
	#together by looking for references to node IDs in the fields of other nodes member/participant
	#fields.  the IDs are frequently prepended with the '#' sign, as in HTML internal referencing.

	#we iterate over every type of class, since most of them can contain subcomponents
	for biopax_class in it.chain(node_classes, edge_classes, graph_classes):
		#use the XML parser to find all objects of the current class
		biopax_objects = soup.find_all(biopax_class)

		#for each object of the current class:
		for biopax_object in biopax_objects:
			#get the ID as above, or use URL of no ID is provided
			node_id = biopax_object.get('rdf:ID')
			if node_id == None:
				node_id = biopax_object.get('rdf:about')

			#if the current object is in the protein, complex, or physicalentity class then
			#if may contain subcomponents (e.g. a complex of multiple proteins)
			if biopax_class in ['Complex', 'Protein', 'PhysicalEntity']:
				#get all labelled components of this object using the XML parser
				components = biopax_object.find_all('component')
				#for each of these components
				for component in components:
					#get the resource.  this corresponds to the ID or URL of the referred object
					resource = component.get('rdf:resource')
					#if the resource is an ID, ir has a # sign in front of it that needs to be removed
					to_node = resource[1:] if resource[0] == '#' else resource
					#we label the edge as contains, since the parent contains the component
					edge_type = 'contains'
					#add the edge from the parent object to the member
					graph.add_edge(node_id, to_node, type=edge_type)

				#there's another type of membership for these classes, and it functions similarly

				#get the memberPhysicalEntity objects
				members = biopax_object.find_all('memberPhysicalEntity')
				#for each member
				for member in members:
					#get the ID/URL
					resource = member.get('rdf:resource')
					#strip off leading '#' sign if it is present
					to_node = resource[1:] if resource[0] == '#' else resource
					#label the edge as 'member'
					edge_type = 'member'
					#add the edge from the parent object to the member
					graph.add_edge(node_id, to_node, type=edge_type)

			#if the current object is in the following list of classes, it represents a form of 
			#regulation.  The controlType field has values like "ACTIVATION" or "INHIBITION" that 
			#define the type of control.  The controller objects are the ones exerting the effect, 
			#and the controlled objects are the ones that are affected
			if biopax_class in ['Control', 'Catalysis', 'TemplateReactionRegulation', 'Modulation']:
				#use the XML parser to get all controller objects
				controllers = biopax_object.find_all('controller')
				#use the XML parser to get all controlled objects
				controlleds = biopax_object.find_all('controlled')

				#for each controller object
				for controller in controllers:
					#get the node ID / URL
					resource = controller.get('rdf:resource')
					#remove the leading # sign if there is one
					from_node = resource[1:] if resource[0] == '#' else resource
					#attempt to get the value of the controlTyoe field
					try:
						edge_type = unicode(biopax_object.find('controlType').string)
					#if it is not present provide a generic "controller" label
					except Exception as e:
						edge_type = unicode('controller')
					#add the edge from the controlled object to the interaction
					graph.add_edge(from_node, node_id, type=edge_type)

				#for each controlled object
				for controlled in controlleds:
					#get the node ID / URL
					resource = controlled.get('rdf:resource')
					#remove leading # sign if it is present
					to_node = resource[1:] if resource[0] == '#' else resource
					#attempt to get the value of the controlType string
					try:
						edge_type = unicode(biopax_object.find('controlType').string)
					#if it is not provided use a generic label
					except Exception as e:
						edge_type = unicode('controlled')
					#add the edge from the interaction to the controlled object
					graph.add_edge(node_id, to_node, type=edge_type)

			#objects in the following list of classes all represent a chemical reaction.  they have 
			#a defined left and right side, with flow conventionally going from left to right. They
			#may also specify RIGHT_TO_LEFT in the reactionDirection field; if so, the roles of left
			#and right are reversed.
			if biopax_class in ['Conversion', 'ComplexAssembly', 'BiochemicalReaction', 'Degradation', 'Transport', 'TransportWithBiochemicalReaction']:
				#get all objects on the left side off this reaction
				lefts = biopax_object.find_all('left')
				#get all objects on the right side of this reaction
				rights = biopax_object.find_all('right')

				#check for RIGHT_TO_LEFT in reaction direction, and reverse right and left if found
				try:
					direction = unicode(biopax_object.find('reactionDirection').string)
					if direction == 'RIGHT_TO_LEFT':
						#reassigns left to right and right to left via tuple packing/unpacking
						(rights, lefts) = (lefts, rights)
				except Exception as e:
					pass
				

				#add edges as appropriate. left edges go from the reactants to the interaction
				#and right edges go from the interaction to the products
				for left in lefts:
					resource = left.get('rdf:resource')
					from_node = resource[1:] if resource[0] == '#' else resource
					edge_type = 'left'
					graph.add_edge(from_node, node_id, type=edge_type)
				for right in rights:
					resource = right.get('rdf:resource')
					to_node = resource[1:] if resource[0] == '#' else resource
					edge_type = 'right'
					graph.add_edge(node_id, to_node, type=edge_type)

			#
			if biopax_class in ['Interaction', 'GeneticInteraction', 'MolecularInteraction', 'TemplateReaction']:
				print "ambiguous direction relationship found"
				participants = biopax_object.find_all('participant')
				for p1 in participants:
					resource = p1.get('rdf:resource')
					to_node = resource[1:] if resource[0] == '#' else resource
					edge_type = 'participant'
					graph.add_edge(node_id, to_node, type=edge_type)
			if biopax_class in ['Pathway']:
				components = biopax_object.find_all('pathwayComponent')
				for component in components:
					resource = component.get('rdf:resource')
					to_node = resource[1:] if resource[0] == '#' else resource
					edge_type = 'component'
					graph.add_edge(node_id, to_node, type=edge_type)

#initial: lines is an array representing strings of each line of the biopax XML+RDF file
#         graph is an (empty) networkX graph 
#final: graph has been populated with the data from the file in lines 
def read_biopax_dev(lines, graph):
	#constuct a BeautifulSoup4 parser object with the libxml backend
	#input in lines is given as an array, the join function converts to 1 long string
	#On windows this fails, try 'lxml' instead of 'lxml-xml'
	soup = BeautifulSoup(''.join(lines), 'xml')
	node_id_to_name = {}

	#in BioPAX, both physical objects (proteins, genes, etc.) and interactions (activation, 
	#inhibition, etc.) are represented as NODES.  The edges in the system represent membership 
	#or participation in these interactions.  Accordingly, we parse all objects into nodes below.

	#for every biopax class name (defined on lines 29-31)
	for biopax_class in it.chain(node_classes, edge_classes, graph_classes):
		#use XML parser to get all objects of the current class
		biopax_objects = soup.find_all(biopax_class)
		
		#for each instance of this class
		for biopax_object in biopax_objects:
			#get the ID of this object.  They should all be unique
			node_id = biopax_object.get('rdf:ID')

			#this is IMPROPER biopax, but some pathwaycommons files have no ID and specify 
			#only a URL to an external resource. in this case we use the URL as the ID
			if node_id == None:
				node_id = biopax_object.get('rdf:about')

			#every biopax object should have a display name that we can make use of 
			#but sometimes it is not present.  in that case, we re-use the ID as a name
			try:
				node_name = unicode(biopax_object.find('displayName').string)
				words = node_name.split(' ')
				if words[-1] == 'mRNA' or words[-1] == 'protein':
					node_name = ' '.join(words[0:-2])
			except Exception as e:
				node_name = node_id

			#add the node to the graph
			graph.add_node(node_name, {'name': node_name, 'type': biopax_class, 'id': node_id})
			node_id_to_name[node_id] = node_name

	#at this point, we have all the nodes in the system.  what's left is to connect these nodes
	#together by looking for references to node IDs in the fields of other nodes member/participant
	#fields.  the IDs are frequently prepended with the '#' sign, as in HTML internal referencing.

	#we iterate over every type of class, since most of them can contain subcomponents
	for biopax_class in it.chain(node_classes, edge_classes, graph_classes):
		#use the XML parser to find all objects of the current class
		biopax_objects = soup.find_all(biopax_class)

		#for each object of the current class:
		for biopax_object in biopax_objects:
			#get the ID as above, or use URL of no ID is provided
			node_id = biopax_object.get('rdf:ID')
			if node_id == None:
				node_id = biopax_object.get('rdf:about')

			node_name = node_id_to_name[node_id]

			#if the current object is in the protein, complex, or physicalentity class then
			#if may contain subcomponents (e.g. a complex of multiple proteins)
			if biopax_class in ['Complex', 'Protein', 'PhysicalEntity']:
				#get all labelled components of this object using the XML parser
				components = biopax_object.find_all('component')
				#for each of these components
				for component in components:
					#get the resource.  this corresponds to the ID or URL of the referred object
					resource = component.get('rdf:resource')
					#if the resource is an ID, ir has a # sign in front of it that needs to be removed
					to_node_id = resource[1:] if resource[0] == '#' else resource
					to_node_name = node_id_to_name[to_node_id]
					#we label the edge as contains, since the parent contains the component
					edge_type = 'contains'
					#add the edge from the parent object to the member
					graph.add_edge(node_name, to_node_name, type=edge_type)

				#there's another type of membership for these classes, and it functions similarly

				#get the memberPhysicalEntity objects
				members = biopax_object.find_all('memberPhysicalEntity')
				#for each member
				for member in members:
					#get the ID/URL
					resource = member.get('rdf:resource')
					#strip off leading '#' sign if it is present
					to_node_id = resource[1:] if resource[0] == '#' else resource
					to_node_name = node_id_to_name[to_node_id]
					#we label the edge as contains, since the parent contains the component
					edge_type = 'member'
					#add the edge from the parent object to the member
					graph.add_edge(node_name, to_node_name, type=edge_type)

			#if the current object is in the following list of classes, it represents a form of 
			#regulation.  The controlType field has values like "ACTIVATION" or "INHIBITION" that 
			#define the type of control.  The controller objects are the ones exerting the effect, 
			#and the controlled objects are the ones that are affected
			if biopax_class in ['Control', 'Catalysis', 'TemplateReactionRegulation', 'Modulation']:
				#use the XML parser to get all controller objects
				controllers = biopax_object.find_all('controller')
				#use the XML parser to get all controlled objects
				controlleds = biopax_object.find_all('controlled')

				#for each controller object
				for controller in controllers:
					#get the node ID / URL
					resource = controller.get('rdf:resource')
					#remove the leading # sign if there is one
					from_node_id = resource[1:] if resource[0] == '#' else resource
					from_node_name = node_id_to_name[from_node_id]

					#attempt to get the value of the controlTyoe field
					try:
						edge_type = unicode(biopax_object.find('controlType').string)
					#if it is not present provide a generic "controller" label
					except Exception as e:
						edge_type = unicode('controller')
					#add the edge from the controlled object to the interaction
					graph.add_edge(from_node_name, node_name, type=edge_type)

				#for each controlled object
				for controlled in controlleds:
					#get the node ID / URL
					resource = controlled.get('rdf:resource')
					#remove leading # sign if it is present
					to_node_id = resource[1:] if resource[0] == '#' else resource
					to_node_name = node_id_to_name[to_node_id]

					#attempt to get the value of the controlType string
					try:
						edge_type = unicode(biopax_object.find('controlType').string)
					#if it is not provided use a generic label
					except Exception as e:
						edge_type = unicode('controlled')
					#add the edge from the interaction to the controlled object
					graph.add_edge(node_name, to_node_name, type=edge_type)

			#objects in the following list of classes all represent a chemical reaction.  they have 
			#a defined left and right side, with flow conventionally going from left to right. They
			#may also specify RIGHT_TO_LEFT in the reactionDirection field; if so, the roles of left
			#and right are reversed.
			if biopax_class in ['Conversion', 'ComplexAssembly', 'BiochemicalReaction', 'Degradation', 'Transport', 'TransportWithBiochemicalReaction']:
				#get all objects on the left side off this reaction
				lefts = biopax_object.find_all('left')
				#get all objects on the right side of this reaction
				rights = biopax_object.find_all('right')

				#check for RIGHT_TO_LEFT in reaction direction, and reverse right and left if found
				try:
					direction = unicode(biopax_object.find('reactionDirection').string)
					if direction == 'RIGHT_TO_LEFT':
						#reassigns left to right and right to left via tuple packing/unpacking
						(rights, lefts) = (lefts, rights)
				except Exception as e:
					pass
				

				#add edges as appropriate. left edges go from the reactants to the interaction
				#and right edges go from the interaction to the products
				for left in lefts:
					resource = left.get('rdf:resource')
					from_node_id = resource[1:] if resource[0] == '#' else resource
					from_node_name = node_id_to_name[from_node_id]

					edge_type = 'left'
					graph.add_edge(from_node_name, node_name, type=edge_type)
				for right in rights:
					resource = right.get('rdf:resource')
					to_node_id = resource[1:] if resource[0] == '#' else resource
					to_node_name = node_id_to_name[to_node_id]

					edge_type = 'right'
					graph.add_edge(node_name, to_node_name, type=edge_type)

			#
			if biopax_class in ['Interaction', 'GeneticInteraction', 'MolecularInteraction', 'TemplateReaction']:
				print "ambiguous direction relationship found"
				participants = biopax_object.find_all('participant')
				for p1 in participants:
					resource = p1.get('rdf:resource')
					to_node_id = resource[1:] if resource[0] == '#' else resource
					to_node_name = node_id_to_name[to_node_id]

					edge_type = 'participant'
					graph.add_edge(node_name, to_node_name, type=edge_type)
			if biopax_class in ['Pathway']:
				components = biopax_object.find_all('pathwayComponent')
				for component in components:
					resource = component.get('rdf:resource')
					to_node_id = resource[1:] if resource[0] == '#' else resource
					to_node_name = node_id_to_name[to_node_id]

					edge_type = 'component'
					graph.add_edge(node_name, to_node_name, type=edge_type)


def simplify_biopax_graph(graph):
	protein_names_and_nodes = {}

	for node in graph.nodes(data=True):
		node_class = node[1]['type']
		if node_class in edge_classes:
			predecessors = graph.predecessors(node[0])
			successors = graph.successors(node[0])
			for p in predecessors:
				for s in successors:
					if not (graph.has_edge(p,s) or graph.has_edge(s,p)) and p != s:
						graph.add_edge(p, s, node[1])
			graph.remove_node(node[0])

		if node_class in ['Protein', 'Complex']:
			protein_names_and_nodes[node[1]['name']] = node[0]

	# remove nodes we don't care about
	for node in graph.nodes(data=True):
		node_class = node[1]['type']
		if node_class in ['SmallMolecule', 'BiochemicalReaction', 'Pathway']:
			predecessors = graph.predecessors(node[0])
			successors = graph.successors(node[0])
			for p in predecessors:
				for s in successors:
					if not (graph.has_edge(p,s) or graph.has_edge(s,p)) and p != s:
						graph.add_edge(p, s, node[1])
			graph.remove_node(node[0])

		#to be looked at in the future
		if node_class in ['Gene', 'Rna', 'Dna']:
			gene_name = node[1]['name']
			if gene_name in protein_names_and_nodes.keys():
				p = protein_names_and_nodes[gene_name]
				successors = graph.successors(node[0])
				for s in successors:
					if not (graph.has_edge(p,s) or graph.has_edge(s,p)) and p != s:
						graph.add_edge(p, s)
				graph.remove_node(node[0])

		if node_class in ['Complex']:
			components = graph.successors(node[0])
			for c in components:
				all_nodes = graph.nodes(data=True)
				for n in all_nodes:
					if n[0] == c:
						c_data = n[1]
				if c_data['type'] == 'Protein':
					p = node[0]
					successors = graph.successors(c)
					for s in successors:
						if not (graph.has_edge(p,s) or graph.has_edge(s,p)) and p != s:
							graph.add_edge(p, s)
					graph.remove_node(c)




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

# PC = Pathway Commons
def download_PC_codes(codelist, graph):
	for code in codelist:
		url = urllib2.urlopen('http://www.pathwaycommons.org/pc2/graph?source='+code+'&kind=neighborhood')
		text=url.readlines()
		read_biopax(text, graph)
		#uncomment below and comment above to enable beta parser with names as indices
		# read_biopax_dev(text, graph)
		print(code)
	simplify_biopax_graph(graph)

# def download_PC_codes(codelist, graph):
# 	for code in codelist:
# 		url = urllib2.urlopen('http://www.pathwaycommons.org/pc2/graph?source='+code+'&kind=neighborhood')
# 		text=url.readlines()
# 		read_biopax_dev(text, graph)
# 		simplify_biopax_graph(graph)
# 		print(code)

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
