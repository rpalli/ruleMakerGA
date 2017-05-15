#import external code
import operator
import networkx as nx
import re
import urllib2 
import csv 
import itertools as it
import sys
from bs4 import BeautifulSoup
from random import randint 
# import our code
import simulation as sim
import utils as utils

#definitions from BioPAX level 3 reference manual (http://www.biopax.org/mediawiki/index.php?title=Specification)
#these biopax classes are iteracted over in the biopax methods
edge_classes = ['Interaction', 'GeneticInteraction', 'MolecularInteraction', 'TemplateReaction', 'Control', 'Catalysis', 'TemplateReactionRegulation', 'Modulation', 'Conversion', 'ComplexAssembly', 'BiochemicalReaction', 'Degradation', 'Transport', 'TransportWithBiochemicalReaction']
node_classes = ['PhysicalEntity', 'Complex', 'Dna', 'DnaRegion', 'Protein', 'Rna', 'RnaRegion', 'SmallMolecule', 'Gene']
graph_classes = ['Pathway']

def rewireNetwork(graph):
	rewired=graph.copy()
	for i in range(3*len(rewired.edges())):
		a= randint(0,len(rewired.edges())-1)
		b=a
		while a==b:
			b= randint(0,len(rewired.edges())-1)
		edge1=rewired.edges()[a]
		edge2=rewired.edges()[b]
		activity1=rewired.get_edge_data(edge1[0],edge1[1])['signal']
		activity2=rewired.get_edge_data(edge2[0],edge2[1])['signal']
		if (edge2[0]!=edge1[1] and edge1[0]!=edge2[1]):
			rewired.remove_edge(edge1[0],edge1[1])
			rewired.remove_edge(edge2[0],edge2[1])
			rewired.add_edge(edge2[0],edge1[1],signal=activity1)
			rewired.add_edge(edge1[0],edge2[1],signal=activity2)
	return rewired

def simplifyNetwork(graph, ss):
	#network simplification algorithm. 
	# # 1. remove nodes which are neither part of input nor have input to them...
	# # 2. remove straigth paths. 
	# # 3. 


	#collapse complexes of nodes already in the list
	for node in graph.nodes():
		predlist=graph.predecessors(node)
		for pred in predlist:
			if '-' in pred:
				# # print(pred)
				genes=pred.split('-')
				flag=True
				for gene in genes:
					if not gene in predlist:
						flag=False
				if flag:
					graph.remove_edge(pred,node)
	flag=True
	
	# remove nodes with no predecessors or value in given steady state data
	while(flag):
		flag=False	
		# # print(len(graph.nodes()))
		newNodes = [x for x in graph.nodes() if not (len(graph.predecessors(x))==0 and (not (x in ss.keys())))]
		if(len(newNodes)<len(graph.nodes())):
			flag=True
		graph=graph.subgraph(newNodes)
		
		#print(len(graph.nodes()))
			
	#  collapse straight lines
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1 and (len(graph.successors(x))==1))]
	for rm in removeNodeList:
		before=graph.predecessors(rm)[0]
		after=graph.successors(rm)[0]
		edge1=graph.get_edge_data(before,rm)['signal']
		edge2=graph.get_edge_data(rm,after)['signal']
		inhCount=0
		if edge1=='i':
			inhCount=inhCount+1
		if edge2=='i':
			inhCount=inhCount+1
		if inhCount==1:
			graph.add_edge(before,after,signal='i')
		else:
			graph.add_edge(before,after,signal='a')
		graph.remove_node(rm)
	flag=True

	#rewire nodes that have only one upstream node
	#print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1) ]
	for rm in removeNodeList:
		if len(graph.predecessors(rm))==1:
			for after in graph.successors(rm):
				before=graph.predecessors(rm)[0]
				edge1=graph.get_edge_data(before,rm)['signal']
				edge2=graph.get_edge_data(rm,after)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(before,after,signal='i')
				else:
					graph.add_edge(before,after,signal='a')
			graph.remove_node(rm)
	flag=True

	#rewire nodes that have only one downstream node
	#print(len(graph.nodes()))
	removeNodeList= [x for x in graph.nodes() if (len(graph.successors(x))==1) ]
	for rm in removeNodeList:
		if len(graph.successors(x))==1:
			for start in graph.predecessors(rm):
				finish=graph.successors(rm)[0]
				edge1=graph.get_edge_data(start,rm)['signal']
				edge2=graph.get_edge_data(rm,finish)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(start,finish,signal='i')
				else:
					graph.add_edge(start,finish,signal='a')
			graph.remove_node(rm)
	#print(graph.nodes())
	flag=True
	return graph

#Upload KEGG codes modified for human pathways
def uploadKEGGcodes_hsa(codelist, graph, hsaDict, KEGGdict):
#queries the KEGG for the pathways with the given codes then uploads to graph. Need to provide the KEGGdict so that we can name the nodes with gene names rather than KO numbers
	for code in codelist:
		url=urllib2.urlopen('http://rest.kegg.jp/get/'+code+'/kgml')
		text=url.readlines()
		readKEGGhsa(text, graph, hsaDict, KEGGdict)
		#print(code)
		#print(graph.nodes())

def readKEGGhsa(lines, graph, hsaDict, KEGGdict):
	#read all lines into a bs4 object using libXML parser
	soup = BeautifulSoup(''.join(lines), 'lxml')
	groups = {} # store group IDs and list of sub-ids
	id_to_name = {} # map id numbers to names

	for entry in soup.find_all('entry'):
		entry_split= entry['name'].split(':')
		if len(entry_split)>2:
			if entry_split[0]=='hsa' or entry_split[0]=='ko':
				if entry_split[0]=='hsa':
					useDict=hsaDict
				else:
					useDict=KEGGdict
				nameList=[]
				entry_name=''
				namer=entry_split.pop(0)
				namer=entry_split.pop(0)
				namer=namer.split()[0]
				entry_name=entry_name+useDict[namer] if namer in useDict.keys() else entry_name+namer
				for i in range(len(entry_split)):
					nameList.append(entry_split[i].split()[0])
				for namer in nameList:
					entry_name=entry_name+'-'+useDict[namer] if namer in useDict.keys() else entry_name+'-'+namer
				entry_type = entry['type']
			else:
				entry_name=entry['name']
				entry_type = entry['type']
		else:
			if entry_split[0]=='hsa':
				entry_name=entry_split[1]
				entry_type = entry['type']
				entry_name = hsaDict[entry_name] if entry_name in hsaDict.keys() else entry_name
			elif entry_split[0]=='ko':
				entry_name=entry_split[1]
				entry_type = entry['type']
				entry_name = KEGGdict[entry_name] if entry_name in KEGGdict.keys() else entry_name
			elif entry_split[0]=='path':
				entry_name=entry['name']
				entry_type='path'
			else:
				entry_name=entry['name']
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

	for node in graph.nodes():
		if graph.degree(node)==0:
			graph.remove_node(node)

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

def parseKEGGdicthsa(filename, aliasDict, dict1):
	#reads KEGG dictionary of identifiers between orthology numbers and actual protein names and saves it to a python dictionary
	#extends above method to also keep track of gene names that are identical so we can recover those from input data as well
	#again, read in lines from the file, save ko number and gene name, identify them in the dictionary. 
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]=='D':
			kline=line.split('      ')
			kstring=kline[1]
			kline=kstring.split(' ')
			k=kline[0]
			nameline=kline[1]
			nameline=nameline.split(';')
			name=nameline[0]
			if ',' in name:
				nameline=name.split(',')
				name=nameline[0]
				for entry in range(1,len(nameline)):
					aliasDict[nameline[entry].strip()]=name
			dict1[k]=name
	return dict1

def parseKEGGdict(filename, aliasDict, dict1):
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
			dict1[k]=name
	return dict1

def readKEGG(lines, graph, KEGGdict):
	#read all lines into a bs4 object using libXML parser
	soup = BeautifulSoup(''.join(lines), 'lxml')
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
	#treating edge classes as nodes for now, will simplify later
	#constuct a BeautifulSoup4 parser object with the libxml backend
	#input in lines is given as an array, the join function converts to 1 long string
	#On windows this fails, try 'lxml' instead of 'lxml'
	soup = BeautifulSoup(''.join(lines), 'lxml')

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
				print("ambiguous direction relationship found")
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
	#On windows this fails, try 'lxml' instead of 'lxml'
	soup = BeautifulSoup(''.join(lines), 'lxml')
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
				if node_name==' ' or node_name=='':
					print("node with blank name")
					node_name=node_id
			except Exception as e:
				node_name = node_id

			#add the node to the graph
			graph.add_node(node_name, {'name': node_name, 'type': biopax_class, 'id': node_id})
			node_id_to_name[node_id] = node_name
			print(node_name)

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
					if from_node_id in node_id_to_name.keys():
						from_node_name = node_id_to_name[from_node_id]
					else:
						from_node_name= from_node_id
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
					if to_node_id in node_id_to_name.keys():
						to_node_name = node_id_to_name[to_node_id]
					else:
						to_node_name= to_node_id
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
					if from_node_id in node_id_to_name.keys():
						from_node_name = node_id_to_name[from_node_id]
					else:
						from_node_name= from_node_id
					edge_type = 'left'
					graph.add_edge(from_node_name, node_name, type=edge_type)
				for right in rights:
					resource = right.get('rdf:resource')
					to_node_id = resource[1:] if resource[0] == '#' else resource
					if to_node_id in node_id_to_name.keys():
						to_node_name = node_id_to_name[to_node_id]
					else:
						to_node_name= to_node_id
					edge_type = 'right'
					graph.add_edge(node_name, to_node_name, type=edge_type)

			#
			if biopax_class in ['Interaction', 'GeneticInteraction', 'MolecularInteraction', 'TemplateReaction']:
				print("ambiguous direction relationship found")
				participants = biopax_object.find_all('participant')
				for p1 in participants:
					resource = p1.get('rdf:resource')
					to_node_id = resource[1:] if resource[0] == '#' else resource
					if to_node_id in node_id_to_name.keys():
						to_node_name = node_id_to_name[to_node_id]
					else:
						to_node_name= to_node_id
					edge_type = 'participant'
					graph.add_edge(node_name, to_node_name, type=edge_type)
			if biopax_class in ['Pathway']:
				components = biopax_object.find_all('pathwayComponent')
				for component in components:
					resource = component.get('rdf:resource')
					to_node_id = resource[1:] if resource[0] == '#' else resource
					if to_node_id in node_id_to_name.keys():
						to_node_name = node_id_to_name[to_node_id]
					else:
						to_node_name= to_node_id
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
		#print(code)
		#print(graph.nodes())

# PC = Pathway Commons
def download_PC_codes(codelist, graph):
	for code in codelist:
		url = urllib2.urlopen('http://www.pathwaycommons.org/pc2/graph?source='+code+'&kind=neighborhood')
		text=url.readlines()
		read_biopax_dev(text, graph)
		#uncomment below and comment above to enable beta parser with names as indices
		# read_biopax_dev(text, graph)
		print(code)
	simplify_biopax_graph(graph)


def ifngStimTestSetup(params):

	aliasDict={}
	dict1={}
	parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)

	# read in list of codes then load them into network
	#inputfile = open('ID_filtered.c2.cp.kegg.v3.0.symbols.txt', 'r')
	inputfile = open('inputData/ID_filtered_IFNGpathways.txt', 'r')

	lines = inputfile.readlines()
	data=dict(utils.loadFpkms('inputData/Hela-C-1.count'))


	lines.pop(0)
	# lines.pop(0)
	# lines=[0]
	# lines=[lines[0]]
	# lines=['04110\n']
	codelist=[]
	for code in lines:
		graph=nx.DiGraph()
		coder=str('ko'+code[:-1])
		uploadKEGGcodes([coder], graph, dict2)
		coder=str('hsa'+code[:-1])
		uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		if(len(graph.edges())>1):
			graph=simplifyNetwork(graph, data)
		
		#graph = utils.simpleNetBuild()
		#coder='unrolled'

		#graph=utils.LiuNetwork1Builder()
		# coder='liu'

		# if the code has an interesting logic rule to find, run the algorithm... 
		checker=False
		for x in graph.in_degree():
			if x>1:
				checker=True
		if(checker):
			print(coder)
			codelist.append(coder)			
			nx.write_graphml(graph,coder+'.graphml')
			nx.write_gpickle(graph,coder+'.gpickle')
if __name__ == '__main__':
	params=sim.paramClass()
	ifngStimTestSetup(params)