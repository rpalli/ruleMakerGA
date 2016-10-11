import fuzzyNetworkConstructor as fnc
import networkx as nx

dict={}
aliasDict={}

#reads in the dictionary mapping of KEGG ID numbers to human-readable 
# gene names
KEGGdict=fnc.parseKEGGdict('ko00001.keg', aliasDict, dict)

#let us compare to some of the values from Rohith's original code:
assert len(KEGGdict) == 10319
assert KEGGdict['K12345'] == 'SRD5A3'
assert KEGGdict.keys()[0] == 'K01369'
assert KEGGdict.values()[0] == 'LGMN'
assert KEGGdict == dict # these two are identical for some reason
print "parseKEGGdict tests passed"
#OK, moving on

# test the reading of a KEGG file
graph = nx.DiGraph()
currentfile='testfiles/ko04060.xml'
inputfile = open(currentfile, 'r')
lines = inputfile.readlines()
fnc.readKEGG(lines, graph, KEGGdict)
assert len(graph.nodes()) == 212
assert graph.nodes()[123] == 'TNFSF10'
assert graph.nodes()[0] == 'LIF'
assert graph.edges()[-1] == ('VEGFC_D', 'KDR')
print "readKEGG", currentfile, "tests passed"

graph = nx.DiGraph()
fnc.readKEGGnew(lines, graph, KEGGdict)
assert len(graph.nodes()) == 212
assert graph.nodes()[123] == 'TNFSF10'
assert graph.nodes()[0] == 'LIF'
assert graph.edges()[-1] == ('VEGFC_D', 'KDR')
print "readKEGGnew", currentfile, "tests passed"

graph = nx.DiGraph()
currentfile='testfiles/ko04062.xml'
inputfile = open(currentfile, 'r')
lines = inputfile.readlines()
fnc.readKEGG(lines, graph, KEGGdict)
assert len(graph.nodes()) == 61
assert graph.nodes()[-1] == 'CDC42'
assert graph.edges()[0] == ('RASGRP2', 'RAP1A-RAP1B')
print "readKEGG", currentfile, "tests passed"

graph = nx.DiGraph()
currentfile = 'testfiles/glycolysis-biopax'
inputfile = open(currentfile, 'r')
lines = inputfile.readlines()
fnc.read_biopax(lines, graph)


