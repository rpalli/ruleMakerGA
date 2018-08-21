
from utils import genInitValueList, setupEmptyKOKI
from GA import GAsearchModel, localSearch

if __name__ == '__main__':
	import time
	start_time = time.time()
	
	# read in arguments from shell scripts
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("iterNum")
	results = parser.parse_args()
	graphName=results.graph
	iterNum=int(results.iterNum)
	name=graphName[:-8]+'_'+results.iterNum
	graph = nx.read_gpickle(graphName)
	
	# load data
	sampleList=pickle.Unpickler(open( graphName[:-8]+'_sss.pickle', "rb" )).load()
	
	# set up parameters of run, model
	params=paramClass()
	model=modelClass(graph,sss, False)
	samples=params.samples

	storeModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	
	# put lack of KOs, initial values into correct format
	knockoutLists, knockinLists= setupEmptyKOKI(samples)
	newInitValueList=genInitValueList(sampleList,model)
	model.initValueList=newInitValueList

	# find rules
	model1, dev, bruteOut =GAsearchModel(model, sss, params, knockoutLists, knockinLists, name) # run GA
	bruteOut1, equivalents = localSearch(model1, bruteOut, sss, params, knockoutLists, knockinLists) # run local search
	storeModel3=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	outputList=[bruteOut1,dev,storeModel, storeModel3, equivalents]
	pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) ) # output rules


	# calculate importance scores and output
	scores1=calcImportance(bruteOut1,params,model1, sss)
	pickle.dump( scores, open( name+"_scores1.pickle", "wb" ) )

	# write rules
	with open(name+"_rules.txt", "w") as text_file:
		text_file.write(utils.writeModel(bruteOut1, model1))
	print("--- %s seconds ---" % (time.time() - start_time))