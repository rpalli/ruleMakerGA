
import pickle
for fileName in ['hsa04380_unsimplified_1_output.pickle','hsa04612_unsimplified_1_output.pickle','hsa04630_unsimplified_1_output.pickle','hsa04650_unsimplified_1_output.pickle']:
	results=pickle.Unpickler(open( fileName, "rb" )).load()
	print(results[0])