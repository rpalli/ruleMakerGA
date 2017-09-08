# import python packages
import matplotlib.pyplot as plt
import numpy as numpy
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from sets import Set
import os as os
# import other parts of our code
import utils as utils

def analyzeExperiment(codes):
	importanceScores=[]
	# for each gpickle
	for code in codes:
		print(code)	
		importanceScore=pickle.Unpickler(open( 'scores/'+code+'_1'+'_scores.pickle', "rb" )).load()	
		print(importanceScore)
		importanceScores.extend(importanceScore)

		n, bins, patches = plt.hist(importanceScore,   facecolor='green')
		print(numpy.mean(importanceScore))

		print(numpy.std(importanceScore))
		print(numpy.median(importanceScore))
		plt.xlabel('Importance Score')
		plt.ylabel('Count')
		plt.title('Histogram of Importance Scores: '+code)
		# plt.axis([40, 160, 0, 0.03])
		plt.grid(True)

		plt.savefig(code+'_importance_hist.png', bbox_inches='tight')
		plt.clf()
	n, bins, patches = plt.hist(importanceScores,    facecolor='green')
	print(numpy.mean(importanceScores))
	print(numpy.std(importanceScores))
	print(numpy.median(importanceScores))
	plt.xlabel('Importance Score')
	plt.ylabel('Count')
	plt.title('Histogram of All Importance Scores')
	# plt.axis([40, 160, 0, 0.03])
	plt.grid(True)
	plt.savefig('total_importance_hist.png', bbox_inches='tight')




if __name__ == '__main__':
	codes=[]
	for file in os.listdir("gpickle"):
		if file.endswith(".gpickle"):
			codes.append(file[:-8])
	analyzeExperiment(codes)
