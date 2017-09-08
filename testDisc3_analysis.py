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
	dfsetuplist=[]

	# for each gpickle
	for code in codes:
		print(code)		

		#FOR EACH RUN WITH THAT GPICKLE
		for i in range(1,6):
			edgeDegree=[]
			nodesensitivities=[[],[],[],[]]
			nodespecificities=[[],[],[],[]]


			[[dev1,dev3],[bruteOut1,bruteOut3],[model1,model3]]=pickle.Unpickler(open( 'pickles/'+code+'_'+str(i)+'_output.pickle', "rb" )).load()
			dfsetuplist.append({'type':'eliminate max','value':numpy.mean(dev3),'net':code,'nodenum':str(len(model3.nodeList))})
			dfsetuplist.append({'type':'keep max','value':numpy.mean(dev1),'net':code,'nodenum':str(len(model1.nodeList))})

	df=pd.DataFrame(dfsetuplist)


	# sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})
	sns.set_style("dark")

	# plt.xlabel('In-degree', color='gray')
	# plt.title('Sensitivity by in-degree')

	sns.set(font_scale=2) 

	sns.set_style("dark")
	ax=sns.swarmplot(data=df,x='net', y='value', hue='type')
	# ax.set_xticklabels(nodeNums, color='white')
	# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')
	sns.set_style("dark")
	plt.xticks( rotation='vertical')
	plt.ylabel('Node-wise Error')
	plt.xlabel('Network')
	plt.title('Murphy Data Fits by Network', color='black')
	plt.ylim([0,1])
	plt.tight_layout()
	# sns.despine()

	plt.savefig('division_prob_fits.png', bbox_inches='tight', transparent=True)


	plt.clf()




	sns.set_style("dark")


	sns.set(font_scale=2) 

	sns.set_style("dark")
	ax=sns.swarmplot(data=df,x='nodenum', y='value', hue='type')
	sns.set_style("dark")
	plt.ylabel('Node-wise Error')
	plt.xlabel('Network Size (Node Number)')
	plt.title('Murphy Data Fits by Network', color='black')
	plt.ylim([0,1])
	plt.tight_layout()
	# sns.despine()

	plt.savefig('node_num_fits.png', bbox_inches='tight', transparent=True)


	plt.clf()

if __name__ == '__main__':
	codes=[]
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			codes.append(file[:-8])
	print(codes)
	analyzeExperiment(codes)
