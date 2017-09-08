import numpy as np
import random as rand
import scipy.stats as stat
import sklearn.mixture as mix

def readCounts(geneDict, filename, samplesAdded):
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		newline=line.split('\t')
		if newline[0] not in geneDict.keys():
			geneDict[newline[0]]=[]
			for i in range(samplesAdded):
				geneDict[newline[0]].append(-1)
		geneDict[newline[0]].append(int(newline[1]))
	samplesAdded+=1
	return samplesAdded

def readData():
	tempGeneDict={}
	samples=0
	for j in range(1,5):
		samples=readCounts(tempGeneDict, 'murphyData/JAR-C-'+str(j)+'.count', samples)
	for j in [1,2,4]:
		samples=readCounts(tempGeneDict, 'murphyData/JAR-d-'+str(j)+'.count', samples)
	for j in range(1,5):
		samples=readCounts(tempGeneDict, 'murphyData/JAR-d_PV-'+str(j)+'.count', samples)
	for j in range(1,5):
		samples=readCounts(tempGeneDict, 'murphyData/JAR-PV-'+str(j)+'.count', samples)
	for j in range(1,5):
		samples=readCounts(tempGeneDict, 'murphyData/JAR-V-'+str(j)+'.count', samples)
	for j in range(1,5):
		samples=readCounts(tempGeneDict, 'murphyData/JAR-d_V-'+str(j)+'.count', samples)
	
	# subset those keys with high values and high coefficient of variation
	keyList=[]


	for key in tempGeneDict.keys():
		mixer= mix.GaussianMixture(n_components=2, covariance_type='full').fit([[a] for a in tempGeneDict[key]])
		difference=abs(mixer.means_[0]-mixer.means_[1])
		if difference>50 and np.mean(tempGeneDict[key])>100:
			keyList.append(key)
	print(len(keyList))
	geneDict= dict((k, tempGeneDict[k]) for k in keyList)

	# return tempGeneDict
	return geneDict

def constructOmicsInput(geneDict, maxRem, divider):
	ssDict={}
	ones=[]
	zeros=[]
	valueLister=[]
	for key in geneDict.keys():
		values=list(geneDict[key])	
		if maxRem:
			values.remove(max(values))
		CAND2= [[entry] for entry in values]
		if divider:
			maxVal=max(values)
			if maxVal>0:
				props=[max(0,min(1,val/maxVal)) for val in geneDict[key]]
			else:
				props=[0 for val in geneDict[key]]
				print('zeros only found:')
				print(key)
		else:
			mixer= mix.GaussianMixture(n_components=2, covariance_type='full').fit(CAND2)
			if (mixer.means_[0]>mixer.means_[1]):
				mew1=mixer.means_[0][0]
				mew2=mixer.means_[1][0]
			else:
				mew2=mixer.means_[0][0]
				mew1=mixer.means_[1][0]
			props=[max(0,min(1,(val-mew2)/(mew1-mew2))) for val in geneDict[key]]
		ssDict[key]=list(props)
		valueLister.extend(props)
		ones.append(props.count(1))
		zeros.append(props.count(0))
	return ssDict, valueLister, ones, zeros
