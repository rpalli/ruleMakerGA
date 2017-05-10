import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
networkNodeNums=[18,32,7,23,68,78,72,28]

finalNodeData=pickle.Unpickler(open( "node_by_node_data.pickle", "rb" )).load()
extendedNodeData=pickle.Unpickler(open( "extended_node_by_node_data.pickle", "rb" )).load()
[sensitivities,sensStd, specificities, specsStd]=pickle.Unpickler(open( "network_by_network_data.pickle", "rb" )).load()



sensitivities=[sens[3] for sens in sensitivities]
specificities=[spec[3] for spec in specificities]
senStd=[sens[3] for sens in sensStd]
specStd=[sens[3] for sens in specsStd]
specStd.pop(0)
senStd.pop(0)
sensitivities.pop(0)
specificities.pop(0)
print(sensitivities)
sns.set_style("ticks")
sns.set(font_scale=2) 
sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
threshold=.8

y_pos = np.arange(len(networkNodeNums))
plt.bar(y_pos, sensitivities, align='center', alpha=0.5,yerr=senStd, capsize=20)
plt.xticks(y_pos, networkNodeNums)
plt.ylabel('Sensitivity')
plt.xlabel('Node number')
plt.title('Sensitivity by Nodes in Pathway')
plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
plt.tight_layout()
sns.despine()
plt.savefig('IFNG_network_sensitivities.png', bbox_inches='tight')

plt.clf()

plt.bar(y_pos, specificities, align='center', alpha=0.5,yerr=specStd)
plt.xticks(y_pos, networkNodeNums)

plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
plt.ylabel('Specificity')
plt.xlabel('Node number')
plt.title('Specificity by Nodes in Pathway')
sns.despine()
plt.tight_layout()
plt.savefig('IFNG_network_specificities.png', bbox_inches='tight')







plt.clf()

# nodeNums=[2, 3, 4, 5, 6, 7, 9]
nodeNums=[0,1,2,3,4,5,6]

nodeArray=[]
extendedNodeData.pop(0)
for tempExtended in extendedNodeData:
	nodeArray.append(tempExtended[3])
for array in nodeArray:
	if len(array)==0:
		array.append(0)
# nodeArray.pop(0)
nodeArray.pop(0)
nodeNums.pop(0)
# nodeArray.pop(0)

dfdictlist=[]
for i in range(len(nodeArray)):
	for element in nodeArray[i]:
		dfdictlist.append({'In-degree':nodeNums[i],'Sensitivity':element})
df=pd.DataFrame(dfdictlist)

# plt.title('Sensitivity by in-degree')

#print(df)
sns.set(font_scale=2) 
sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})

ax=sns.swarmplot(data=df,x='In-degree', y='Sensitivity', color='orange')

plt.ylabel('Sensitivity')
plt.xlabel('In-degree')
plt.title('Sensitivity by in-degree', color='white')
plt.ylim([0,1])
plt.tight_layout()
plt.ylim([0,1])
plt.savefig('IFNG_clear_sens.png', bbox_inches='tight', transparent=True)
plt.clf()


nodeNums=[0,1,2,3,4,5,6]

nodeArray=[]
for tempExtended in extendedNodeData:
	nodeArray.append(tempExtended[7])
for array in nodeArray:
	if len(array)==0:
		array.append(0)
nodeArray.pop(0)
nodeNums.pop(0)
nodeArray.pop(0)
nodeNums.pop(0)
# nodeArray.pop(0)


dfdictlist=[]
for i in range(len(nodeArray)):
	for element in nodeArray[i]:
		dfdictlist.append({'In-degree':nodeNums[i],'Specificity':element})
df=pd.DataFrame(dfdictlist)


sns.set_style("dark",{'axes.edgecolor': 'white', 'axes.facecolor': 'white',  'axes.grid': True,'axes.labelcolor': 'white',  'figure.facecolor': 'white',  'grid.color': '.2',  'text.color': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})


plt.xlabel('In-degree', color='white')
# plt.title('Sensitivity by in-degree')

sns.set(font_scale=2) 


ax=sns.swarmplot(data=df,x='In-degree', y='Specificity', color='orange')
# ax.set_xticklabels(nodeNums, color='white')
# ax.set_yticklabels(np.arange(0,1.2,.2),color='white')

plt.ylabel('Specificity')
plt.xlabel('In-degree')
plt.title('Specificity by in-degree', color='white')
plt.ylim([0,1])
plt.tight_layout()
# sns.despine()

plt.savefig('IFNG_clear_spec.png', bbox_inches='tight', transparent=True)


plt.clf()