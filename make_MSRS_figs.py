import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
networkNodeNums=[104,144,8,110,33,4]

finalNodeData=pickle.Unpickler(open( "node_by_node_data.pickle", "rb" )).load()
extendedNodeData=pickle.Unpickler(open( "extended_node_by_node_data.pickle", "rb" )).load()
[sensitivities,sensStd, specificities, specsStd]=pickle.Unpickler(open( "network_by_network_data.pickle", "rb" )).load()
#pickle.dump( truthholder, open( "expt_true_corrected_bits.pickle", "rb" ) )



# sensitivities=[sens[0] for sens in sensitivities]
# specificities=[spec[0] for spec in specificities]
# senStd=[sens[0] for sens in sensStd]
# specStd=[sens[0] for sens in specsStd]




# sns.set_style("ticks")
# sns.set(font_scale=2) 
# sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
# threshold=.8

# y_pos = np.arange(len(networkNodeNums))
# plt.bar(y_pos, sensitivities, align='center', alpha=0.5,yerr=senStd, capsize=20)
# plt.xticks(y_pos, networkNodeNums)
# plt.ylabel('Sensitivity')
# plt.xlabel('Node number')
# plt.title('Sensitivity by Nodes in Pathway')
# plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
# plt.tight_layout()
# sns.despine()s
# plt.show()



# plt.bar(y_pos, specificities, align='center', alpha=0.5,yerr=specStd)
# plt.xticks(y_pos, networkNodeNums)

# plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
# plt.ylabel('Specificity')
# plt.xlabel('Node number')
# plt.title('Specificity by Nodes in Pathway')
# sns.despine()
# plt.tight_layout()
# plt.show()


# nodeNums=[2, 3, 4, 5, 6, 7, 9]
nodeNums=[0,1,2]

nodeArray=[]
for tempExtended in extendedNodeData:
	nodeArray.append(tempExtended[0])
for array in nodeArray:
	if len(array)==0:
		array.append(0)
# nodeArray.pop(0)
# nodeArray.pop(0)
for element in nodeArray:
	print(len(element))


sns.set_style("ticks")
sns.set(font_scale=2) 
threshold=.8
sns.set_style("whitegrid")
sns.set(font_scale=2) 
# sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
ax=sns.boxplot(nodeArray,nodeNums, color=sns.light_palette("purple/blue", input="xkcd", reverse=True))
plt.tight_layout()
plt.ylabel('Sensitivity')
plt.xlabel('In-degree')
plt.title('Sensitivity by in-degree')
ax.set_xticklabels(nodeNums)
plt.tight_layout()
plt.show()



nodeArray=[]
for tempExtended in extendedNodeData:
	nodeArray.append(tempExtended[3])
for array in nodeArray:
	if len(array)==0:
		array.append(0)
# nodeArray.pop(0)
# nodeArray.pop(0)
plt.ylabel('Specificity')
plt.xlabel('In-degree')
plt.title('Specificity by in-degree')

ax=sns.boxplot(nodeArray,nodeNums, color=sns.light_palette("purple/blue", input="xkcd", reverse=True))
ax.set_xticklabels(nodeNums)

plt.tight_layout()
plt.show()
