import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
networkNodeNums=[104,144,8,110,33,4]

finalNodeData=pickle.Unpickler(open( "node_by_node_data.pickle", "rb" )).load()
extendedNodeData=pickle.Unpickler(open( "extended_node_by_node_data.pickle", "rb" )).load()
[sensitivities,sensStd, specificities, specsStd]=pickle.Unpickler(open( "network_by_network_data.pickle", "rb" )).load()
# pickle.dump( truthholder, open( "expt_true_corrected_bits.pickle", "rb" ) )


sensitivities=[sens[0] for sens in sensitivities]
specificities=[spec[0] for spec in specificities]
senStd=[sens[0] for sens in sensStd]
specStd=[sens[0] for sens in specsStd]


sns.set_style("ticks")
sns.set_palette(sns.light_palette("purple/blue", input="xkcd", reverse=True))
threshold=.8

y_pos = np.arange(len(networkNodeNums))
plt.bar(y_pos, sensitivities, align='center', alpha=0.5,yerr=senStd, capsize=20)
plt.xticks(y_pos, networkNodeNums)
plt.ylabel('Sensitivity')
plt.xlabel('Node number')
plt.title('Sensitivity by Nodes in Pathway')
plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
sns.despine()
plt.show()
# Draw a nested barplot to show survival for class and sex
# g = sns.factorplot(x=nodeNums, y=nodesens, hue="sex", kind="bar", palette="muted")
# g.despine(left=True)
# g.set_ylabels("survival probability")
plt.bar(y_pos, specificities, align='center', alpha=0.5,yerr=specStd)
plt.xticks(y_pos, networkNodeNums)

plt.plot( [-.5,len(networkNodeNums)-.5],[threshold, threshold], "k--")
plt.ylabel('Specificity')
plt.xlabel('Node number')
plt.title('Specificity by Nodes in Pathway')
sns.despine()
plt.show()




nodeNums=[0, 1, 2, 3, 4, 5, 6, 7, 9]

nodesens=[sens[0] for sens in finalNodeData]
nodespec=[spec[3] for spec in finalNodeData]
nodeNums.pop(0)
nodesens.pop(0)
nodespec.pop(0)
nodeNums.pop(0)
nodesens.pop(0)
nodespec.pop(0)



# y_pos = np.arange(len(nodeNums))
# plt.bar(y_pos, nodesens, align='center', alpha=0.5)
# plt.xticks(y_pos, nodeNums)
# plt.ylabel('Sensitivity')
# plt.xlabel('In-degree')
# plt.title('Sensitivity by in-degree')
# plt.plot( [-.5,len(nodeNums)-.5],[threshold, threshold], "k--")
# sns.despine()
# plt.show()
# # Draw a nested barplot to show survival for class and sex
# # g = sns.factorplot(x=nodeNums, y=nodesens, hue="sex", kind="bar", palette="muted")
# # g.despine(left=True)
# # g.set_ylabels("survival probability")
# plt.bar(y_pos, nodespec, align='center', alpha=0.5)
# plt.xticks(y_pos, nodeNums)
# plt.plot( [-.5,len(nodeNums)-.5],[threshold, threshold], "k--")
# plt.ylabel('Specificity')
# plt.xlabel('In-degree')
# plt.title('Specificity by in-degree')
# sns.despine()
# plt.show()