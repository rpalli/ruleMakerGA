def propBoolModel(individualParse, nodeList, inputOrderList, evaluateNodes, possibilityNumList, inputOrderInvertList,graph,sss1):
	async, iters, simSteps, generations, popSize, bitFlipProb, crossoverProb,mutationProb, mu,lambd, complexPenalty, genSteps=returnSimParams()
	bestList=[]
	initValueList=[]
	for j in range(0,len(sss1)):
		initValueList.append([])
	for i in range(0,len(nodeList)):
		for j in range(0,len(sss1)):
			ss=sss1[j]
			if  nodeList[i] in sss1[0].keys():
				initValueList[j].append(ss[nodeList[i]])
			else:
				initValueList[j].append(0.5)
	for i in range(0,len(nodeList)):
		print(nodeList[i])
		print(individualParse[i])
		currentDev=10*len(sss1)
		best=[]
		if possibilityNumList[i]>0:
			for j in range(0,int(2**int(individualParse[i][2]-individualParse[i][0]+1))):
				#if possibilityNumList[i]==1:
					# print('single possibility')
					# print(j)
				bits=[]
				bits=bitList(j)
				while len(bits)<(individualParse[i][2]-individualParse[i][0]+1):
					bits.insert(0,0)
				deviation=0
				for steadyStateNum in range(0,len(sss1)):
					derivedVal=propUpdateNode(i,initValueList[steadyStateNum],bits, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList, )
					deviation=deviation+(derivedVal-sss1[steadyStateNum][nodeList[i]])**2
				print(bits)
				print(deviation)
				if(deviation<currentDev):
					print("best")
					best=bits
					currentDev=deviation	
		# print(i)
		# print(nodeList[i])
		bestList.append(best)
		# print(best)
	return [item for sublist in bestList for item in sublist]

def propUpdateNet(currentNode,oldValue,individual, triple, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	inputOrder=inputOrderList[currentNode] # find the list of possible input combinations for the node we are on 
	inputOrderInvert=inputOrderInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if possibilityNumList[currentNode]>0:
		logicOperatorFlags=list(individual[triple[1]:triple[2]]) # find demarcations on individual for bits we need
		inputOrder=list(inputOrder[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]]) # determine what order of inputs is
		inputOrderInvert=inputOrderInvert[bit2int(individual[triple[0]:triple[1]])%possibilityNumList[currentNode]] # lookup which inputs need to be inverted ('not')
		if len(inputOrder)==0:
			value=oldValue[currentNode] #if no inputs, maintain value
			return value
		elif len(inputOrder)==1:
			#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
			if individual[triple[0]]==1:
				value=Inv(oldValue[inputOrder[0]],inputOrderInvert[0])
			else:
				value=oldValue[currentNode]
			return value
		else:
			#update nodes with more than one input
			# update and then or
			upstreamVals=[]
			for upstream in range(0,len(inputOrder)):
				upstreamVals.append(Inv(oldValue[inputOrder[upstream]],inputOrderInvert[upstream]))
			counter =0
			# print(upstreamVals)
			# print(logicOperatorFlags)
			while counter < len(logicOperatorFlags) and counter+1<len(inputOrder):
				if logicOperatorFlags[counter]==0:
					tempVal=fuzzyAnd(upstreamVals[counter],upstreamVals[counter+1])
					inputOrder.pop(counter)
					inputOrder.pop(counter)
					logicOperatorFlags.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.pop(counter)
					upstreamVals.insert(counter,tempVal)
				else:
					counter=counter+1
				# print(upstreamVals)

			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
			while len(upstreamVals)>1:
				tempVal=fuzzyOr(upstreamVals.pop(0),upstreamVals.pop(0))
				upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
			return upstreamVals[0]
	else:
		#returns savme value if now inputs
		return oldValue[currentNode]						

def propUpdateNode(currentNode,oldValue,individual, individualParse, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList):
	
	oldTriple=individualParse[currentNode]
	triple=[]
	triple.append(0)
	triple.append(oldTriple[1]-oldTriple[0])
	triple.append(oldTriple[2]-oldTriple[0])

	retinue = propUpdateNet(currentNode,oldValue,individual, triple, nodeList,inputOrderList, inputOrderInvertList, possibilityNumList)
	return retinue
					

def propOrE1(num1, num2, index1, index2,corrMat):
	return (num1+num2)*(2-corrMat[index1][index2])/2

def propAndE1(num1,num2,index1,index2,corrMat):
	return (num1+num2)*corrMat[index1][index2]/2
