import re
inputfile = open('noiseTest.txt', 'r')
linup=inputfile.read()
table=re.findall('10 \t[0-9]........................................',linup)
stuff=[]
for item in table:
	temp=re.split('\t',item)
	if len(temp)>0:
		stuff.append(temp[4])
print(stuff)