#sort_and_fill_nodekstats.py
import sys
node_k_stats= sys.argv[1] #See Analysis.sh
outFile=node_k_stats[:-4]+".sortedFilled.txt"
node_max=0
node_dict={}
for line in open(node_k_stats,"r"):
	vals=line.split("\t")
	if int(vals[0])>node_max:
		node_max=int(vals[0]) # assumes very last node isn't missing
	node_dict[int(vals[0])]=line

for i in range(0,node_max+1,1):
	if i not in node_dict:
		#print(i)
		node_dict[i]=str(i)+"\t"+"0\t0\t0.0\t0.0\t\n" #add a trailing tab here so that pandas can read it in as a dataframe even if first line is 0s
out=open(outFile,"w+")
for l in sorted(node_dict):
	out.write(node_dict[l])


