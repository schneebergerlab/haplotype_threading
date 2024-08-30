#convert_nodekmer_to_fasta.py
import sys

inFile=sys.argv[1]
outFile=inFile[:-3]+".tsv"
out=open(outFile,"w+")
starting=True
started=False
oldline1=""
oldline2=""
for line in open(inFile,"r"):
	if starting and started: 
		out.write(oldline2+"\t"+oldline1+"\n")
		started=False
	#get string from thisline	
	val=line.rstrip()
	if starting:
		oldline1=val[1:] # remove fasta header bit
		starting=False
		started=True
	elif started and not starting:
		oldline2=val
		starting=True

#tidy up
if starting and started: 
	out.write(oldline2+"\t"+oldline1+"\n")
