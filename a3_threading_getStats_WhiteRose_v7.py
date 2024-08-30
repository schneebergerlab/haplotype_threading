#getStatsTable.py
import sys

chrs=sys.argv[1]
chr_size=sys.argv[2]
inDir=sys.argv[3]
outFile=open(sys.argv[4]+"/WhRStats_V7.tsv","w+")
print("Pericentromeres","Chromosome","Total WhR Nodes","Visible WhR Nodes","Visible WhR(%)","Thread(2+)Nodes Accurate","Thread(2+)Nodes Total","Thread(2+) Accuracy(%)","Thread(2+)Nodes_noInversions Accurate","Thread(2+)Nodes_noInversions Total","Thread(2+) Accuracy(%)_noInversions","HapSwitches(WhR-to-WhR,2nodes+)","HapSwitches_nodesAffected","N50","thread(2+) coverage","windowCov(0)","windowCov(1)","windowCov(2)","windowCov(3)","windowCov(4)",sep="\t",file=outFile)
chrMbs=[]
for m in chr_size.rstrip().split(","):
	chrMbs.append(float(m))
for cdx,c in enumerate(chrs.rstrip().split(",")):
	chr_this_Mb=chrMbs[cdx]
	visibleStat=0
	Accuracy2plusStat=0
	N50=0
	maxCoverage2plusSeen=0.0
	HapswitchCount=0
	HapswitchNodes=0
	acc_file=open(inDir+"/"+c+"_div0.1.WhRaccuracy.tsv","r")
	for ldx,line in enumerate(acc_file):
		if ldx==0:
			totalWhRStat=int(line.rstrip().split("\t")[1])
			visibleWhRStat=int(line.rstrip().split("\t")[2])
			visibleWhR_pc=float(line.rstrip().split("\t")[3])
		if ldx==1:
			Accurate2plusStat=int(line.rstrip().split("\t")[1])
			Total2plusStat=int(line.rstrip().split("\t")[2])
			Accuracy2plus_pc=float(line.rstrip().split("\t")[3])
		if ldx==2:
			Accurate2plusStat_nonInverted=int(line.rstrip().split("\t")[1])
			Total2plusStat_nonInverted=int(line.rstrip().split("\t")[2])
			Accuracy2plus_nonInverted_pc=float(line.rstrip().split("\t")[3])
		if ldx==4:
			HapswitchCount=int(line.rstrip().split("\t")[1])
			HapswitchNodes=int(line.rstrip().split("\t")[3])
	thread_file=open(inDir+"/"+c+"_div0.1.threads.tsv","r")
	foundN50=False
	for ldx,line in enumerate(thread_file):
		if foundN50==False:
			if float(line.rstrip().split("\t")[8]) >=0.5:
				foundN50=True
				N50=float((line.rstrip().split("\t")[1]))/10
		if int(line.rstrip().split("\t")[1]) >=2:
			maxCoverage2plusSeen=float(line.rstrip().split("\t")[7])
			maxCoverage2plus_pc=float(line.rstrip().split("\t")[7])/chr_this_Mb


	windowCov_file=open(inDir+"/"+c+"_div0.1.windowCov.tsv","r")
	windowCov0=0
	windowCov1=0
	windowCov2=0
	windowCov3=0
	windowCov4=0
	for ldx,line in enumerate(windowCov_file):
		cov,num=line.rstrip().split("\t")
		#print(c,cov,num)
		if int(cov)==0:
			windowCov0=int(num)
		elif int(cov)==1:
			windowCov1=int(num)
		if int(cov)==2:
			windowCov2=int(num)
		elif int(cov)==3:
			windowCov3=int(num)
		if int(cov)==4:
			windowCov4=int(num)

	print("excluded",c,totalWhRStat,visibleWhRStat,round(visibleWhR_pc,3),Accurate2plusStat,Total2plusStat,round(Accuracy2plus_pc,3),Accurate2plusStat_nonInverted,Total2plusStat_nonInverted,round(Accuracy2plus_nonInverted_pc,3), HapswitchCount,HapswitchNodes,round(N50,1),round(maxCoverage2plusSeen,1),round(chr_this_Mb,1),round(maxCoverage2plus_pc,2),windowCov0,windowCov1,windowCov2,windowCov3,windowCov4,sep="\t",file=outFile)
