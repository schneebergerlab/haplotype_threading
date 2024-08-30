#Script for running custom EM approach to calling copy number of nodes in a haplotype graph using kmers

from collections import Counter
import copy
import numpy as np
import argparse
import itertools
import scipy as scipy
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker

version=0.4

class Window:
	def __init__(self, start, end):
		self.start=int(start)
		self.end=int(end)
		self.Nodes=[] # make this a set later

	def get_pos(self):
		return (self.start,self.end)
	def addNode(self, N):
		self.Nodes.append(N)
	def getNodes(self):
		return self.Nodes

class Node:
	def __init__(self, ID, totalMarkers, coveredMarkers, KmerCounts, genomes):
		self.ID=ID
		self.totalMarkers=totalMarkers
		self.coveredMarkers=coveredMarkers
		self.KmerCounts=KmerCounts # A counter dictionary
		self.genomes=set(genomes)
		self.EMPriors=None # Will hold the priors of each step of the EM algorithm, and can also be accessed to compare likelihoods of two competing blocks.
		self.isConverged=False # tells us whether the EM is done here
		self.iterationsDone=None
		self.Contributions=None #Will hold the first-pass contributions of all kmers within this block, to this block, which can then be accessed int he second pass to update EM priors.
		self.neighbour_Contributions=None #will hold sum of contributions from neighbouring nodes in a given round of EM. kept separate here for a while so that neighbour updates don't influence this nodes contributions to its own neighbours after this value is set
		self.relative_amps={}
	def get_ID(self):
		return self.ID

	def get_KmerCounts(self):
		return self.KmerCounts
	
	def get_EMpriors(self):
		return copy.copy(self.EMPriors)	

	def get_zeroCount(self):
		return self.totalMarkers-self.coveredMarkers

	def get_covFraction(self):
		if self.coveredMarkers !=0:
			return self.coveredMarkers/self.totalMarkers
		else:
			return 0

	def get_totalMarkers(self):
		return self.totalMarkers
	#takes a list of relative amplitudes and stores a dictionary. Assumes list length is equal to number of ploidies we are considering
	def set_relativeAmps(self,ra, Ploidy_order):
		self.relative_amps={}
		for p in Ploidy_order:
			self.relative_amps[p]=ra[p]

	def get_relativeAmps(self):
		return copy.copy(self.relative_amps)	

	def getGenomes(self):
		return copy.copy(self.genomes)

	def check_isConverged(self):
		return self.isConverged

	def set_EMpriors(self, p):
		self.EMPriors=copy.deepcopy(p)

	def update_EMpriors_withChecks(self,newPriors,minUpdate,iterationCount, Ploidy_order):
		isDone=True
		oldPriors=self.EMPriors
		for p in Ploidy_order:
			if abs(oldPriors[p]-newPriors[p]) >= minUpdate: # If our update step is bigger than the minimum update (still haven't converged), then cancel the signal that we're done
				isDone=False
			self.EMPriors[p] = newPriors[p]
		if isDone:
			self.isConverged=True
			self.iterationsDone=iterationCount

	def set_Contributions(self,c):
		self.Contributions = copy.deepcopy(c)

	def get_Contributions(self):
		return copy.copy(self.Contributions)

	def set_NeighbourContributions(self,c):
		self.neighbour_Contributions = copy.deepcopy(c)

	def get_NeighbourContributions(self):
		return copy.copy(self.neighbour_Contributions)

	def calc_neighbour_contributions(self,Ploidy_order):
		#if self.ID=="2232":
		neighbour_l_vals=[]
		for w in self.neighbour_Contributions: # for each window distance from this node
			w_sum=0
			these_conts=self.neighbour_Contributions[w]

			for c in these_conts:
				w_sum+=sum(c.values()) # build up sum of all contributions for this window.
			#now we want to calculate the likelihood of copy number at this window given the ?product? of likelihoods from nodes of this window.
			likelihood={}
			for p in Ploidy_order:
				likelihood[p]=0

			combinations=list(itertools.product(*these_conts))
			for c in combinations:
				if sum(c) in Ploidy_order: #if this is a realistic combination of ploidy-counts
					vals=[]
					for xdx, x in enumerate(c): #eg we are interating through (1,3,0) giving, 0,1 then 1,3 then 2,0
						vals.append(these_conts[xdx][x]) # eg these_conts[0] {0:10, 1:5,2:0,3:0,4:0} we will pull out 10
					likelihood[sum(c)]+=np.prod(vals) # multiply the contributions of these different nodes together




			#normalise by the maximum contributions of this window
			likelihood_sum=sum(likelihood.values())
			for l in likelihood:
				likelihood[l]= (likelihood[l]/likelihood_sum)*w_sum
			#print(w,w_sum, likelihood)#,these_conts)
			neighbour_l_vals.append(likelihood)

		#summarise and return sum of contributions for each likelihood
		output={}
		for p in Ploidy_order:
			output[p]=0
		for ldict in neighbour_l_vals:
			for p in ldict:
				output[p]+=ldict[p]


		return output

def div_by_zero(a,b):
	if b != 0 and b != 0.0:
		return a/b
	else:
		return 0.0
def NormaliseDictionaryValues(d):

	norm_factor=div_by_zero(1.0,sum(d.values()))
	for k in d:
		d[k] = d[k]*norm_factor
	return d
	#Thanks StackOverflow user Benedict & anilbey

#dist is the distance of the two blocks we are analysing
#maxNeighbourDist is the point at which our decay line should hit 1% of yint
#yint is the value at which the line should intersect x=1 (basically 0 distance, or maximum possible weight)
def calc_weight_function_ExpDecay(dist, maxNeighbourDist, yint):
	if yint==0.0:
		return 0.0 # returns nan otherwise??
	#if dist<1: 
	#	dist=1
	#Assumes y=a*b^x
	#we'll model the first point as (1,yint) instead of (0,yint) for convenience of my math
	#we'll model the other point as (maxNeighbourDist,0.01)
	#solve for a and b
	yint=float(yint)
	b=((0.01*yint)/yint)**(1/maxNeighbourDist)
	a=yint/b
	#calculate weight at this distance
	weight=a*b**abs(dist)
	
	if weight >0:
		return weight
	else:
		return 0.0

def _4gaussian_fixedCen_plusError_fixedSigma(x_array, k2,b2, amp1,cen1,sigma1,amp2,amp3,amp4):
	cen2=cen1*2
	cen3=cen1*3
	cen4=cen1*4
	sigma2=(sigma1*2)/np.sqrt(2)
	sigma3=(sigma1*3)/np.sqrt(3)
	sigma4=(sigma1*4)/np.sqrt(4)
	return [k2*(b2**(x-0)) + \
	amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
	amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-(cen1*2))/sigma2)**2))) + \
	amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-(cen1*3))/sigma3)**2))) + \
	amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-(cen1*4))/sigma4)**2))) for x in x_array]
     #0 replaces c2

def _1gaussian(x_array, amp1,cen1,sigma1):
    return [amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) for x in x_array]

def _1exponential(x_array,k2,b2,c2):
    return [k2*(b2**(x-c2)) for x in x_array]


def get_probs_fromModel(count, model_params):
	tmp_probs={}
	for p in model_params: #model_params is a dictionary
		if p ==0:
			k= model_params[p][0]
			b= model_params[p][1]
			c= model_params[p][2]
			x=count
			if k * (b**(x-c)) < 1.0:
				tmp_probs[p]=1.0 # there is a baseline here such that any curve prdeicting less than 1 kmer count in the whole dataset (something like 1 in a million)
			else:
				tmp_probs[p]= k * (b**(x-c))
		else: #calcualte a gaussian model for it
			amp1= model_params[p][0]
			cen1= model_params[p][1]
			sigma1= model_params[p][2]
			x=count
			tmp_probs[p]= amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
	#Normalise Dictionary values
	norm_factor=1.0/sum(tmp_probs.values())
	for t in tmp_probs:
		tmp_probs[t] = tmp_probs[t]*norm_factor
	return tmp_probs

#######SOF############################################################################################################################################
def call(graphFile,nodeStatsFile,modelsFile,outFile,plotFolder,is_nodePlots,is_calcPriors,minUpdate,iterMax,neighbour_distance_max,neighbour_weight_max,maxK):
	print("neighbour_distance_max:",neighbour_distance_max)
	weights_example={}
	ex_s=""
	for r in range(0,10,1):
		weights_example[r*100000]=calc_weight_function_ExpDecay(r*100000, neighbour_distance_max, neighbour_weight_max)
	#print(weights_example)
	for w in sorted(weights_example):
		ex_s+="{0}:{1}\t".format(w, round(weights_example[w],5))
	print("weights ({0}):".format(neighbour_weight_max), ex_s,sep="\t")
	Ploidy_order=[0,1,2,3,4]	
	
	#Read model parameters into memory
	ploidy_levels=set()
	#Read in parameters from modelsFile
	Pars = {} #dictionary of parameters for models
	for ldx, l in enumerate(open(modelsFile,"r")):
		if ldx>0: # skip header
			vals=l.rstrip().split("\t")
			ploidy=int(vals[0])
			ploidy_levels.add(ploidy)
			Pars[ploidy]=[float(vals[1]),float(vals[2]),float(vals[3])]
	print("probability model \n",Pars)
	
	
	#Keep a list of windows, each window has a set of nodes
	
	Windows=[]
	Nodes=[]
	currentWin=-1
	windx=-1
	num_nodes=0

	for g, n in zip(open(graphFile),open(nodeStatsFile)):
		g=g.rstrip()
		n=n.rstrip()
		#print("{0}\n{1}".format(g,n[:10]))
		graphVals=g.split("\t")
		nodeVals=n.split("\t")
		win_start,win_end=graphVals[0:2]
		if win_start != currentWin: # if new window make new window
			Windows.append(Window(win_start,win_end))
			windx +=1 # in first iter, windx will now point to 0
			currentWin=win_start
		genomes=graphVals[2].split(",")
		kmerCounts=Counter()
		if len(nodeVals)==6: #if this is a full kmer containing line
			for k in nodeVals[5].split(","):
				if k != "" and int(k) <= maxK:
					kmerCounts[int(k)]+=1
		kmerCounts[0]+=int(nodeVals[1])-int(nodeVals[2]) # add zeroes for everything not recorded in kmer list
		newnode=Node(ID=nodeVals[0],totalMarkers=int(nodeVals[1]),coveredMarkers=int(nodeVals[2]),KmerCounts=kmerCounts,genomes=genomes)

		Windows[windx].addNode(newnode)

		num_nodes+= 1



	#Find good priors for each node!##############################################################
	if is_calcPriors:
		print("Finding priors for each node")
		initial_prob=1.0/len(Ploidy_order)
		#for each node, fit amplitudes of the model 
		x_array=list(range(0,maxK,1))
	
		windowSize=Windows[0].get_pos()[1]-Windows[0].get_pos()[0]
		out=open(outFile,"w+")
		for windx,currentWin in enumerate(Windows):
			if windx%50==0:
				print(windx,"/",len(Windows))
			for ndx, node in enumerate(Windows[windx].getNodes()):
				node_id=node.get_ID()
				if node_id==node_id:#"1244":
					Kc=node.get_KmerCounts()
					#print(Kc)
					y_array=[]
					for x in x_array:
						y_array.append(Kc[x])
					#print(node_id,y_array)
		
					def _helper(xarr,k2,b2,amp1,amp2,amp3,amp4): #calls the guassian fit function with parameters fixed
						return _4gaussian_fixedCen_plusError_fixedSigma(xarr,k2,b2,amp1,Pars[1][1],Pars[1][2],amp2,amp3,amp4)
	
					a1_guess=max(y_array[:4])
					a2_guess=max(y_array[:4])/3
					a3_guess=max(y_array[:4])/15
					a4_guess=max(y_array[:4])/1000
					popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_helper, x_array, y_array, p0=[y_array[0],Pars[0][1],a1_guess,a1_guess,a1_guess,a1_guess],bounds=(0,10000000)) #ignore zero-three for fitting the first peak
					#print(popt_gauss)
					pa_0 = popt_gauss[0:2]
					a1 = popt_gauss[2]
					a2 = popt_gauss[3]
					a3 = popt_gauss[4]
					a4 = popt_gauss[5]
					amps=[pa_0[0],a1,a2,a3,a4]
					relative_amps=[round(a/sum(amps),3) for a in amps]
					#print(amps)
					#print(relative_amps)
					node.set_relativeAmps(relative_amps, Ploidy_order)
					#print(node.get_relativeAmps())
	
	
					#Plot tmp1
					if is_nodePlots:
						ras=node.get_relativeAmps()
	
						#calculate expected values
						gauss_peak_4_fixed = _4gaussian_fixedCen_plusError_fixedSigma(x_array, pa_0[0],pa_0[1],a1,Pars[1][1],Pars[1][2],a2,a3,a4)
						#calculate goodness of fit
						fig = plt.figure(figsize=(6,4))
						gs = gridspec.GridSpec(1,1)
						ax1 = fig.add_subplot(gs[0])
						ax1.set_xlim([0,maxK])
						yplotmax=max(y_array[1:])
						if yplotmax < 1:
							yplotmax=y_array[0]*1.2
						if yplotmax <1:
							yplotmax=1/1.2
						ax1.set_ylim([0,yplotmax*1.2])
						ax1.text(3,yplotmax*1.1,"y_intercept ={}".format(round(pa_0[0],2)))
						ax1.text(maxK/2-10,yplotmax*1.22,"Node: {}".format(node.get_ID()),size=15)
		
						ax1.plot(x_array, y_array, "ro")
						ax1.plot(x_array, gauss_peak_4_fixed,lw=1, color="k")
						#plt.show()
		
		
						noise_decay= _1exponential(x_array,*pa_0,0.0)
						gauss_peak_hap= _1gaussian(x_array, a1,Pars[1][1],Pars[1][2])
						gauss_peak_dip= _1gaussian(x_array, a2,Pars[2][1],Pars[2][2])
						gauss_peak_trip=_1gaussian(x_array, a3,Pars[3][1],Pars[3][2])
						gauss_peak_tet= _1gaussian(x_array, a4,Pars[4][1],Pars[4][2])
		
						ax1.text(maxK-20,yplotmax*0.5,"0-copy: {}".format(round(ras[0],2)))
						ax1.text(maxK-20,yplotmax*0.6,"1-copy: {}".format(round(ras[1],2)))
						ax1.text(maxK-20,yplotmax*0.7,"2-copy: {}".format(round(ras[2],2)))
						ax1.text(maxK-20,yplotmax*0.8,"3-copy: {}".format(round(ras[3],2)))
						ax1.text(maxK-20,yplotmax*0.9,"4-copy: {}".format(round(ras[4],2)))
						ax1.text((maxK/2)-10,yplotmax*1.1,"Kmers = {}".format(sum(y_array)))
		
		
						#overcount= [1.0]*len(xplot)
						colours=["grey","blue","orange","green","purple"]
						ax1.plot(x_array, noise_decay,color=colours[0],ls="--",lw=2)
						ax1.plot(x_array, gauss_peak_hap,color=colours[1],ls="--",lw=2)
						ax1.plot(x_array, gauss_peak_dip,color=colours[2],ls="--",lw=2)
						ax1.plot(x_array, gauss_peak_trip,color=colours[3],ls="--",lw=2)
						ax1.plot(x_array, gauss_peak_tet,color=colours[4],ls="--",lw=2)
						fig.tight_layout()
						#print(outFile+"."+str(node.get_ID())+".pdf")
						plt.savefig(plotFolder+"."+str(node.get_ID())+".pdf", transparent=False, dpi=300, format="pdf")
						plt.close()
						#plt.show()
	else:
		print("Using equal probabilities as priors")

	###Do the thing!################################################################
	initial_prob=1.0/len(Ploidy_order) # at the start, all ploidy levels for each block are equally likely
	
	windowSize=Windows[0].get_pos()[1]-Windows[0].get_pos()[0]
	print("number of widows:",len(Windows),"   windowSize:",windowSize)
	
	#initialise priors
	for windx,currentWin in enumerate(Windows):
		for ndx, node in enumerate(Windows[windx].getNodes()):
			node_id=node.get_ID()
			tmp_prior={}
			if is_calcPriors:
				relative_amps=node.get_relativeAmps()
				for p in Ploidy_order:
					tmp_prior[p]=relative_amps[p]
			else:
				for p in Ploidy_order:
					tmp_prior[p]=initial_prob
			node.set_EMpriors(tmp_prior)

	
	
	converged=False
	iter_tracker=0 # stores the currect ineration number for later logs
	convergerCount=0 # Tracks if all blocks have converged or not
	
	#Perform the First step. Calculate weighted contributions of kmers within each block to that block
	
	for iterationCount in range(1,iterMax+1,1):
		iter_tracker=iterationCount
	
		if convergerCount == num_nodes: #Check if all nodes have converged. If so, end the loop
			iter_tracker=iterationCount-1
			converged=True
			print("All nodes converged in iteration",iter_tracker)
			break

		convergerCount=0 # Tracks if all blocks have converged or not
		#print('priors')
		for windx,currentWin in enumerate(Windows):
			for ndx, node in enumerate(Windows[windx].getNodes()): #for each node in each window
				if not node.check_isConverged():
					tmp_prob_tally={} #define a temporary dict that will hold average
					tmp_prob_avg={} #define a temporary dict that will hold average
					tmp_prob_sum={}
					
					for p in Ploidy_order: #for each ploidy level
						tmp_prob_tally[p]=[] #Fill a dummy value 
						tmp_prob_sum[p]=0.0 #Fill a dummy value 
						tmp_prob_avg[p]=node.get_EMpriors()[p]
					#print(windx, node.get_ID(), tmp_prob_avg)
					#print('kmers, block',bdx)
	
					for kmer_count in sorted(node.get_KmerCounts()):#SPEEDUP- can drop the sorted here
						
						kmer_num=node.get_KmerCounts()[kmer_count] # get the number of times we have seen a kmer of this count for this node
						
						kmer_weight=1.0
						kmer_prob=get_probs_fromModel(kmer_count,Pars)
						#print(kmer_count, kmer_num,kmer_prob)
	
						for p in kmer_prob: #basically, for each ploidy level p
							tmp_contribution=kmer_prob[p]*tmp_prob_avg[p]*kmer_weight*kmer_num #the contribution = probability from the count X the probability of the block already X the kmer weighting X the number of times we've seen this kmer count
							#?tmp_prob_tally[p].append(tmp_contribution)#Make a list of all probabilities observed for this block at each ploidy level
							tmp_prob_sum[p]+=tmp_contribution #Fill a dummy value 
	
					node.set_Contributions(tmp_prob_sum)
	
	
	
		store_tmpContributions={} # store the unnormalised EM prior values for each block, only ipdate after they have all updated eachother (so that block 1 doesn't update block 0, and then get updated by the updated block 0)
		#print('Summing Weighted neighbour contributions and normalising contributions')
		for windx,currentWin in enumerate(Windows):
			for ndx, node in enumerate(Windows[windx].getNodes()): #for each node in each window
				if not node.check_isConverged():
					#initialise a tmp dictionary to store contributions from this node and weighted contributions from neighbouring nodes.
					neighbour_Contributions={}
					#for p in Ploidy_order: #for each ploidy level
					#	neighbour_Contributions[p]=0.0
					og_genomes=node.getGenomes()
					#Get influence from neighbours
					window_distance=int(neighbour_distance_max/windowSize) # how far up and downstream on the graph should we look?
					for idx in range(windx-window_distance,windx+window_distance,1):
						if idx >=0 and idx < len(Windows): #If the neighbouring nodes are within bounds of this haplotype
							if idx==windx: #If we are at the window with the node of interest
								localConts=node.get_Contributions() #store contributions of this node
								og_conts=localConts
								#for p in Ploidy_order:
								#	tmp_Contributions[p] += localConts[p] # update our contributions tracker with the unedited contributions which we calculated before
							
							else: #If this is a neighbouring window. TODO:TEST. UNTESTED!!!
								dist=Windows[idx].get_pos()[0]-currentWin.get_pos()[0]
								neighbour_Contributions[dist]=[] #this now holds list of contribution dictionaries (already intersection-weighted). keys are distance to node 
								#We want to check each node in this window, assess the overlap of their genomes with our current node, and update according to weight and overlap. 
								for mdx,m_node in enumerate(Windows[idx].getNodes()): #get nodes from this neighbouring window
									#calculate the overlap of genomes
									m_genomes=m_node.getGenomes()
									#Let the intesrect weight be 
									intersect_weight=(len(og_genomes.intersection(m_genomes))/len(og_genomes.union(m_genomes)))**10
									#print(windx, node.get_ID(),idx, m_node.get_ID(),og_genomes, m_genomes,len(og_genomes.intersection(m_genomes)),len(og_genomes.union(m_genomes)),intersect_weight)
	
									localConts=m_node.get_Contributions()
									
									#Apply distance function later in the piece, probably within the Node class.
									dist_weight=calc_weight_function_ExpDecay(dist, neighbour_distance_max,neighbour_weight_max)
									#print("Window",windx, "node",node.get_ID(),"will be updated by node", m_node.get_ID(),"in window",idx,dist,"bp distant, with", localConts, "with distance weight", dist_weight, "and intersect weight", intersect_weight)
									#for p in Ploidy_order:

									for key in localConts: # update contributions by the intersect weight (edge weight between this node and the node of interest)
										localConts[key]*=intersect_weight
										localConts[key]*=dist_weight
									#if node.get_ID()=="2232":
									#	print(dist, m_node.get_ID(), intersect_weight, dist_weight,localConts)
									local_sum=sum(localConts.values())
									if local_sum >1.0: #if sum effect of this node is greater than 1 kmers worth of contribution
										neighbour_Contributions[dist].append(localConts)
	
					node.set_NeighbourContributions(neighbour_Contributions)
				else:
					#print("window",windx, "node", node.get_ID(),"is converged")
					convergerCount+=1
	
		#See how the weighted neighbour contributions measure up to the nodes own contributions
		#for n in Windows[2].getNodes():
		#	print(n.get_ID())
		#	print(n.get_Contributions())
		#	print(n.get_NeighbourContributions())
		print('Iteration:',iterationCount, "  ",convergerCount,"/",num_nodes,"nodes converged","({0}%)".format(round(convergerCount/num_nodes*100,2)))

	
		for windx,currentWin in enumerate(Windows):
			for ndx, node in enumerate(Windows[windx].getNodes()): #for each node in each window
				if not node.check_isConverged():



					#release this later
					neighbour_conts=node.calc_neighbour_contributions(Ploidy_order)
					newContributions={}
					for p in Ploidy_order:
						#check if this node has no kmers, then ignore for now
						#newCont=np.log10(node.get_Contributions()[p])+np.log10(neighbour_conts[p]) # use log10 to soften impact of kmer differences.. butttt also softens differences between ploidies duh
						if node.get_totalMarkers() <1: # don't update empty nodes with neighbours
							newCont=node.get_Contributions()[p]
						else:
							newCont=node.get_Contributions()[p]+neighbour_conts[p]
						if newCont > 0.0000000001:
							newContributions[p]=newCont
						else:
							newContributions[p]=0.0
					

					newPriors= NormaliseDictionaryValues(copy.copy(newContributions))#can drop this copy later, just need for debugging now
					#Very handy debug code block
					#if node.get_ID()=="2232":
					#	print(node.get_ID(),round(node.get_covFraction(),2))
					#	print("contributions(self)",node.get_Contributions())
					#	print("contributions(neighbours)",neighbour_conts)
					#	print("contibutions(sum)",newContributions)
					#	print("newPriors",newPriors)

					node.set_NeighbourContributions(None)
					node.update_EMpriors_withChecks(newPriors,minUpdate,iterationCount,Ploidy_order)


	
	if not converged:
		print("Did not converge, stopped after",iter_tracker, "iterations")

	#output finale
	print("Final EM probs")
	out=open(outFile,"w+")
	for windx,currentWin in enumerate(Windows):
		for ndx, node in enumerate(Windows[windx].getNodes()): #for each node in each window
			EMPs=node.get_EMpriors()
			outEMPstr=""
			for e in EMPs:
				outEMPstr= outEMPstr+str(round(EMPs[e],2))+"\t"
			outEMPstr=outEMPstr[:-1]
			if is_calcPriors:
				RAs=node.get_relativeAmps()
				outRAstr=""
				for r in RAs:
					outRAstr= outRAstr+str(round(RAs[r],2))+"\t"
				outRAstr=outRAstr[:-1]
			else:
				outRAstr="\t".join([str(initial_prob),str(initial_prob),str(initial_prob),str(initial_prob),str(initial_prob)])
			out.write(str(currentWin.get_pos()[0])+"\t"+str(currentWin.get_pos()[1])+"\t"+str(node.get_ID())+"\t"+outRAstr+"\t"+str(node.check_isConverged())+"\t"+outEMPstr+"\n")


	
if __name__ == "__main__":
	print("\nEM.py v:"+str(version)+"\n")
	# this won't be run when imported
	#Parse arguments form the command line
	parser = argparse.ArgumentParser(description="EM.py  - Neighbour-aware Expectation Maximisation for nodes of a Haplotype Graph")
	subparsers = parser.add_subparsers(dest="command")

	

	#Parser for arguments when user calls command 'call'
	#call (graphFile,nodeStatsFile,modelsFile,minUpdate,iterMax,Neighbour_distance_max,neighbour_weight_max)
	parser_call = subparsers.add_parser('call')
	parser_call.add_argument('-G', '--GraphFile', dest='graphFile', help="PanToHap Graph File (haplotype_graph.txt)", required=True)
	parser_call.add_argument('-N', '--NodeFile', dest='nodeStatsFile', help="PanToHap Node stats File (node_k_stats.txt)", required=True)
	parser_call.add_argument('-M', '--Model', dest='modelsFile', help="The tsv file output from fit_peaks fit command, with parameters of models explaining short-read coverage", required=True)

	parser_call.add_argument('-o', '--outFile', dest='outFile',type=str, default="/EM_results.tsv", help="full path to  tsv file, including file name (/EM_results.tsv)", required=False)
	parser_call.add_argument('--nodePlots',dest='is_nodePlots',action='store_true')
	parser_call.add_argument('-p', '--plotFolder', dest='plotFolder',type=str, default="/EM_results.tsv", help="full path to  directory where you want node plots to be saved", required=False)
	parser_call.add_argument('--CalcPriors',dest='is_calcPriors',action='store_true')

	parser_call.add_argument('-u', '--minUpdate', dest='minUpdate', type=float, default=0.001, help="The minimum EM-step update before we say this block has converged (0.00001)", required=False)
	parser_call.add_argument('-i', '--maxIterations', dest='iterMax', type=int, default=100, help="The maximum niber of EM steps before stopping(100)", required=False)

	#parser_call.add_argument('-w', '--windowSize_Kb', dest='windowSize_Kb', type=int, help="The size of the windows to break genomes into (in Kilobases)", required=True) # now extracts from graph file - assumes all are same size
	parser_call.add_argument('-d', '--maxNeighbourDist', dest='neighbour_distance_max', type=int, default=600000, help="the distance from which an adjacent block will influence the descisions of another 0.5 as much as itself (essentially the 50% reduction point of the LD-decay), default is 600Kb  (600000)", required=False)
	parser_call.add_argument('-w', '--neighbourWeightMax', dest='neighbour_weight_max', type=float, default=1.0, help="the maximum influence that an an adjacent block can exert on the decisions of another 0.0-1.0 (1.0)", required=False)
	#parser_call.add_argument('-p', '--minProb', dest='min_prob', type=float, help="The probabilty that a kmer not present on an assembled genome did in fact come from that genome (0.01)", required=True)

	parser_call.add_argument('-K', '--maxK', dest='maxK', type=int, default=1000, help="the maximum kmer count that will be used for updating EM probabilities, recommend to set this at the tail end of your 4-ploid kmer peaks (1000)", required=False)


	#Parse arguments
	kwargs = vars(parser.parse_args())
	#print(kwargs)
	command = kwargs.pop('command')

	if not 0:
		# call appropriate functions
		globals()[command](**kwargs)
