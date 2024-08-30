#Script to generate relative probabilities of different kmers coming from a haplotype of 0/1/2/3/4-copy number
version=0.7
#v7 fixes standard deviation for 2,3,4-copy peaks as 1/sqrt(n) of 1-copy peak
import sys
import time

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
import math
import pandas as pd
import argparse




def _1gaussian(x_array, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))

def _subtractGaussian(x_array, y_array, amp1,cen1,sigma1):
    return y_array - (amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))))

def _4gaussian(x_array, amp1,cen1,sigma1, amp2,cen2,sigma2, amp3,cen3,sigma3, amp4,cen4,sigma4):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen2)/sigma2)**2))) + \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen3)/sigma3)**2))) + \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen4)/sigma4)**2)))


def _4gaussian_fixed(x_array, amp1,cen1,sigma1,amp2,sigma2,amp3,sigma3,amp4,sigma4):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*2))/sigma2)**2))) + \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*3))/sigma3)**2))) + \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*4))/sigma4)**2)))


def _4gaussian_fixed_plusError(x_array, k2,b2, amp1,cen1,sigma1,amp2,sigma2,amp3,sigma3,amp4,sigma4):
    return [k2*(b2**(x-0)) for x in x_array] + \
    		amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*2))/sigma2)**2))) + \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*3))/sigma3)**2))) + \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*4))/sigma4)**2)))
     #0 replaces c2


def _4gaussian_fixedCen_plusError_fixedSigma(x_array, k2,b2, amp1,cen1,sigma1,amp2,amp3,amp4):
	cen2=cen1*2
	cen3=cen1*3
	cen4=cen1*4
	sigma2=(sigma1*2)/np.sqrt(2)
	sigma3=(sigma1*3)/np.sqrt(3)
	sigma4=(sigma1*4)/np.sqrt(4)
	return [k2*(b2**(x-0)) for x in x_array] + \
	amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
	amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*2))/sigma2)**2))) + \
	amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*3))/sigma3)**2))) + \
	amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-(cen1*4))/sigma4)**2)))
     #0 replaces c2



def _1exponential(x_array,k2,b2,c2):
    return [k2*(b2**(x-c2)) for x in x_array]


def prob_display(modelsFile, xmin):
	d=pd.DataFrame()
	x_vals=[] #list of x values for plotting
	y_vals={} #dictionary that will hold lists of y values
	min_prob=0.0001
	ploidy_levels=set()

	#Read in parameters from modelsFile
	if(os.path.isfile(modelsFile)):
		print(modelsFile)
	else:
		print("models file not found")
		return 0

	Pars = {} #dictionary of parameters for models
	for ldx, l in enumerate(open(modelsFile,"r")):
		print("here",l)
		if ldx>0: # skip header
			vals=l.rstrip().split("\t")
			ploidy=int(vals[0])
			ploidy_levels.add(ploidy)
			Pars[ploidy]=[float(vals[1]),float(vals[2]),float(vals[3])]

	ploidy_order=sorted(ploidy_levels)
	print(ploidy_order)
	xmax=int(Pars[1][1])*(max(ploidy_order)+3)
	print("xmax:",xmax)
	outTSV=open(modelsFile+".prob_display.tsv","w+")
	outTSV.write("kmer_ID"+"\t"+"\t".join([str(x) for x in ploidy_order])+"\n")

	for kmer_count in range(xmin,xmax+1,1):
		x_vals.append(kmer_count)
		kmer_ploidy_probs={}
		for p in ploidy_order:
			if p==0: #if we are looking at a noise model
				pars_tmp= Pars[p]
				#kmer_ploidy_probs[p]= _1reciprocal(int(kmer_count),*pars_tmp)
				kmer_ploidy_probs[p]= _1exponential([int(kmer_count)],*pars_tmp)[0] # save as list to satisfy list comprehension in function
				#print(p, kmer_ploidy_probs)
				if kmer_ploidy_probs[p]<min_prob: #if this model is predicting less than 1 kmer, then lets set an error baseline of 1
					kmer_ploidy_probs[p] = 1.0
			elif p>0:
				pars_tmp= Pars[p]
				kmer_ploidy_probs[p]= _1gaussian(int(kmer_count),*pars_tmp)
				if kmer_ploidy_probs[p]<min_prob: #if this model is predicting less than 1 kmer, then lets consider it to be basically 0
					kmer_ploidy_probs[p] = min_prob
		#print(kmer_ploidy_probs)

		prob_sum = sum(kmer_ploidy_probs.values())
		outString=str(kmer_count)
		#normalise by sum of models at this position to get relative probabilities
		for p in ploidy_order:
			y=kmer_ploidy_probs[p]/prob_sum
			if p not in y_vals:
				y_vals[p]=[]
			y_vals[p].append(y)
			outString=outString+"\t"+f"{y:.5f}"
		outTSV.write(outString+"\n")


	#Plot all relative prob values
	plotmax=xmax
	colours=["grey","blue","orange","green","purple"]
	fig = plt.figure(figsize=(10,3))
	gs = gridspec.GridSpec(1,1)
	ax1 = fig.add_subplot(gs[0])
	ax1.set_xlim([0,plotmax])
	for adx,a in enumerate(y_vals):
		ax1.plot(x_vals, y_vals[a], color=colours[adx], marker='o', linestyle='solid', linewidth=2, markersize=5)
	fig.tight_layout()
	plt.savefig(modelsFile+".prob_display.pdf", transparent=False, dpi=300, format="pdf")
	#plt.show()



def fit(histogram_txt, ploidy):
	outTSV=open(histogram_txt+".gaussianFit.tsv","w+")
	df=pd.read_csv(histogram_txt,sep=" ",header=None)
	df.columns=['count','frequency']
	print(df)

	xseries=df.iloc[:,0:1].max()
	xmax=xseries.max()
	xseries=df.iloc[:,0:1].min()
	xmin=xseries.min()

	ymax_index=df.iloc[3:,1:2].idxmax().max()
	est_cov=int(df.iloc[ymax_index,0:1])
	print("estimated haploid coverage:",est_cov)
	xmax=est_cov*7
	#create new list with corresponding values
	x_array = list(range(xmin,xmax+1,1))
	#create y array with corresponsing values or zeros
	y_array=[]
	for x in x_array:
		try:
			idx=df.loc[df['count'] == x].index[0]
			y_array.append(df.iloc[idx]['frequency'])
		except Exception:
			#print(x,"0")
			y_array.append(0)

	gauss="4_fixedCen_plusError_fixedSigma"

	#Estimate params for haplomer peak using a single gaussian curve.
	max_yish = max(y_array[3:])*10
	est_cov=est_cov
	sigma_0 = est_cov/4 
	print(max_yish, est_cov, sigma_0)
	
	popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_1gaussian, x_array[est_cov-int(est_cov/2):est_cov+int(est_cov/2)], y_array[est_cov-int(est_cov/2):est_cov+int(est_cov/2)], p0=[max_yish, est_cov, sigma_0]) #ignore zero-three for fitting the first peak
	
		
	pars_1 = popt_gauss[0:3]
	gauss_peak_1 = _1gaussian(x_array, *pars_1)
	#Plot tmp1
	plotmax=est_cov*7

	fig = plt.figure(figsize=(10,3))
	gs = gridspec.GridSpec(1,1)
	ax1 = fig.add_subplot(gs[0])
	ax1.set_xlim([0,plotmax])
	ax1.set_ylim([0,max(y_array[3:])*1.1])
	ax1.plot(x_array, y_array, "ro")
	ax1.plot(x_array, gauss_peak_1)
	ax1.text(est_cov*3,max(y_array[3:]),"check to see if curve roughly fits haploid peak",size=10)
	fig.tight_layout()
	plt.savefig(inFile+".hap_check.pdf", transparent=False, dpi=300, format="pdf")
	#plt.show()
	#set up for actual run
	guess1=pars_1

	if gauss=="4_fixedCen_plusError_fixedSigma":
		max_yish = guess1[0]
		est_cov=est_cov
		sigma_0 = guess1[2]*0.8
		y_int=y_array[0]
		exp_b=0.3
		print(max_yish, est_cov, sigma_0)
	
		popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_4gaussian_fixedCen_plusError_fixedSigma, x_array[0:], y_array[0:], p0=[ y_int, exp_b, max_yish, est_cov, sigma_0, max_yish/3, max_yish/15,max_yish/50])

		print(popt_gauss)
		print(pcov_gauss)
		pars_0 = popt_gauss[0:2]
		pars_1 = popt_gauss[2:5]
		pars_2 = [popt_gauss[5],(popt_gauss[3]*2),(popt_gauss[4]*2)/np.sqrt(2)]
		pars_3 = [popt_gauss[6],(popt_gauss[3]*3),(popt_gauss[4]*3)/np.sqrt(3)]
		pars_4 = [popt_gauss[7],(popt_gauss[3]*4),(popt_gauss[4]*4)/np.sqrt(4)]

		gauss_peak_4_fixed = _4gaussian_fixedCen_plusError_fixedSigma(x_array, *popt_gauss)
		#Plot tmp1
		plotmax=int(pars_1[1])*7
		print(plotmax)
		fig = plt.figure(figsize=(10,6))
		gs = gridspec.GridSpec(1,1)
		ax1 = fig.add_subplot(gs[0])
		ax1.set_xlim([0,plotmax])
		ax1.set_ylim([0,max(y_array[4:])*1.2])

		ax1.plot(x_array, y_array, "ro")
		ax1.plot(x_array, gauss_peak_4_fixed,lw=1, color="k")

		xplot=list(range(0,plotmax,1))
		noise_decay= _1exponential(xplot,*pars_0,0.0)
		gauss_peak_hap= _1gaussian(xplot, *pars_1)
		gauss_peak_dip= _1gaussian(xplot, *pars_2)
		gauss_peak_trip=_1gaussian(xplot, *pars_3)
		gauss_peak_tet= _1gaussian(xplot, *pars_4)
		#overcount= [1.0]*len(xplot)

		colours=["grey","blue","orange","green","purple"]
		ax1.plot(xplot, noise_decay,color=colours[0],ls="--",lw=2)
		ax1.plot(xplot, gauss_peak_hap,color=colours[1],ls="--",lw=2)
		ax1.plot(xplot, gauss_peak_dip,color=colours[2],ls="--",lw=2)
		ax1.plot(xplot, gauss_peak_trip,color=colours[3],ls="--",lw=2)
		ax1.plot(xplot, gauss_peak_tet,color=colours[4],ls="--",lw=2)

		fig.tight_layout()
		plt.savefig(inFile+".gaussianFit.pdf", transparent=False, dpi=300, format="pdf")

		#plt.show()


		#Write the output
		outTSV.write("ploidy\tamplitude/k\tcentre/b\tsigma/c\n")
		outTSV.write(("0\t"+"\t".join([f"{x:.4f}" for x in pars_0])+"\t0.0"+"\n"))
		outTSV.write(("1\t"+"\t".join([f"{x:.4f}" for x in pars_1])+"\n"))
		outTSV.write(("2\t"+"\t".join([f"{x:.4f}" for x in pars_2])+"\n"))
		outTSV.write(("3\t"+"\t".join([f"{x:.4f}" for x in pars_3])+"\n"))
		outTSV.write(("4\t"+"\t".join([f"{x:.4f}" for x in pars_4])+"\n"))
		#outTSV.write(("5+\t"+"\t".join([f"{x:.4f}" for x in pars_5plus])+"\n"))



inFile=sys.argv[1] # see Analysis.sh
ploidy=int(sys.argv[2]) # should be 4 for tetraploid potato

fit(inFile, ploidy)
prob_display(inFile+".gaussianFit.tsv",0)