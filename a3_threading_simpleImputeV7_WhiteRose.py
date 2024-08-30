#simpleImputeV2 fixes bug with Subtractions not being included for imputes when checking if ploidys are same either side of impute.
#Fixed bug where last_ploidy was being used for checks instead of this_ploidy. Meaning threads couldn't meaningfully pass through ploidy changes.
#V3 tracks the most seen genomes for each path and won't follow forks unless the main players (most seen nodes) are carried through in thread. 
#V4 fixes bug where V3 was dipping into (subtracted) empty paths
#V5 treats nodes without informative kmers differently from nodes with EM calls as 0-copy. This should relax the imputation a bit and maybe lead to higher coverage.
#   So it draws blind nodes differently
#   Also gets more conservative around forks, to avoid haplotype switches.
#V6 considers inversions and calculates a second accuracy score accordingly
#V7 fixes big haplotype switches by adding a heuristic to only follow an edge (even a single edge) if most of the most-seen genomes are going that way.
    #After some burn-in period (2 nodes) where all of the seen genomes must be the same. 
    #Also allows imputation at end of thread recursion if nothing further was found
    #Also fixes bug in haplotype accuracy calling where two threads could both have a max gneome of A_hap3. Meaning one would be called a hap-switch error when its really not.

import sys
from collections import deque 
from collections import Counter
from itertools import product
from itertools import chain
import copy
import string
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import random
import ast

emoutfin=sys.argv[1]
hapfin=sys.argv[2]
outPrefix=sys.argv[3]
#chrsize = 46102915
pericentromere=[int(x) for x in sys.argv[4].rstrip().split(",")]
print(pericentromere)
chrom=sys.argv[5]
inversionsFile=sys.argv[6]
targetHaps=["A_hap1","A_hap2","A_hap3","A_hap4"]

ReferenceLines=[]
for t in targetHaps:
    ReferenceLines.append([])


samples = ['dm'] + [f'{i}_hap{j}' for i, j in product(string.ascii_uppercase[:10], range(1, 5))]

#read in inversions file
inversions={} #dictionary, each key contains a list of tuples (dm_start,dm_end)
for line in open(inversionsFile,"r"):
        if line.split("/")[0]==chrom: # if we are looking a the right chromosome
            vals=line.rstrip().split("\t")
            hap=vals[0].split("/dm_")[1].split("_")[0]+"_"+vals[0].split("/dm_")[1].split("_")[2].split("syri")[0]
            dm_start=vals[1]
            dm_end=vals[2]
            if hap not in inversions:
                inversions[hap]=[]
            inversions[hap].append((dm_start,dm_end))
inversions["dm"]=[] # dummy value because dm can't be inverted relative to itself

outFile_log=open(outPrefix+".log","w+")
outFile_threads=open(outPrefix+".threads.tsv","w+")
outFile_accuracy=open(outPrefix+".WhRaccuracy.tsv","w+")
outFile_WindowCov=open(outPrefix+".windowCov.tsv","w+")


def plot_threads():
    import igraph as ig
    def find_all_paths(graph, start, end):
        def dfs(current, path):
            path.append(current)
            if current == end:
                all_paths.append(list(path))
            else:
                for neighbor in graph[current]:
                    if neighbor not in path:
                        dfs(neighbor, path)
            path.pop()

        all_paths = []
        dfs(start, [])
        return all_paths

    class hapobject:
        """
        A haplotype block object
        """
        def __init__(self, id, start, end, genomes):
            self.id = id
            self.start = start
            self.end = end
            self.genomes = genomes
        # END

        def hasgen(self, gen):
            return gen in self.genomes
        # END
    # END

    Positions=[]
    RefPoints=[]
    scale=4
    hapoblist = deque()
    lastWindow=1-100000
    numWindows=0
    addedMissingNode=False
    WindowsPos={} #dictionary holding which nodes are in which windows
    with open(hapfin, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            startpos=int(line[0])
            if int(startpos/100000) not in WindowsPos:
                WindowsPos[int(startpos/100000)]=[]
            if startpos != lastWindow: # If we are updating window 
                while startpos != lastWindow+100000:
                    print("filling dummy value for pos,",lastWindow+100000)
                    print("filling dummy value for pos,",lastWindow+100000,file=outFile_log)
                    for tdx, t in enumerate(targetHaps):
                        ReferenceLines[tdx].append((tdx*scale)-0.4-(0.22*tdx)) #fill a dummy
                    lastWindow=lastWindow+100000
                    numWindows+=1
            endpos=int(line[1])
            this_haps=sorted(line[2].split(','))
            hapoblist.append(hapobject(int(line[3]), startpos, endpos, this_haps))
            ypos=5+(random.randint(0,20)/6) #10 rows "other" row of nodes
            for tdx,target in enumerate(targetHaps): #for each target haplotype
                if target in this_haps: 
                    ypos=tdx
                    break
            WindowsPos[int(startpos/100000)].append(len(Positions))
            Positions.append([(startpos/100000),ypos*scale])
            if "dm" in this_haps:
                RefPoints.append([startpos/100000,ypos*scale,1])
            else:
                RefPoints.append([startpos/100000,ypos*scale,0])
            RefFound=False
            mintx=5
            for tdx,target in enumerate(targetHaps): #for each target haplotype
                if target in this_haps: 
                    if tdx<mintx: #updated minimum target seen in thise node(in terms of index in list)
                        mintx=tdx
                    ReferenceLines[tdx].append((mintx*scale)-0.4-(0.22*tdx))
            if startpos != lastWindow:
                lastWindow=startpos
                numWindows+=1
    print(numWindows, "windows")
    print(numWindows, "windows",file=outFile_log)

    # Create a Graph object for the haplotype blocks and use it to determine the Y-coordinate for visualisation of the blocks
    #print(Positions)
    G = ig.Graph()
    G = ig.Graph.as_undirected(G) #changed for plotting
    G.add_vertices(range(len(hapoblist)))

    addededge = set()
    startnodes = deque()
    for hap in samples:
        hc = -1
        for h in hapoblist:
            if h.hasgen(hap):
                if hc == -1:
                    hc = h.id
                    startnodes.append(h.id)
                    continue
                if (hc, h.id) not in addededge:
                    G.add_edge(hc, h.id)
                    addededge.add((hc, h.id))
                hc = h.id

    # Add ploidy levels in the graph
    ploidy = deque()
    emout = pd.read_csv(emoutfin, header=None, sep='\t')
    def get_ploidy(row):
        #print(row[9:])
        if all([r == 0 for r in row[9:]]):
            return -2
        # TODO: for a node, I assign based on the maximum probablity from EM.
        #  So, if for a node, EM outputs probablity of 51% for 1 copy and 49%
        #  for 2 copy, I assign ploidy of 1 to the node. Probably, this can be
        #  fine-tuned.
        return row[9:].to_list().index(max(row[9:]))
    # END

    ploidy = emout.apply(get_ploidy, axis=1)
    #print(ploidy)
    G.vs['ploidy'] = ploidy

    names= []
    for h in hapoblist:
        names.append(h.id)
    G.vs['name'] = names
    genomes=[]
    for h in hapoblist:
        genomes.append(h.genomes)
    G.vs['genomes']=genomes

    inv_list=[]
    for hdx,h in enumerate(G.vs['genomes']):
        #print(hdx, h)
        inv_cov='circle'
        for g in h:#for each genome in genome list for this node
            for i in inversions[g]: # for each tuple of inversion boundary positions associated with this haplotype
                if Positions[hdx][0]*100000 >= int(i[0]) and (Positions[hdx][0]*100000)+100000 <= int(i[1]): #if this window fits inside the inversion
                    #print(g,Positions[hdx][0]*100000,(Positions[hdx][0]*100000)+100000,i[0],i[1])
                    inv_cov='triangle-down'
        inv_list.append(inv_cov)
    G.vs['inversionShape']=inv_list

    # Delete edges between nodes with unequal ploidy, and edges to blind nodes, and 
    to_del = [e for e in G.get_edgelist() if (G.vs[e[0]]['ploidy'] != G.vs[e[1]]['ploidy'] and (G.vs[e[0]]['ploidy']==0 or G.vs[e[1]]['ploidy']== 0)) or (G.vs[e[0]]['ploidy']==0) or (G.vs[e[1]]['ploidy']== 0) or (G.vs[e[0]]['ploidy']==-2) or (G.vs[e[1]]['ploidy']== -2)]
    #print(to_del)
    G.delete_edges(to_del)
    #print(G.get_edgelist())
    #Nodes don't change! just edges go

    # Check whether the predicted ploidy is correct
    ploidy_color = {
        -3: "darkgrey",
        -2: "lightgrey",
        -1: "white",
        0: "black",
        1: "orange",
        2: "blue",
        3: "red",
        4: "green"
    }
    cor_plo = deque()   
    for v in G.vs:
       p = v['ploidy']
       gs = vars(hapoblist[v['name']])['genomes']
       if p == len([g for g in gs if 'A_' in g ]):
           cor_plo.append(ploidy_color[v['ploidy']])
       elif p==-2:
           cor_plo.append(ploidy_color[-1])
           #v['ploidy'] = 0
       else:
        if p==0:
            cor_plo.append(ploidy_color[-3])
        else:
            cor_plo.append(ploidy_color[-2])

    G.vs['cor'] = cor_plo

    # Save artsy figure
    #ig.plot(G, target=outPrefix+'.tmp.pdf', layout='kk', vertex_size=[p+4 for p in G.vs['ploidy']], vertex_color=G.vs['cor'], vertex_frame_color="black", vertex_frame_width=0.1,  edge_width=0.1, edge_arrow_size=0.2, bbox=[800, 800])

    #Try get it in genomic coordinates
    chromlayout=G.layout("grid")
    chromlayout=Positions

    #get genric plot output
    ig.plot(G, target='tmp.grid.pdf',layout=chromlayout, vertex_size=np.array(G.degree())+4, vertex_color=G.vs['cor'], vertex_frame_color="black", vertex_frame_width=0,  edge_width=0.1, edge_arrow_size=0.2, bbox=[2000, 300])

    #Get output with reflines
    fig,ax=plt.subplots()
    fig.set_figwidth(60)
    fig.set_figheight(5)
    ax.set_xlim(-2,numWindows+2)
    ax.set_ylim(-2,9*scale)
    #plt.hlines(1,0,230)

    #get pericentromere block
    peri_rect=plt.Rectangle(((pericentromere[0]/100000)-0.5,-1),pericentromere[1]/100000-pericentromere[0]/100000,9*scale)
    coll=[peri_rect]
    pc=PatchCollection(coll,facecolor="green",alpha=0.08,zorder=-2)
    ax.add_collection(pc)
    
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[0],color="green",zorder=-1)
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[1],color="purple",zorder=-1)
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[2],color="firebrick",zorder=-1)
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[3],color="pink",zorder=-1)

    ig.plot(G, target=ax,layout=chromlayout, vertex_shape=G.vs['inversionShape'],vertex_size=[p+4 for p in G.vs['ploidy']], vertex_color=G.vs['cor'], vertex_frame_color="black", vertex_frame_width=0.1,  edge_width=0.1, edge_arrow_size=0.0, bbox=[3000, 300])


    #Final touches and plot
    ax.set_xticks(list(range(0,numWindows+1,1)))
    ax.set_xticklabels(list(range(0,numWindows+1,1)),rotation=90, size=7)
    plt.tight_layout()
    fig.savefig(outPrefix+".gridRefs.pdf",dpi=300) # plot without threads
    #plt.show()



    #Now try finding accuracy
    totalCall=0
    correctCall=0
    for window in chain(range(0,int(pericentromere[0]/100000),1),range(int(pericentromere[1]/100000),numWindows+1,1)):
        if window >=0 and (window*100000<=pericentromere[0]-1 or window*100000 >=pericentromere[1]): 
            copiesCalled=0
            for pdx,p in enumerate(Positions):
                #print(p)
                if int(p[0])==window:
                    copiesCalled+=G.vs['ploidy'][pdx]

            #print(window,":",copiesCalled)
            if copiesCalled==4:
                correctCall+=1
            totalCall+=1

    
    outFile_log.write(str(numWindows)+" windows\n")
    tmp_str=str(correctCall/totalCall)+" non-pericentromeric windows called as 4-copy from raw EM output: "+str(correctCall)+" out of "+str(totalCall)
    print(tmp_str)
    outFile_log.write(tmp_str+"\n")
    tmp_str=str(totalCall/10*4)+" Mb to cover, across 4 haps, outside of pericentromere"
    print(tmp_str)
    outFile_log.write(tmp_str+"\n")
    
    #Now try finding N50s

    Subtractions=[0]*len(Positions)
    WindowsThreadCalls=[0]*(numWindows+1)

    #Fancy lookahead to increase edge count if there is a likely path forward for one of the most seen genomes that is blocked only by a blind node.
    artificialEdgeTracker={} # a dictionary holding all genomes for a given node which connect to a blind node - later these will be assessed if they are a most-seen genome for a thread.
    for currentVtx, p in enumerate(Positions): #currentVtx now acts as our vertex reference number
        artificialEdgeTracker[currentVtx]=[]
        for w in WindowsPos:
            if currentVtx in WindowsPos[w]:
                this_window=int(w)
        nextWindow=this_window+1
        if nextWindow in WindowsPos:
            for ndx in WindowsPos[nextWindow]:
                relevant=False
                for g in G.vs['genomes'][currentVtx]:
                    if g in G.vs['genomes'][ndx]: #if this blind node also contains a genome in the current node, then there is a connection from this node to a blind node
                        if G.vs['ploidy'][ndx]==-2:#if this is a blind node
                            artificialEdgeTracker[currentVtx].append(g)
                        
    #Function for recursive thread searching
    def maxPathLenfindMaxPath1PloidyChange(currentVtx,lastPloidy,ploidchangesSoFar,imputesSoFar,pathSoFar,SeenGenomes,EdgeList):
        #print(currentVtx,pathSoFar)
        this_window=int(Positions[currentVtx][0])
        this_Ploidy=int(G.vs['ploidy'][currentVtx]) #get ploidy of the current vertex, raw from EM
        this_edgeList=[]
        longestLength=len(pathSoFar)
        longestPath=pathSoFar
        longestPloidyChanges=ploidchangesSoFar
        longestPathImputes=imputesSoFar

        this_Genomes=G.vs['genomes'][currentVtx]
        for g in this_Genomes:
            SeenGenomes[g]+=1

        maxSighting=max(SeenGenomes.values())
        mostSeenGenomes=[]
        for s in SeenGenomes:
            if SeenGenomes[s]>=0.6*maxSighting:
                mostSeenGenomes.append(s)

        #Break if we enter a node with 0 ploidy remaining
        if (int(this_window)>=int(pericentromere[0]/100000) and int(this_window)<=int(pericentromere[1]/100000)) or WindowsThreadCalls[this_window]==4 or this_Ploidy-Subtractions[currentVtx] ==0: #if we've entered the pericentromere, or saturated this window, or entered a 0-state
            return True,longestLength-1,longestPath[:-1],ploidchangesSoFar,longestPathImputes #return the previous values, basically cut short the search
        #Break if we change genomes in the first 3 nodes of the thread
        if len(pathSoFar)<2:# and not all(x in mostSeenGenomes for x in this_Genomes):
            for t in this_Genomes:
                if t not in mostSeenGenomes:
                    return True,longestLength-1,longestPath[:-1],ploidchangesSoFar,longestPathImputes #return the previous values, basically cut short the search
        
        
        #Fancy lookahead to increase edge count if there is a likely path forward for one of the most seen genomes that is blocked only by a blind node.
        artificialEdgeCount=0
        blindConnection=False
        for m in mostSeenGenomes:
            if m in artificialEdgeTracker[currentVtx]:
                blindConnection=True
        if blindConnection:
            artificialEdgeCount+=1

        edgeCount=0
        for e in EdgeList:
            if int(e[0])==currentVtx: #if theres is an edge starting with this node
                tpar=int(e[1]) #get the partner
                if int(G.vs['ploidy'][tpar])-Subtractions[tpar] >0:
                    this_edgeList.append(tpar)
                edgeCount+=1
        #for t in this_edgeList:
        #recursively follow all edges starting from this node, to find maximum path with no more than 1 ploidy change
        for t in this_edgeList: 
            if edgeCount==1 and edgeCount+artificialEdgeCount==1 or all([m in G.vs['genomes'][t] for m in mostSeenGenomes]): #len(set(G.vs['genomes'][currentVtx]).intersection(set(G.vs['genomes'][t])))/len(set(G.vs['genomes'][currentVtx]))>0.5:
                if int(G.vs['ploidy'][t]) ==this_Ploidy and int(G.vs['ploidy'][t])-Subtractions[t]!=0: # if we are not passing through a ploidy change
                    #print(currentVtx,t)
                    windowFull,maxPathLength,maxPath,maxPloidyChanges,maxImputes=maxPathLenfindMaxPath1PloidyChange(t,this_Ploidy,ploidchangesSoFar,imputesSoFar,pathSoFar+[str(t)],SeenGenomes.copy(),EdgeList)
                    if not windowFull and maxPathLength > longestLength:
                        longestLength=maxPathLength
                        longestPath=maxPath
                        longestPloidyChanges=maxPloidyChanges
                        longestPathImputes=maxImputes
                elif int(G.vs['ploidy'][t]) !=this_Ploidy and int(G.vs['ploidy'][t])-Subtractions[t]!=0 and ploidchangesSoFar==0: # else if there is a difference and it does not equal 0
                    windowFull,maxPathLength,maxPath,maxPloidyChanges,maxImputes=maxPathLenfindMaxPath1PloidyChange(t,this_Ploidy,ploidchangesSoFar+1,imputesSoFar,pathSoFar+[str(t)],SeenGenomes.copy(),EdgeList)
                    if not windowFull and maxPathLength > longestLength:
                        longestLength=maxPathLength
                        longestPath=maxPath
                        longestPloidyChanges=maxPloidyChanges  
                        longestPathImputes=maxImputes
        if len(this_edgeList)==0 or longestLength==len(pathSoFar): #if we didn't find anything, then look 2 windows ahead.
            ubernextWindow=this_window+2
            for udx, u in enumerate(Positions):
                if ubernextWindow==int(u[0]): #if this node is in this window
                    if len(set(G.vs['genomes'][currentVtx]).intersection(set(G.vs['genomes'][udx])))/len(set(G.vs['genomes'][currentVtx]).union(set(G.vs['genomes'][udx])))>=0.8: #if they share 90% of their haplotypes
                        if G.vs['ploidy'][udx]==G.vs['ploidy'][currentVtx] and G.vs['ploidy'][udx]-Subtractions[udx]!=0: #if nodes are same ploidy
                            for pdx, p in enumerate(Positions):
                                if int(p[0])==this_window+1: # if this is an intermediate node
                                    #If this node shares enough haplotypes to be imputed through
                                    if len(set(G.vs['genomes'][currentVtx]).intersection(set(G.vs['genomes'][pdx])))/len(set(G.vs['genomes'][currentVtx]))>0.5:
                                        tmpPath=pathSoFar+["imputed_"+str(pdx)]
                                        windowFull,maxPathLength,maxPath,maxPloidyChanges,maxImputes=maxPathLenfindMaxPath1PloidyChange(udx,this_Ploidy,ploidchangesSoFar,imputesSoFar+1,tmpPath+[str(udx)],SeenGenomes.copy(),EdgeList)
                                        if not windowFull and maxPathLength > longestLength:
                                            longestLength=maxPathLength
                                            longestPath=maxPath
                                            longestPloidyChanges=maxPloidyChanges  
                                            longestPathImputes=maxImputes
            
    
        return False, longestLength,longestPath,longestPloidyChanges,longestPathImputes

    #Start the recursive thread search
    EdgeList=G.get_edgelist()
    outer=[]
    lastThreadlenSum=-1
    threadLenSum=0
    threadCounter=0
    lastMax=0
    lasthits={}
    while lastThreadlenSum!=threadLenSum:

        maxmaxPath=[]
        maxmaxWindowStart=0
        maxmaxPloidyChanges=0
        maxmaxPathImputes=0

        for currentVtx, p in enumerate(Positions): #currentVtx now acts as our vertex reference number
            #print(currentVtx)
            if currentVtx not in lasthits:
                lasthits[currentVtx]=0 #dummy value to track size of best thread seen at node last iteration
            #try a little speedup assuming not a big drop in thread sizes 
            if lasthits[currentVtx] >=lastMax*0.5 or lastMax <4: # SKIPPING FOR NOW
    
                this_window=int(p[0])
                this_Ploidy=int(G.vs['ploidy'][currentVtx])
                lastPloidy=this_Ploidy #dummy variable to start
                ploidchangesSoFar=0
                imputesSoFar=0
                pathSoFar=[str(currentVtx)]
                longestLength=1
                longestPath=pathSoFar
                longestPloidyChanges=0
                longestPathImputes=0
                #print(G.vs['genomes'][currentVtx])
                SeenGenomes=Counter()
                this_Genomes=G.vs['genomes'][currentVtx]
                for g in this_Genomes:
                    SeenGenomes[g]+=1
                
                maxSighting=max(SeenGenomes.values())
                mostSeenGenomes=[]
                for s in SeenGenomes:
                    if SeenGenomes[s]==maxSighting:
                        mostSeenGenomes.append(s)
    
                if this_Ploidy-Subtractions[currentVtx] >0 and (int(this_window)<=int(pericentromere[0]/100000) or int(this_window)>=int(pericentromere[1]/100000)) and WindowsThreadCalls[this_window]<4: #int(this_window)==284 #only thake things right of the pericentrome
                    #get all edges starting from this node
                    this_edgeList=[]
                    edgeCount=0
                    for e in EdgeList:
                        if int(e[0])==currentVtx: #if theres is an edge starting with this node
                            tpar=int(e[1]) #get the partner
                            if int(G.vs['ploidy'][tpar])-Subtractions[tpar] >0:
                                this_edgeList.append(tpar)
                            edgeCount+=1
                    #recursively follow all edges starting from this node, to find maximum path with no more than 1 ploidy change
                    for t in this_edgeList:#for each possible partner
                        if edgeCount==1 or all([m in G.vs['genomes'][t] for m in mostSeenGenomes]): #len(set(G.vs['genomes'][currentVtx]).intersection(set(G.vs['genomes'][t])))/len(set(G.vs['genomes'][currentVtx]))>0.5:
                            if int(G.vs['ploidy'][t]) ==this_Ploidy and int(G.vs['ploidy'][t])-Subtractions[t]!=0: # if we have no ploidy change from the last pick
                                #print(currentVtx,t)
                                windowFull,maxPathLength,maxPath,maxPloidyChanges,maxImputes=maxPathLenfindMaxPath1PloidyChange(t,this_Ploidy,ploidchangesSoFar,imputesSoFar,pathSoFar+[str(t)],SeenGenomes.copy(),EdgeList)
                                if not windowFull and maxPathLength > longestLength:
                                    longestLength=maxPathLength
                                    longestPath=maxPath
                                    longestPloidyChanges=maxPloidyChanges
                                    longestPathImputes=maxImputes
                            #could replace this lookahead elif with a simple if this:ploidy != lastPloidy inside the recursion.
                            elif int(G.vs['ploidy'][t])!=this_Ploidy and int(G.vs['ploidy'][t])-Subtractions[t]!=0 and ploidchangesSoFar==0: # else if there is a difference and it does not equal 0
                                #print(currentVtx,t)
                                windowFull,maxPathLength,maxPath,maxPloidyChanges,maxImputes=maxPathLenfindMaxPath1PloidyChange(t,this_Ploidy,ploidchangesSoFar+1,imputesSoFar,pathSoFar+[str(t)],SeenGenomes.copy(),EdgeList)
                                if not windowFull and maxPathLength > longestLength:
                                    longestLength=maxPathLength
                                    longestPath=maxPath
                                    longestPloidyChanges=maxPloidyChanges  
                                    longestPathImputes=maxImputes
                    if len(this_edgeList)==0: #if we didn't find anything, then look 2 windows ahead.
                        ubernextWindow=this_window+2
                        for udx, u in enumerate(Positions):
                            if ubernextWindow==int(u[0]): #if this node is in this window
                                if len(set(G.vs['genomes'][currentVtx]).intersection(set(G.vs['genomes'][udx])))/len(set(G.vs['genomes'][currentVtx]).union(set(G.vs['genomes'][udx])))>=0.8: #if they share 90% of their haplotypes
                                    if G.vs['ploidy'][udx]==G.vs['ploidy'][currentVtx] and G.vs['ploidy'][udx]-Subtractions[udx] !=0: #if nodes are same ploidy
                                        for pdx, p in enumerate(Positions):
                                           if int(p[0])==this_window+1: # if this is an intermediate node
                                               #If this node shares enough haplotypes to be imputed through
                                               if len(set(G.vs['genomes'][currentVtx]).intersection(set(G.vs['genomes'][pdx])))/len(set(G.vs['genomes'][currentVtx]))>0.5: #Not so stringent on the imputed node. It just has to contain more than half of the haplotypes from the node we are leaving
                                                   tmpPath=pathSoFar+["imputed_"+str(pdx)]
                                                   windowFull,maxPathLength,maxPath,maxPloidyChanges,maxImputes=maxPathLenfindMaxPath1PloidyChange(udx,this_Ploidy,ploidchangesSoFar,imputesSoFar+1,tmpPath+[str(udx)],SeenGenomes.copy(),EdgeList)
                                                   if not windowFull and maxPathLength > longestLength:
                                                       longestLength=maxPathLength
                                                       longestPath=maxPath
                                                       longestPloidyChanges=maxPloidyChanges  
                                                       longestPathImputes=maxImputes
    
                    if longestLength >= len(maxmaxPath):
                        maxmaxPath=longestPath
                        maxmaxWindowStart=this_window
                        maxmaxPloidyChanges=longestPloidyChanges
                        maxmaxPathImputes=longestPathImputes
                    lasthits[currentVtx]=longestLength

        lastThreadlenSum=threadLenSum
        threadLenSum+=len(maxmaxPath)
        lastMax=len(maxmaxPath)
        #print(threadLenSum)
        for mdx,m in enumerate(maxmaxPath):
            if "imputed" not in m: #if this is a regular node call
                Subtractions[int(m)]+=1
                WindowsThreadCalls[maxmaxWindowStart+mdx]+=1 #track how many times we've seen a thread through this window
            else: #if this a window that we have imputed through - then increase the thread count for this window.
                imputed_window=int(Positions[int(m.split("_")[1])][0])
                WindowsThreadCalls[imputed_window]+=1
        threadCounter+=1
        if len(maxmaxPath) > 0:
            outer.append([str(threadCounter),str(len(maxmaxPath)),str(maxmaxWindowStart),str(maxmaxPath)])


    #reset ploidy
    #G.vs['ploidy'] = ploidy
    genomes={}
    for v in G.vs:
        genomes[v['name']] = vars(hapoblist[v['name']])['genomes']

    #See how much of WhR haps are visible
    WhR_seen=0
    WhR_visible=0
    for pdx,p in enumerate(Positions):
        if p[0]<=(pericentromere[0]-1)/100000 or p[0] >=(pericentromere[1]-1)/100000: 
            for g in genomes[pdx]:
                if "A_" in g:
                    WhR_seen+=1
                    if G.vs['ploidy'][pdx] != -2: #if this is not a blank node
                        WhR_visible +=1

    print("")
    print("WhR nodes visible",WhR_visible,WhR_seen,WhR_visible/WhR_seen,sep="\t")
    print("WhR nodes visible",WhR_visible,WhR_seen,WhR_visible/WhR_seen,sep="\t",file=outFile_log)
    print("WhR nodes visible",WhR_visible,WhR_seen,WhR_visible/WhR_seen,sep="\t",file=outFile_accuracy)
    #SeenNodes=Counter() # Counter dictionary which stores how many times we've seen a thread traverse a node.
    ThreadedAccurateCount=0
    Threaded_AccurateCount_NonInverted_offset=0
    GlobalAccurateCount=0
    HapSwitchCount=0
    HapSwitchNodeCount=0
    twocount=0
    for o in outer:
        if int(o[1]) >1:
            twocount+=1
    print(len(outer),"threads, of which",twocount,"are length 2+")
    print(len(outer),"threads, of which",twocount,"are length 2+", file=outFile_log)
    ThreadedCount=0
    GlobalCount=0 # global here just means including threads of length 1
    Threaded_Count_NonInverted_offset=0


    freebies=Counter() #keep track of how many times we've seen an inverted WhiteRose Node in a thread window. These can then be used to "discount" a false positive in another thread passing through this window via a non-white-rose node
    pos_freebie_blacklist={}
    for odx, o in enumerate(outer):
        this_nodes=ast.literal_eval(o[3])
        for n in this_nodes:
            WhR_inv=False
            if "imputed" in n:
                n = n.split("_")[1]
            window_pos=Positions[int(n)][0] #get window pos # find the window 
            #print(n, window_pos,gens, inv_list)
            for pdx,p in enumerate(Positions): #
                if p[0]==window_pos:#for each node within this window
                    inv_list=G.vs['inversionShape'][pdx]
                    if "triangle-down" in inv_list: #if this node has an inversion
                        gens=G.vs['genomes'][pdx]
                        for g in gens:
                            if "A_" in g:   #if this inverted node contains a white rose haplotype
                                if g not in pos_freebie_blacklist:
                                    pos_freebie_blacklist[g]=[]
                                if p[0] not in pos_freebie_blacklist[g]:
                                    freebies[Positions[int(n)][0]]+=1
                                    pos_freebie_blacklist[g].append(Positions[int(n)][0])

    print("Windows with inverted White Rose nodes:",freebies)
    print("Windows with inverted White Rose nodes:",freebies,sep="\t",file=outFile_log)

    #getMaxGenome
    for odx,o in enumerate(outer): # for each thread line
        this_nodes=ast.literal_eval(o[3])
        #Calculate the dominant genomes for this thread
        threadWeight=Counter()
        for n in this_nodes:
            if "imputed" in n:
                n = n.split("_")[1]
            this_genomes=genomes[int(n)]
            for t in this_genomes:
                threadWeight[t]+=1
        #print(odx+1,threadWeight)
        maxLen=max(list(threadWeight.values()))
        maxGenomes=[]
        for w in threadWeight:
            if threadWeight[w]==maxLen:
                maxGenomes.append(w)   
        outer[odx].append(str(maxGenomes))

    #reset genomes
    genomes={}
    for v in G.vs:
        genomes[v['name']] = vars(hapoblist[v['name']])['genomes']

    #Sort threads by length to account for bug where longer thread will very occasionally be called later than a shorter thread
    outer=sorted(outer,reverse=True,key=lambda x: int(x[1]))
    ####Now calculate hapSwitches
    for odx,o in enumerate(outer): # for each thread line

        has_hapSwitch=False
        this_nodes=ast.literal_eval(o[3])
        
        #Calculate the dominant genome for this thread
        threadWeight=Counter()
        for n in this_nodes:
            if "imputed" in n:
                n = n.split("_")[1]
            this_genomes=genomes[int(n)]
            for t in this_genomes:
                threadWeight[t]+=1
        #print(odx+1,threadWeight)
        maxLen=0
        maxGenome="NA"
        for w in threadWeight:
            if threadWeight[w]>=maxLen:
                if "A_" in w:
                    maxGenome=w
                    maxLen=threadWeight[w] 
        
        #Now subtract this thread from genomes and score accuracy
        correctNodes=[]
        lastnodehasHapSwitch=False
        thisnodehasHapSwitch=False
        firstInSwitch=False
        for n in this_nodes:
            
            if "imputed" in n:
                n = n.split("_")[1]

            if maxGenome in genomes[int(n)]:
                GlobalAccurateCount+=1
                if int(o[1]) >1: # only consider a thread greater than 1
                    ThreadedAccurateCount+=1
                    if G.vs['inversionShape'][int(n)]=="triangle-down": #if this is an inverted node #THIS IS WRONG.
                        Threaded_AccurateCount_NonInverted_offset-=1
                idx= genomes[int(n)].index(maxGenome)
                genomes[int(n)].pop(idx) # drop that genome out of consideration for later threads
                correctNodes.append(True)
                thisnodehasHapSwitch=False
                firstInSwitch=False
            else:
                correctNodes.append(False)
                if freebies[Positions[int(n)][0]]>0: #if we have a freebie here form an inversion
                    Threaded_Count_NonInverted_offset-=1
                    freebies[Positions[int(n)][0]] -=1
                if "A_" in ",".join(genomes[int(n)]):
                    thisnodehasHapSwitch=True
                    print("hapswitch node at thread",o[0],"node",n)
                    print("hapswitch node at thread",o[0],"node",n,sep="\t",file=outFile_log)
                    if not lastnodehasHapSwitch:
                        firstInSwitch=True
                    if thisnodehasHapSwitch and lastnodehasHapSwitch:
                        HapSwitchNodeCount+=1
                        if firstInSwitch: # if this is the beginning of a bigger hapswitch,.
                            has_hapSwitch=True
                            
            lastnodehasHapSwitch=thisnodehasHapSwitch

            GlobalCount+=1
            if int(o[1]) >1: # only consider a thread greater than 1
                ThreadedCount+=1
                if G.vs['inversionShape'][int(n)]=="triangle-down": #if this is an inverted node #THIS IS WRONG.
                        Threaded_Count_NonInverted_offset-=1
        
        outer[odx].append(maxGenome)
        outer[odx].append(str(correctNodes))
        if has_hapSwitch:
            print("hapSwitch in thread:", odx+1)
            HapSwitchCount+=1
            HapSwitchNodeCount+=1 # count the first node that we would have otherwise ignored



    print("Thread 2+ Accuracy:",ThreadedAccurateCount,ThreadedCount,ThreadedAccurateCount/ThreadedCount)
    print("Thread 2+ Accuracy(non-inverted):",ThreadedAccurateCount+Threaded_AccurateCount_NonInverted_offset ,ThreadedCount+Threaded_Count_NonInverted_offset,(ThreadedAccurateCount+Threaded_AccurateCount_NonInverted_offset)/(ThreadedCount+Threaded_Count_NonInverted_offset))
    print("Thread 1+ Accuracy:",GlobalAccurateCount,GlobalCount,GlobalAccurateCount/GlobalCount)
    outFile_log.write("Thread (2+) Accuracy:"+str(ThreadedAccurateCount)+" "+str(ThreadedCount)+" "+str(ThreadedAccurateCount/ThreadedCount)+"\n")
    print("Thread 2+ Accuracy(non-inverted):",ThreadedAccurateCount+Threaded_AccurateCount_NonInverted_offset ,ThreadedCount+Threaded_Count_NonInverted_offset,(ThreadedAccurateCount+Threaded_AccurateCount_NonInverted_offset)/(ThreadedCount+Threaded_Count_NonInverted_offset),sep="\t",file=outFile_log)

    outFile_log.write("Thread (1+) Accuracy:"+str(GlobalAccurateCount)+" "+str(GlobalCount)+" "+str(GlobalAccurateCount/GlobalCount)+"\n")
    print("HapSwitch-looking errors in",HapSwitchCount,"threads, node-wise:", HapSwitchNodeCount, ThreadedCount,HapSwitchNodeCount/ThreadedCount)
    outFile_log.write("HapSwitch-looking errors in "+str(HapSwitchCount)+" threads, node-wise: "+ str(HapSwitchNodeCount)+" "+str(ThreadedCount)+" "+str(HapSwitchNodeCount/ThreadedCount)+"\n")
    print("Thread 2+ Accuracy",ThreadedAccurateCount,ThreadedCount,ThreadedAccurateCount/ThreadedCount,sep="\t",file=outFile_accuracy)
    print("Thread 2+ Accuracy(non-inverted):",ThreadedAccurateCount+Threaded_AccurateCount_NonInverted_offset ,ThreadedCount+Threaded_Count_NonInverted_offset,(ThreadedAccurateCount+Threaded_AccurateCount_NonInverted_offset)/(ThreadedCount+Threaded_Count_NonInverted_offset),sep="\t",file=outFile_accuracy)
    print("Thread 1+ Accuracy",GlobalAccurateCount,GlobalCount,GlobalAccurateCount/GlobalCount,sep="\t",file=outFile_accuracy)
    print("HapSwitch errors (WhR-WhR)",HapSwitchCount,"affecting", HapSwitchNodeCount, "nodes",sep="\t",file=outFile_accuracy)
    threadCount1=Counter()
    for wdx,w in enumerate(WindowsThreadCalls):
        if (wdx*100000)+1 <= pericentromere[0] or (wdx*100000)+1 >= pericentromere[1]:#if this window is not in the pericentromere
            threadCount1[w]+=1
    print("Number of threads called per window (thread length 1+, non-pericentromeric):",threadCount1)
    print("Number of threads called per window (thread length 1+, non-pericentromeric):",threadCount1,file=outFile_log)
    for c in sorted(threadCount1):
        print(threadCount1[c])
        print(threadCount1[c],file=outFile_log)
        print(c,threadCount1[c],sep="\t",file=outFile_WindowCov)

    totalSeenMb=0
    for odx,o in enumerate(outer):
        totalSeenMb+=float(o[1])/10
        fraction_covered=totalSeenMb/(totalCall/10*4)
        outer[odx].append(str(totalSeenMb))
        outer[odx].append(str(fraction_covered))
    print("Thread\tLength\tStart_Window\tNode_Path\tMostSeenGenomes\tCalledHap\tWhR_hap_Node_Accuracy\tGenome_covered(Mb)\tFraction_of_total")
    for o in outer:
        #print(outer)
        print(int(o[0]),int(o[1]),int(o[2]),",".join(ast.literal_eval(o[3])),",".join(ast.literal_eval(o[4])),o[5],",".join([str(x) for x in ast.literal_eval(o[6])]),float(o[7]),float(o[8]),sep="\t")
        print(int(o[0]),int(o[1]),int(o[2]),",".join(ast.literal_eval(o[3])),",".join(ast.literal_eval(o[4])),o[5],",".join([str(x) for x in ast.literal_eval(o[6])]),float(o[7]),float(o[8]),sep="\t", file=outFile_log)
        print(int(o[0]),int(o[1]),int(o[2]),",".join(ast.literal_eval(o[3])),",".join(ast.literal_eval(o[4])),o[5],",".join([str(x) for x in ast.literal_eval(o[6])]),float(o[7]),float(o[8]),sep="\t", file=outFile_threads)


    #Draw plot with threads
    #Get positions
    Threads=[]
    for o in outer:
        xposes=[]
        yposes=[]
        this_nodes=ast.literal_eval(o[3])
        for n in this_nodes:
            if "imputed" in n:
                n = n.split("_")[1]
            p=Positions[int(n)]
            xposes.append(p[0])
            yposes.append(p[1]+(random.randint(0,10))*0.1)
        Threads.append((xposes,yposes))



    #Get output with reflines
    fig,ax=plt.subplots()
    fig.set_figwidth(60)
    fig.set_figheight(5)
    ax.set_xlim(-2,numWindows+2)
    ax.set_ylim(-2,9*scale)
    #plt.hlines(1,0,230)

    #get pericentromere block
    peri_rect=plt.Rectangle(((pericentromere[0]/100000)-0.5,-1),pericentromere[1]/100000-pericentromere[0]/100000,9*scale)
    coll=[peri_rect]
    pc=PatchCollection(coll,facecolor="green",alpha=0.05,zorder=-2)
    ax.add_collection(pc)
    
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[0],color="green",zorder=-1)
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[1],color="purple",zorder=-1)
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[2],color="firebrick",zorder=-1)
    ax.plot(list(range(0,numWindows,1)),ReferenceLines[3],color="pink",zorder=-1)

    for t in Threads:
        ax.plot(t[0],t[1])

    ig.plot(G, target=ax,layout=chromlayout, vertex_shape=G.vs['inversionShape'], vertex_size=[p+4 for p in G.vs['ploidy']], vertex_color=G.vs['cor'], vertex_frame_color="black", vertex_frame_width=0.1,  edge_width=0.1, edge_arrow_size=0.0, bbox=[3000, 300])
    

    #Final touches and plot
    ax.set_xticks(list(range(0,numWindows+1,1)))
    ax.set_xticklabels(list(range(0,numWindows+1,1)),rotation=90, size=7)
    plt.tight_layout()
    fig.savefig(outPrefix+".gridRefs.ThreadsV7.pdf",dpi=300) # plot with threads
    #plt.show()




# END

# <editor-fold desc="OBSOLETE FUNCTIONS">

#MAIN
plot_threads()