wd=your/working/directory
##Requires paired-end short-reads .fq.gz of tetraploid cultivar: r1 and r2
r1=foobar.R1.fq.gz
r2=foobar.R2.fq.gz
kmers_to_count_fasta=Revisions_package/WhiteRose_chr06/a1_count_kmers/A.0.15_51mer_chr06_div0.1_counts.nodekmers.dumps.fa
nodekmers=${wd}/WhiteRose_chr06/a0_haplotype_graph/nodekmers_k51_2024_08_07_haplotype_graph_chr06_div0.1_2024_08_02.txt

jellyfish count -m 51 -s 240M -C -o A_51mer_chr06_div0.1_counts.nodekmers.jf --if ${kmers_to_count_fasta} <(zcat ${r1}) <(zcat ${r2}) # wont work without shor read files
#Get histogram of counts
jellyfish histo -f A_51mer_chr06_div0.1_counts.nodekmers.jf> ${wd}/WhiteRose_chr06/a1_count_kmers/A.0.15_51mer_chr06_div0.1_counts.nodekmers.histo.txt #wont work without output froma above
#dump counts - scripts in this package can be run from here
jellyfish dump A_51mer_chr06_div0.1_counts.nodekmers.jf > ${wd}/WhiteRose_chr06/a1_count_kmers/A.0.15_51mer_chr06_div0.1_counts.nodekmers.dumps.fa
#convert dumps.fa to a two-column count file.
python ${wd}/scripts/a0_convert_dumpsFasta_to_tsv.py ${wd}/WhiteRose_chr06/a1_count_kmers/A.0.15_51mer_chr06_div0.1_counts.nodekmers.dumps.fa

#Model probabilities of different counts
python ${wd}/scripts/a1_fitPeaks_v7.py ${wd}/WhiteRose_chr06/a1_count_kmers/A.0.15_51mer_chr06_div0.1_counts.nodekmers.histo.txt 4

####MOVE FOCUS TO EM###################################
cd ${wd}/a2_EM/
kmerCounts=${wd}/a1_countKmers/A_51mer_chr06_div0.1_counts.nodekmers.dumps.tsv
out=${wd}/a2_EM/A_51mer_chr06_div0.1_node_k51_stats.txt
python ${wd}/scripts/a2_get_kmer_stats.py --out ${out} $kmerCounts $nodekmers 
#convert nodekmers to a fasta file
bash ${wd}/scripts/a0_convert_nodekmers_to_fasta.sh ${nodekmers}

#Sort and fill node k stats for EM
python ${wd}/scripts/a2_sort_and_fill_nodekstats.py ${out}

#Run EM
graphFile=${wd}/WhiteRose_chr06/a0_haplotype_graph/haplotype_graph_chr06_div0.1_2024_08_02.txt
nodesStatsFile=${wd}/WhiteRose_chr06/a2_EM/A.0.15_51mer_chr06_div0.1_node_k51_stats.sortedFilled.txt
modelFile=${wd}/WhiteRose_chr06/a1_countKmers/A_51mer_chr06_div0.1_counts.nodekmers.histo.txt.gaussianFit.tsv
mkdir ${wd}/a2_EM/${graphName}/EM_results
outDirEMtsv=${wd}/WhiteRose_chr06/a2_EM/
python ${wd}/scripts/a2_EM.v0.4.py call -G $graphFile -N $nodesStatsFile -M $modelFile --CalcPriors -o ${outDirEMtsv}/EM.v04.WhR.w1_d600k.chr06_div0.1.results.tsv -i 100 -d 600000 -w 1 -K 112

####MOVE FOCUS TO THREADING############################
#Calculate threads
thrd_script=${wd}/scripts/a3_threading_simpleImputeV7_WhiteRose.py
chromarray=(chr06)
peristart=("06473683") # from WhiteRose_chr06/a3_threading/DM_pericentromeric_coords.txt
periend=("34511084") # from WhiteRose_chr06/a3_threading/DM_pericentromeric_coords.txt
#from /netscratch/dep_mercier/grp_schneeberger/projects/pantohap/a1_ChromosomeLevelGraphPlotting/allChroms_div0.1_2024_08_02/DM_pericentromeric_coords.txt)
for i in {0..1}; do
    #if [ $i -gt -1 ] && [ $i -lt 12 ]; then
    echo $i ${chromarray[$i]} ${peristart[$i]} ${periend[$i]}

    chr=${chromarray[$i]}
    div=0.1
    this_pstart=${peristart[$i]}
    this_pend=${periend[$i]}

    graph=${wd}/WhiteRose_chr06/a0_haplotype_graph/haplotype_graph_chr06_div0.1_2024_08_02.txt
    emfin=${wd}/WhiteRose_chr06/a2_EM/EM.v04.WhR.w1_d600k.chr06_div0.1.results.tsv
    inversions=${wd}/WhiteRose_chr06/a3_threading/all_inversions_coordinates.txt
    outprefix=${wd}/WhiteRose_chr06/a3_threading/results
    python $thrd_script $emfin $graph $outprefix ${this_pstart},${this_pend} $chr $inversions
done

#Get summary stats
chrs=chr06
chr_target_Mbs=124.0
inDir=${wd}/WhiteRose_chr06/a3_threading/results
outDir=${wd}/WhiteRose_chr06/a3_threading/stats
python ${wd}/a3_threading_getStats_WhiteRose_v7.py $chrs $chr_target_Mbs $inDir $outDir
