#!/bin/bash

sampleid=$1 # HEL_2
draft=$2 # analysis/HEL_2/hifiasm-draft_HEL_2/HEL_2.p_utg.fa
hap2=$3 # alt haplotype
draftdir=$4 # analysis/HEL_2/purge_HEL_2
raw_reads=$5 #hifi_reads/HEL_2-reads.fastq
topdir=$6

mkdir -p ${draftdir}/stats/
genome_size=$(grep -e "Genome Haploid Length" ${topdir}/${sampleid}/raw_genomescope/summary.txt | sed 's/,//g ; s/bp//g' | awk '{print $5}')

############ genometools-genometools
echo "Screening ${sampleid} assemblies with genometools"
echo $(gt --version)
gt seqstat -contigs yes -genome ${genome_size} ${draft} ${hap2} > ${draftdir}/stats/${sampleid}.seqstat
echo "finished processing ${sampleid}!"
############ busco
echo "Assessing gene completeness for ${sampleid} using lepidopteran BUSCO for primary/alternate or unitig assembly..."
singularity exec docker://ezlabgva/busco:v5.4.7_cv1 busco -i ${draft} -o ${sampleid} -f --cpu 20 -l lepidoptera_odb10 --out_path ${draftdir}/stats/hap1 -m geno
singularity exec docker://ezlabgva/busco:v5.4.7_cv1 busco -i ${hap2} -o ${sampleid} -f --cpu 20 -l lepidoptera_odb10 --out_path ${draftdir}/stats/hap2 -m geno
echo $(date -u) "Lepidopteran BUSCOs calulated!"
cat ${draftdir}/stats/hap1/${sampleid}/short_summary.specific.lepidoptera_odb10.${sampleid}.txt ${draftdir}/stats/hap2/${sampleid}/short_summary.specific.lepidoptera_odb10.${sampleid}.txt > ${draftdir}/stats/busco.summary
echo "Wrapping up genome-qc.sh script"
