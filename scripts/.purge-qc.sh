#!/bin/bash

sampleid=$1 # HEL_2
draft=$2 # analysis/HEL_2/purge-draft_HEL_2/HEL_2.p_utg.purge.fa
draftdir=$3 # analysis/HEL_2/purge-draft_HEL_2/
raw_reads=$4

mkdir -p ${draftdir}/stats/
genome_size=$(grep -e "Genome Haploid Length" analysis/${sampleid}/raw_genomescope_${sampleid}/summary.txt | sed 's/,//g ; s/bp//g' | awk '{print $5}')

############ gfastats
echo "Screening assemblies with gfastats..."
echo "processing ${sampleid}..."
echo $(gfastats --version)
gfastats -j 4 --stats -s c ${draft} ${genome_size} > ${draftdir}/stats/${sampleid}.gfastats
echo "finished processing ${sampleid} with gfastats!"
############ genometools-genometools
echo "Screening ${sampleid} assemblies with genometools"
echo $(gt --version)
gt seqstat -contigs yes -genome ${genome_size} ${draft} > ${draftdir}/stats/${sampleid}.seqstat
echo "finished processing ${sampleid}!"
############ busco
echo "Assessing gene completeness for ${sampleid} using lepidopteran BUSCO for primary/alternate or unitig assembly..."
echo $(busco --version)
singularity exec docker://ezlabgva/busco:v5.4.7_cv1 busco -i ${draft} -o ${sampleid} -f --cpu 20 -l lepidoptera_odb10 --out_path ${draftdir}/stats -m geno
echo $(date -u) "Lepidopteran BUSCOs calulated!"
############ merqury
echo "Evaluating assemblies with Merqury..."
echo $(merqury.sh --version)
cd ${draftdir}/stats/
merqury.sh ../../raw_genomescope_${sampleid}/hifi_reads.meryl ${draft} merqury
