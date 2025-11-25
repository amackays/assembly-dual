#!/bin/bash

sampleid=$1 # HEL_2
draft=$2 # analysis/HEL_2/hifiasm-draft_HEL_2/HEL_2.p_utg.fa
hap=$3 #hap1
raw_reads=$4 # hifi_reads/{sample}-reads.fastq
workdir=$5 # analysis-ec

echo $(date -u) "Aligning hifi reads to primary assembly with minimap2..."

hap=$(echo ${draft} | cut -d "." -f 3)

#Make output directory
mkdir -p ${workdir}/$sampleid/purged-draft/$hap
outdir=${workdir}/$sampleid/purged-draft/$hap

##Align hifi reads to primary assembly with minimap2 with recommended hifiread parameters set
minimap2 -t 16 -x map-hifi $draft $raw_reads | gzip -c - > $outdir/${hap}_asm.paf.gz

echo $(date -u) "Primary alignment complete!"
echo $(date -u) "Getting base coverage and statistics files..."

##bounds
##transition
lower_bound=`awk '{print $2}' ${workdir}/${sampleid}/raw_genomescope/bound_stats.txt`
##upper bound
upper_bound=`awk '{print $3}' ${workdir}/${sampleid}/raw_genomescope/bound_stats.txt`

##Split assembly 
split_fa $draft > $outdir/${hap}_asm.split

##Move into analysis directory for all later steps to make purge_dups happy
cd $outdir

## calculate cutoffs
pbcstat ${hap}_asm.paf.gz #produces PB.base.cov and PB.stat files
calcuts -m $lower_bound -u $upper_bound -d 0 PB.stat > ${hap}_asm_cutoffs

echo $(date -u) "Base coverage and statistics files generated!"
echo $(date -u) "Generating self alignments..."

##align split assembly to itself
minimap2 -xasm5 -DP -k19 -w 19 -m200 ${hap}_asm.split ${hap}_asm.split | gzip -c - > ${hap}_asm.split.self.paf.gz

echo $(date -u) "Self alignments generated!"
echo $(date -u) "Purging haplotigs and overlaps..."

##purge duplicates
purge_dups -2 -a 70 -T ${hap}_asm_cutoffs -c PB.base.cov ${hap}_asm.split.self.paf.gz > ${hap}_asm_dups.bed

## if BUSCO is negatively affected, try adjusting -a parameter from default of 70 up to 75 or 80
## Based on our tests, setting -a above 80 did not improve BUSCO gains, so we chose -a 80

echo $(date -u) "Haplotigs and overlaps purged!"
echo $(date -u) "Retrieving purged primary and haplotig sequences from draft assembly!"

## extract purged assembly from primary assembly
get_seqs -p ${sampleid}.${hap}.p_ctg -e ${hap}_asm_dups.bed ../../hifiasm-draft/${sampleid}.bp.${hap}.p_ctg.fa

echo $(date -u) "Sequences retreieved from draft assembly!"
