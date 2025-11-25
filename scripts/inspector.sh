#!/bin/bash

sampleid=$1 # HEL_2
hap=$2 #hap1
draft=$3 # analysis/HEL_2/hifiasm-draft_HEL_2/hap/HEL_2.p_utg.fa
raw_reads=$4 # hifi_reads/HEL_2-reads.fastq
workdir=$4

mkdir -p ${workdir}/$sampleid/corrected/${hap}
cp $draft ${workdir}/$sampleid/corrected/${hap}/${sampleid}.${hap}.contig_corrected_round0.fa

iteration=0
(
while [ $iteration -lt 3 ]; do
    echo "running counter round $iteration.."
    inspector.py -c ${workdir}/$sampleid/corrected/${hap}/${sampleid}.${hap}.contig_corrected_round${iteration}.fa -r $raw_reads -o ${workdir}/$sampleid/corrected/${hap}/inspector_$iteration --datatype hifi
    inspector-correct.py -i ${workdir}/$sampleid/corrected/${hap}/inspector_$iteration --datatype pacbio-hifi -o ${workdir}/$sampleid/corrected/${hap} -t 16
    iteration=$((iteration + 1))
    cp ${workdir}/$sampleid/corrected/${hap}/contig_corrected.fa ${workdir}/$sampleid/corrected/${hap}/${sampleid}.${hap}.contig_corrected_round$iteration.fa
    echo "finished that round! on to round $iteration"
done
)
mv ${workdir}/$sampleid/corrected/${hap}/${sampleid}.${hap}.contig_corrected_round2.fa ${workdir}/$sampleid/corrected/${sample}.${hap}.contig_corrected_inspectorFINAL.fa
