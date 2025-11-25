#!/bin/bash

#load in variables
stats=$1
reads=$2
assemble_type=$3
path=$4
sample=$5
bound=`awk '{print $3}' $stats`

# assembly with disabled bloom filter, light purging based on genomescope

echo "beginning genome assembly...printing tool version" 

hifiasm --version
hifiasm -t 100 -f0 $assemble_type -l 1 --purge-max $bound -o $path/$sample $reads 2> logs/$sample-assemble.log

echo "wrapped up hifiasm, continuing to generate gfas and index file"

# gfa generation
awk '/^S/{print ">"$2;print $3}' $path/$sample.bp.hap1.p_ctg.gfa > $path/$sample.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' $path/$sample.bp.hap2.p_ctg.gfa > $path/$sample.bp.hap2.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' $path/$sample.a_ctg.gfa > $path/$sample.a_ctg.fa

# index
samtools faidx $path/$sample.bp.hap1.p_ctg.fa
samtools faidx $path/$sample.bp.hap2.p_ctg.fa
#samtools faidx $path/$sample.a_ctg.gfa  > $idx-a

echo "completed gfa generation and indexing!"
