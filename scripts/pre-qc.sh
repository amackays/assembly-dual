#!/bin/bash
reads=$1
sample=$2
dir=$3
mkdir -p $dir/$sample/raw_genomescope/
mkdir -p $dir/$sample/raw_genomescope/jellyfish/

#convert to fasta
echo "converting .bam file from Revio to .fasta..."
gunzip -c $reads > $sample.tmp.bam
samtools sort -n $sample.tmp.bam -o $sample.tmp.sorted && samtools fastq $sample.tmp.sorted > $sample.tmp.fastq
echo "conversion done"

#jellyfish
echo "Running jellyfish kmer counter for genomescope..."
jellyfish count -C -m 31 -s 100M $sample.tmp.fastq -o $dir/$sample/raw_genomescope/jellyfish/$sample.jf
jellyfish histo $dir/$sample/raw_genomescope/jellyfish/$sample.jf > $dir/$sample/raw_genomescope/jellyfish/$sample.histo
echo "jellyfish done!"

#genomescope
echo "running GenomeScope2 to plot kmer counts and determine purge_dups params..."
genomescope2 -i $dir/$sample/raw_genomescope/jellyfish/$sample.histo -o $dir/$sample/raw_genomescope/ -k 31 --testing
echo "finished GenomeScope2! Cleaning up..."

#cleanup
mv $sample.tmp.fastq hifi_reads/${sample}-reads.fastq
rm $sample.tmp*
echo "Congrats on finishing pre-qc.sh!"

