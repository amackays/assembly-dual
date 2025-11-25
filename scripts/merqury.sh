#!/bin/bash

sampleid=$1 # HEL_2
draftdir=$2 # analysis/HEL_2/purge_HEL_2

############ merqury
echo "Evaluating assemblies with Merqury..."
echo $(merqury.sh --version)
cd ${draftdir}/stats/
cp -r ../../raw_genomescope/${sampleid}-hifi_reads.meryl ./
cp ../${sampleid}.bp.hap*.p_ctg.fa ./
merqury.sh ${sampleid}-hifi_reads.meryl ${sampleid}.bp.hap1.p_ctg.fa ${sampleid}.bp.hap2.p_ctg.fa merqury
echo "finished merqury! check stats/logs/ for errors if this step finished quickly"
