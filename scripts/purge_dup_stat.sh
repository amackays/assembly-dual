#!/bin/bash
statfile=$1
outfile=$2

echo "Extracting parameters for purge_dups..."
avg_cov=`awk 'FNR == 1 {print $3}' $statfile`
lower_bound=$(python -c "import sys; print(round(float(sys.argv[1])*1.5))" "$avg_cov")
upper_bound=$((lower_bound * 3))
printf "${avg_cov}\t${lower_bound}\t${upper_bound}\n" > $outfile
echo "completed stats calculation"
