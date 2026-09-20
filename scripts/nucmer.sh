#PRE SCRIPT:
#conda activate mummer
# make nucmer_circos dir
# cp ref genome and draft genome in decontam folder to nucmer_circos
# move into that dir (?)

#USAGE:
# ./nucmer.sh ../GCA_917862395.2_iHelSar1.2_genomic.fna analysis/demophoon/decontam/demophoon.p_utg.fa demophoon iHelSar1.2
###############################################################

ref_genome=$1 #ref.fa
draft_genome=$2 #draft.fa
sample=$3 #HEL_2
refname=$4 #Hed2.5

#mkdir -p ${dir}/${sample}/nucmer_circos/
#cp ${ref_genome} ${draft_genome} ${dir}/${sample}/nucmer_circos/
#cd ${dir}/${sample}/nucmer_circos/

# align genomes

nucmer --version

#nucmer -l 100 --prefix demophoontoihelsar GCA_917862395.2_iHelSar1.2_genomic.fna analysis/demophoon/decontam/demophoon.p_utg.fa

nucmer -l 100 --prefix ${sample}to${refname} ${ref_genome} ${draft_genome}

delta-filter -l 1000 -q -r ${sample}to${refname}.delta > ${sample}to${refname}.filt.delta
show-coords -c -l -L 1000 -r -T ${sample}to${refname}.filt.delta > ${sample}to${refname}.filt.coords

#Rscript mummerCoordsDotPlotly.R -i ${sample}to${refname}.filt.coords -o ${sample}to${refname} -s -t -m 500 -q 10000 -l -x
