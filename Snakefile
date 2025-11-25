# AMS genome assembly Snakemake pipeline

# this is based on the VGA pipeline, my own experiments with assembly, and https://doi.org/10.1186/s12864-024-10747-8.

# config
import os
import yaml

# functions and setup
with open("config-ec.yaml", "r") as f:
    config = yaml.safe_load(f)

# variable definitions

samples = config["read_loc"]
WORK_DIR = config["workdir"]

#######################################################################################

rule all:
    input:
        expand(WORK_DIR + "/{sample}/decontam/{hap}/cleaned_sequences/{sample}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa", sample=samples.keys(), hap=["hap1", "hap2"]),
        expand(WORK_DIR + "/{sample}/mito/{hap}/{sample}{hap}_blast_out.txt", sample=samples.keys(), hap=["hap1", "hap2"]),
        expand(WORK_DIR + "/{sample}/decontam/{hap}/kraken_sequences.report", sample=samples.keys(), hap=["hap1", "hap2"]),
        expand(WORK_DIR + "/{sample}/full-stat-report.txt", sample=samples.keys())
        #stats_table = expand("path/to/QC/stats.txt", sample=config["sample_id"])

rule genomescope:
    message: "take {input.hifi_zip} and convert to fasta format, then run GenomeScope to determine heterozygosity etc."
    input:
        hifi_zip=lambda wildcards: samples[wildcards.sample],
    output:
        genomescope=WORK_DIR + "/{sample}/raw_genomescope/linear_plot.png",
        fastq="hifi_reads/{sample}-reads.fastq",
        stats=WORK_DIR + "/{sample}/raw_genomescope/SIMULATED_testing.tsv"
    log: "logs/{sample}.genomescope.log"
    params:
        dir=WORK_DIR
    conda:
        "envs/assembly-2.yaml"
    shell:
        """
        mkdir -p hifi_reads
        mkdir -p {params.dir}
        scripts/pre-qc.sh {input.hifi_zip} {wildcards.sample} {params.dir}
        """

rule purge_params:
    message: "using genomescope2 parameters from {input.stats}, calculate bounds for aggressive purging (if required)."
    input:
        stats=WORK_DIR + "/{sample}/raw_genomescope/SIMULATED_testing.tsv"
    output:
        bounds=WORK_DIR + "/{sample}/raw_genomescope/bound_stats.txt"
    log: "logs/{sample}.purge_params.log"
    shell:
        """
        scripts/purge_dup_stat.sh {input.stats} {output.bounds}
        """

rule draft_hifiasm:
    message: "build the first draft assembly from {input.hifi_fastq} with light purging and index the genome."
    input:
        hifi_fastq="hifi_reads/{sample}-reads.fastq",
        stats=WORK_DIR + "/{sample}/raw_genomescope/bound_stats.txt",
    output:
        hap1=WORK_DIR + "/{sample}/hifiasm-draft/{sample}.bp.hap1.p_ctg.fa",
        hap2=WORK_DIR + "/{sample}/hifiasm-draft/{sample}.bp.hap2.p_ctg.fa",
        dir=directory(WORK_DIR + "/{sample}/hifiasm-draft/")
    log: "logs/{sample}.draft_hifiasm.log"
    conda:
        "envs/assembly-2.yaml"
    params:
        assembly=config["assembly_type"]
    shell:
        """
        mkdir -p {output.dir}
        scripts/hifiasm_assemble.sh {input.stats} {input.hifi_fastq} {params.assembly} {output.dir} {wildcards.sample}
        """

rule draft_qc:
    message: "first round of assembly statistics on {input.hap1} and {input.hap2}."
    input:
        hap1=WORK_DIR + "/{sample}/hifiasm-draft/{sample}.bp.hap1.p_ctg.fa",
        hap2=WORK_DIR + "/{sample}/hifiasm-draft/{sample}.bp.hap2.p_ctg.fa",
        hifi_fastq="hifi_reads/{sample}-reads.fastq"
    output:
        genometools=WORK_DIR + "/{sample}/hifiasm-draft/stats/{sample}.seqstat",
        merqury=WORK_DIR + "/{sample}/hifiasm-draft/stats/merqury.qv",
        busco=WORK_DIR + "/{sample}/hifiasm-draft/stats/busco.summary",
        allstats=WORK_DIR + "/{sample}/hifiasm-draft/stats/aggregate.txt"   
    params:
        draft_dir=WORK_DIR + "/{sample}/hifiasm-draft",
    log: "logs/{sample}.draft_qc.log"
    conda:
        "envs/genome-qc.yaml"
    shell:
        """
        meryl --version
        meryl k=31 count {input.hifi_fastq} output {WORK_DIR}/{wildcards.sample}/raw_genomescope/{wildcards.sample}-hifi_reads.meryl
        meryl histogram {WORK_DIR}/{wildcards.sample}/raw_genomescope/{wildcards.sample}-hifi_reads.meryl > {log}

        scripts/genome-qc.sh {wildcards.sample} {input.hap1} {input.hap2} {params.draft_dir} {input.hifi_fastq} {WORK_DIR}
        scripts/merqury.sh {wildcards.sample} {params.draft_dir}
        cat {output.genometools} {output.merqury} {output.busco} > {output.allstats}
        """

rule purge_dups:
    message: "use minimap2 and purge_dups to purge the primary assembly."
    input:
        draft=WORK_DIR + "/{sample}/hifiasm-draft/{sample}.bp.{hap}.p_ctg.fa",
        hifi_fastq="hifi_reads/{sample}-reads.fastq",
        allstats=WORK_DIR + "/{sample}/hifiasm-draft/stats/aggregate.txt"
    output:
        phap=WORK_DIR + "/{sample}/purged-draft/{hap}/{sample}.{hap}.p_ctg.purged.fa"
    log: "logs/{sample}.{hap}.purge_dups.log"
    conda:
        "envs/purgedups.yml"
    shell:
        """
        scripts/purgedups.sh {wildcards.sample} {input.draft} {wildcards.hap} {input.hifi_fastq} {WORK_DIR}
        """

rule second_qc:
    message: "second round of assembly statistics on the purged assembly."
    input:
        hap1=WORK_DIR + "/{sample}/purged-draft/hap1/{sample}.hap1.p_ctg.purged.fa",
        hap2=WORK_DIR + "/{sample}/purged-draft/hap2/{sample}.hap2.p_ctg.purged.fa",
    output:
        genometools=WORK_DIR + "/{sample}/purged-draft/stats/{sample}.seqstat",
        merqury=WORK_DIR + "/{sample}/purged-draft/stats/merqury.qv",
        busco=WORK_DIR + "/{sample}/purged-draft/stats/busco.summary",
        allstats=WORK_DIR + "/{sample}/purged-draft/stats/aggregate.txt"
    params:
        purged_dir=WORK_DIR + "/{sample}/purged-draft",
        hifi_fastq="hifi_reads/{sample}-reads.fastq"
    log: "logs/{sample}.second_qc.log"
    conda:
        "envs/genome-qc.yaml"
    shell:
        """
        scripts/genome-qc.sh {wildcards.sample} {input.hap1} {input.hap2} {params.purged_dir} {params.hifi_fastq} {WORK_DIR}
        cp {input.hap1} {params.purged_dir}/{wildcards.sample}.bp.hap1.p_ctg.fa
        cp {input.hap2} {params.purged_dir}/{wildcards.sample}.bp.hap2.p_ctg.fa
        scripts/merqury.sh {wildcards.sample} {params.purged_dir}
        cat {output.genometools} {output.merqury} {output.busco} > {output.allstats}
        """

rule inspector:
    message: "check purged assembly for errors using raw reads using Inspector."
    input:
        draft=WORK_DIR + "/{sample}/purged-draft/{hap}/{sample}.{hap}.p_ctg.purged.fa",
        reads="hifi_reads/{sample}-reads.fastq",
        allstats=WORK_DIR + "/{sample}/purged-draft/stats/aggregate.txt"
    output:
        inspector=WORK_DIR + "/{sample}/corrected/{hap}/{sample}.{hap}.contig_corrected_inspectorFINAL.fa"
    conda:
        "envs/corrections.yaml"
    log: "logs/inspector.{sample}.{hap}.log"
    shell:
        """
        scripts/inspector.sh {wildcards.sample} {wildcards.hap} {input.draft} {input.reads} {WORK_DIR}
        """

rule ntlink:
    message: "fill gaps after Inspector correction with raw reads using ntLink. Note that this step will fail if the pythonpath is not set to the current conda env!"
    input:
        draft=WORK_DIR + "/{sample}/corrected/{hap}/{sample}.{hap}.contig_corrected_inspectorFINAL.fa",
        reads="hifi_reads/{sample}-reads.fastq"
    output:
        ntlink=WORK_DIR + "/{sample}/corrected/{hap}/{sample}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa"
    log: "logs/ntlink.{sample}-{hap}.log"
    conda:
        "envs/corrections.yaml"
    params:
        dir=WORK_DIR + "/{sample}/corrected/{hap}/"
    shell:
        """
        cp {input.reads} {params.dir}{wildcards.sample}-reads.fastq
        cd {params.dir}
        ntLink_rounds run_rounds_gaps clean rounds=3 target={wildcards.sample}.{wildcards.hap}.contig_corrected_inspectorFINAL.fa reads={wildcards.sample}-reads.fastq
        """

rule fcs_toolkit:
    message: "screen the draft genome for adapter contamination with ncbi fcs-adaptor. NOTE that this rule can be modified to also use ncbi-gx to screen for foreign contaminants, although this should be taken care of with Kraken2 in the previous step."
    input:
        fasta=WORK_DIR + "/{sample}/corrected/{hap}/{sample}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
    params:
        path_to_contam_dir=WORK_DIR + "/{sample}/decontam/{hap}",
        adaptor_sif="scripts/fcs-adaptor.sif"
    output:
        cleaned_seq=WORK_DIR + "/{sample}/decontam/{hap}/cleaned_sequences/{sample}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa"
    conda:
        "envs/contam.yaml"
    log: "logs/fcs_toolkit.{sample}.{hap}.log"
    shell:
        """
        cp {input.fasta} {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/{wildcards.sample}.{wildcards.hap}.ntLink.fa
        scripts/run_fcsadaptor.sh --fasta-input {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/{wildcards.sample}.{wildcards.hap}.ntLink.fa --output-dir {params.path_to_contam_dir} --euk --container-engine singularity --image {params.adaptor_sif}
        """

rule mito:
    message: "search for the mitochondrial genome in the {input.draft} assembly using Heliconius clysonymus as the closest relative complete mitochondrial genome, using MitoFinder and MitFi."
    input:
        draft=WORK_DIR + "/{sample}/corrected/{hap}/{sample}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
        mito_genbank=config["mito_genbank"]
    output:
        mito_assembled=directory(WORK_DIR + "/{sample}/mito/{hap}/{sample}{hap}_MitoFinder_mitfi_Final_Results/"),
        mito_blast=WORK_DIR + "/{sample}/mito/{hap}/{sample}{hap}_blast_out.txt"
    params:
        singularity_mitofinder="scripts/remiallio_default_mitofinder.sif"
    log: "logs/mito.{sample}.{hap}.log"
    shell:
        """
        mkdir -p {WORK_DIR}/{wildcards.sample}/mito/{wildcards.hap}
        singularity run {params.singularity_mitofinder} -j {wildcards.sample}{wildcards.hap} -a {input.draft} -r {input.mito_genbank} -o 5 -p 8 -m 8 --rename-contig yes 
        cp -r {wildcards.sample}{wildcards.hap}/. {WORK_DIR}/{wildcards.sample}/mito/{wildcards.hap}
        """

rule kraken:
    message: "check corrected and gap-filled assembly for contaminants using Kraken2."
    input:
        draft=WORK_DIR + "/{sample}/corrected/{hap}/{sample}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
        reads="hifi_reads/{sample}-reads.fastq"
    output:
        kraken=WORK_DIR + "/{sample}/decontam/{hap}/kraken_sequences.report"
    conda:
        "envs/contam.yaml"
    log: "logs/kraken.{sample}.{hap}.log"
    shell:
        """
        mkdir -p {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}
        tr a-z A-Z < {input.draft} | dustmasker -level 40 -out {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/tmp.fa -outfmt 'fasta'
        tr [:lower:] 'N' < {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/tmp.fa > {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/out.masked
        
        kraken2 --db scripts/k2_core_nt --threads 8 --conf 0.30 --report {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/kraken_sequences.report --use-names \
        --output {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/kraken_sequences.kraken {WORK_DIR}/{wildcards.sample}/decontam/{wildcards.hap}/out.masked
        """

rule third_qc:
    message: "second round of assembly statistics on the draft assembly after correction and decontamination."
    input:
        hap1=WORK_DIR + "/{sample}/corrected/hap1/{sample}.hap1.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
        hap2=WORK_DIR + "/{sample}/corrected/hap2/{sample}.hap2.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
        #hap1=WORK_DIR + "/{sample}/decontam/hap1/cleaned_sequences/{sample}.hap1.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
        #hap2=WORK_DIR + "/{sample}/decontam/hap2/cleaned_sequences/{sample}.hap2.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa",
        hifi_fastq="hifi_reads/{sample}-reads.fastq"
    output:
        genometools=WORK_DIR + "/{sample}/decontam/stats/{sample}.seqstat",
        merqury=WORK_DIR + "/{sample}/decontam/stats/merqury.qv",
        busco=WORK_DIR + "/{sample}/decontam/stats/busco.summary",
        allstats=WORK_DIR + "/{sample}/decontam/stats/aggregate.txt"
    params:
        decontam_dir=WORK_DIR + "/{sample}/decontam/"
    log: "logs/third_qc.{sample}.log"
    conda:
        "envs/genome-qc.yaml"
    shell:
        """
        scripts/genome-qc.sh {wildcards.sample} {input.hap1} {input.hap2} {params.decontam_dir} {input.hifi_fastq} {WORK_DIR}
        cp {input.hap1} {params.decontam_dir}/{wildcards.sample}.bp.hap1.p_ctg.fa
        cp {input.hap2} {params.decontam_dir}/{wildcards.sample}.bp.hap2.p_ctg.fa
        scripts/merqury.sh {wildcards.sample} {params.decontam_dir}
        cat {output.genometools} {output.merqury} {output.busco} > {output.allstats}
        """

rule compile_stats:
    message: "Aggregate genome assembly statistics on raw reads and on all three QC checkpoints in the assembly process - draft, after purging, and after decontamination/correction"
    input:
        genomescope=WORK_DIR + "/{sample}/raw_genomescope/summary.txt",
        draft=WORK_DIR + "/{sample}/hifiasm-draft/stats/aggregate.txt",
        purged=WORK_DIR + "/{sample}/purged-draft/stats/aggregate.txt",
        decontam=WORK_DIR + "/{sample}/decontam/stats/aggregate.txt"
    output:
        stats_report=WORK_DIR + "/{sample}/full-stat-report.txt"
    shell:
        """
        cat {input.genomescope} {input.draft} {input.purged} {input.decontam} > {output.stats_report}
        """

rule quast:
    message: "run QUAST on all three checkpoint assemblies - draft, after purging, after decontamination/correction"
    input:
        hifiasm=expand(WORK_DIR + "/{{sample}}/hifiasm-draft/{{sample}}.bp.{hap}.p_ctg.fa", hap=["hap1", "hap2"]),
        purged=expand(WORK_DIR + "/{{sample}}/purged-draft/{hap}/{{sample}}.{hap}.p_ctg.purged.fa", hap=["hap1", "hap2"]),
        ntlink=expand(WORK_DIR + "/{{sample}}/corrected/{hap}/{{sample}}.{hap}.contig_corrected_inspectorFINAL.fa.k32.w100.z1000.ntLink.gap_fill.3rounds.fa", hap=["hap1", "hap2"]),
    params:
        out_path=WORK_DIR + "/{sample}/all_QUAST"
    threads: 8
    output:
        quast=WORK_DIR + "/{sample}/all_QUAST/report.txt"
    conda:
        "envs/quast.yaml"
    shell:
        """
        mkdir -p {params.out_path}
        quast.py -t 8 --eukaryote --large --k-mer-stats --k-mer-size 31 --est-ref-size 425000000 -o {params.out_path} {input.hifiasm}
        quast.py -t 8 --eukaryote --large --k-mer-stats --k-mer-size 31 --est-ref-size 425000000 -o {params.out_path} {input.purged}
        quast.py -t 8 --eukaryote --large --k-mer-stats --k-mer-size 31 --est-ref-size 425000000 -o {params.out_path} {input.ntlink}
        """
#repeatmasker, workdir, mask

#z/w chromosome identification

# print done

################ REFERENCE GENOMES ONLY: #############

#minimap 2 for old assembly to new

#compare to ihelsar/nucmer/R

#genomescope of old genome?

# print done
