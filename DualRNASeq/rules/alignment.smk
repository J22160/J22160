########################################   alignment.smk  ################################################################################
#
#     Description:
#        This Snakemake pipeline automates the alignment of dualRNA-Seq data against host and pathogen reference genomes. 
#        The script takes paired-end FASTQ files from a "trimmed data" directory, builds genome index for the host 
#        and pathogen reference genomes (using STAR and Bowtie2, respectively), and performs alignments to generate 
#        sorted BAM files. Unmapped reads from each alignment step are also extracted. 
#        The generated BAM files are essential for downstream gene quantification in RNA-Seq and dual-RNASeq pipelines.
#
#     Inputs:
#        - Paired-end FASTQ files located in "output/trimmed_data/"
#        - Host genome FASTA and annotation files
#        - Pathogen genome FASTA file
#
#     Outputs:
#        - Host genome index created using STAR
#        - Pathogen genome index created using Bowtie2
#        - Aligned BAM files for host and pathogen
#        - Unmapped reads in FASTQ format
#
#     Tools used:
#        - STAR: For host genome alignment
#        - Bowtie2: For pathogen genome alignment
#        - Samtools: For SAM to BAM conversion and sorting
#
#     Usage:
#        This script is designed to be part of a larger RNA-Seq or dual-RNASeq pipeline.

####################################################################################################################################



import os

############# Paths from config.yaml #################################################

# HOST_FASTA: Path to the reference genome FASTA file for the host
HOST_FASTA = config["host_fasta"]

# PATHOGEN_FASTA: Path to the reference genome FASTA file for the pathogen
PATHOGEN_FASTA = config["pathogen_fasta"]

# HOST_ANNOTATION: Path to the GTF annotation file for the host
HOST_ANNOTATION = config["host_annotation"]

# PATHOGEN_ANNOTATION: Path to the GTF annotation file for the pathogen
PATHOGEN_ANNOTATION = config["pathogen_annotation"]

#######################################################################################


# Dynamically determine index directories based on FASTA paths 

# HOST_INDEX_DIR: Directory to store STAR index for the host
HOST_INDEX_DIR = os.path.join(os.path.dirname(HOST_FASTA), "star_index")

# PATHOGEN_INDEX_DIR: Directory to store Bowtie2 index for the pathogen
PATHOGEN_INDEX_DIR = os.path.join(os.path.dirname(PATHOGEN_FASTA), "bowtie2_index")

#########################################################################################



# Rule: Create STAR index for host
# This rule generates a STAR index for the host reference genome using the FASTA and GTF files.
# The index enables fast alignment of sequencing reads to the host genome.

rule build_host_index:
    
    input:
        fasta=HOST_FASTA,
        gtf=HOST_ANNOTATION
    
    output:
        index_dir=directory(HOST_INDEX_DIR)
    
    log:
        "logs/reference_genome_indexing/host_indexing.log"
    
    conda:
        "environments/star.yaml"
    
    threads: 16

    shell:
        """
        mkdir -p {output.index_dir} &&
        STAR --runThreadN 16 \
             --runMode genomeGenerate \
             --genomeDir {output.index_dir} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 100 > {log} 2>&1
        """


# Rule: Create Bowtie2 index for pathogen
# This rule generates a Bowtie2 index for the pathogen reference genome.
# The index speeds up the alignment process of sequencing reads to the pathogen genome.

rule build_pathogen_index:
    
    input:
        fasta=PATHOGEN_FASTA
    
    output:
        index_base=directory(PATHOGEN_INDEX_DIR)
    
    log:
        "logs/reference_genome_indexing/pathogen_indexing.log"
    
    conda:
        "environments/bowtie2.yaml"
    
    shell:
        """
        mkdir -p {PATHOGEN_INDEX_DIR} &&
        bowtie2-build {input.fasta} {output.index_base} > {log} 2>&1
        """


# Rule: Align reads to the host genome using STAR
# This rule aligns paired-end trimmed reads to the host genome using the STAR aligner.
# Outputs include a sorted BAM file for downstream gene quantification and unmapped reads for further analysis.

rule host_alignment:
    
    input:
        r1="output/trimmed_data/{sample}_1_trimmed.fastq.gz",
        r2="output/trimmed_data/{sample}_2_trimmed.fastq.gz",
        index_dir=HOST_INDEX_DIR
    
    output:
        bam="output/alignment/mapped/host/{sample}_aligned.bam",
        unmapped_r1="output/alignment/unmapped/host/{sample}_unmapped1.fastq.gz",
        unmapped_r2="output/alignment/unmapped/host/{sample}_unmapped2.fastq.gz"
    
    log:
        "logs/alignment/host/{sample}_host_alignment.log"
    
    threads: 32

    conda:
        "environments/star.yaml"
    
    benchmark:
        "output/benchmarks/alignment/host/{sample}.host_alignment.benchmark.txt"

    resources:
        mem_mb= config["STAR_alignment_host"]["mem_mb"]
    
    shell:
        """
        mkdir -p output/alignment/unmapped/host &&
        STAR --runThreadN {threads} \
             --genomeDir {input.index_dir} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix output/alignment/mapped/host/{wildcards.sample}_ \
             --outReadsUnmapped Fastx > {log} 2>&1
        """


# Rule: Align reads to the pathogen genome using Bowtie2
# This rule aligns paired-end trimmed reads to the pathogen genome using the Bowtie2 aligner.
# Outputs include a sorted BAM file and unmapped reads for further processing.

rule pathogen_alignment:
    
    input:
        r1="output/trimmed_data/{sample}_1_trimmed.fastq.gz",
        r2="output/trimmed_data/{sample}_2_trimmed.fastq.gz",
        index_dir=PATHOGEN_INDEX_DIR
    
    output:
        bam="output/alignment/mapped/pathogen/{sample}_aligned.bam",
        unmapped_r1="output/alignment/unmapped/pathogen/{sample}_unmapped_1.fastq.gz",
        unmapped_r2="output/alignment/unmapped/pathogen/{sample}_unmapped_2.fastq.gz"
    
    log:
        "logs/alignment/pathogen/{sample}_pathogen_alignment.log"
    
    conda:
        "environments/bowtie2.yaml"
    
    threads:32
    
    benchmark:
        "output/benchmarks/alignment/pathogen/{sample}.pathogen_alignment.benchmark.txt"

    resources:
        mem_mb= config["Bowtie2_alignment_pathogen"]["mem_mb"]

    shell:
        """
        mkdir -p output/alignment/unmapped/pathogen &&
        bowtie2 -x {input.index_base} \
                -1 {input.r1} -2 {input.r2} \
                -S output/alignment/mapped/pathogen/{wildcards.sample}_pathogen.sam \
                --threads {threads} \
                --un-conc-gz output/alignment/unmapped/pathogen/{wildcards.sample}_unmapped%.fastq.gz > {log} 2>&1

        ### Convert SAM to BAM and sort ###
        samtools view -bS output/alignment/mapped/pathogen/{wildcards.sample}_pathogen.sam | \
        samtools sort -o {output.bam}

        ### Rename unmapped files ###
        mv output/alignment/pathogen/unmapped/{wildcards.sample}_unmapped1.fastq.gz {output.unmapped_r1} &&
        mv output/alignment/pathogen/unmapped/{wildcards.sample}_unmapped2.fastq.gz {output.unmapped_r2}

        ### Remove intermediate SAM file ###
        rm output/alignment/mapped/pathogen/{wildcards.sample}_pathogen.sam
        """
