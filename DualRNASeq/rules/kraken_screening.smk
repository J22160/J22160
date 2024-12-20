############################################  kraken_screening.smk  #######################################################################
#
#   Description:
#       This Snakemake rule performs taxonomic classification of unmapped reads from the host genome using Kraken2. 
#       The script takes unmapped FASTQ files as input and screens them against a Kraken database to identify and classify 
#       potential pathogen reads. The output is a comprehensive Kraken report, which details the taxonomic composition 
#       of the reads, enabling further analysis of microbial diversity and pathogen presence.
#
#   Input:
#       - Host unmapped paired-end reads:
#       
#   Output:
#       - Kraken report: "output/bacterial_screening/kraken/{sample}_kraken_report.txt"
#
#   Tool used:
#       - Kraken2: A highly efficient tool for assigning taxonomic labels to DNA sequences using exact k-mer matches.
#
#   Parameters provided by Kraken2 output:
#       - Taxonomic classification of reads with detailed lineage information (e.g., domain, phylum, genus, species)
#       - Total number of classified and unclassified reads
#       - Percentage of reads assigned to each taxonomic group
#       - Read counts at each taxonomic level
#       - Comprehensive report format suitable for visualization and downstream analysis
#
#   Usage:
#       This script is part of a dual-RNASeq pipeline and is used to detect pathogen reads among unmapped host reads, 
#       aiding in the identification of microbial or bacterial components in mixed samples.
#
#####################################################################################################################################


rule kraken_screening:
    
    input:
        r1="output/alignment/unmapped/host/{sample}_unmapped1.fastq.gz",
        r2= "output/alignment/unmapped/host/{sample}_unmapped2.fastq.gz",
    
    output:
        kraken_output="output/bacterial_screening/kraken/{sample}_kraken_report.txt"        
    
    params:
        db=config["kraken2"]["database"]

    log:
        "logs/kraken/{sample}_kraken_screening.log"
    
    conda:
        "environments/kraken2.yaml"
    
    threads: 16

    benchmark:
        "output/benchmarks/kraken_screening/{sample}.kraken.benchmark.txt"

    shell:
        """
        mkdir -p output/bacterial_screening/kraken
        kraken2 --db {params.db} \
            --threads {threads}} \
            --gzip-compressed \
            --use-names \
            --paired {input.r1} {input.r2} \
            --output {output.kraken_output} \
            --report {output.kraken_report} > {log} 2>&1
        """
