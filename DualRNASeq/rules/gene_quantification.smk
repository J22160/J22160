#######################################   gene_quantification.smk   ######################################################################
#
#    Description:
#       This Snakemake pipeline performs gene quantification using HTSeq-count for both host and pathogen samples. 
#       The script processes aligned BAM files and utilizes genome annotation files (GTF format) to count the number 
#       of genes expressed in the input samples. It generates count tables for downstream differential expression analysis.
#
#   Inputs:
#       - BAM files for host and pathogen:
#       - GTF annotation files for host and pathogen
#           
#   Outputs:
#       - Gene count tables for host and pathogen. These tables contain counts of reads assigned to genes for each sample, formatted for downstream analysis.
#
#   Tools used:
#       - HTSeq-count: A tool for counting reads in features (e.g., genes) defined in GTF annotation files.
#
#   Parameters used for HTSeq-count:
#       - `-f bam`: Input file format is BAM.
#       - `-r pos`: Input reads are sorted by position.
#       - `-s no`: Strand-specific counting is disabled (non-strand-specific library).
#       - `-t gene`: Specifies the feature type to count as "gene."
#       - `-i gene_id`: Uses the "gene_id" attribute from the GTF file to assign reads to genes.
#       - `--with-header`: Includes the header in the output count table for metadata and clarity.
#
#   Usage:
#       This rule quantifies gene expression for the host and pathogen genome using aligned host and pathogen BAM files and the provided GTF annotation.
# 
#   Applications:
#       - The resulting count tables are directly usable for differential gene expression analysis (e.g., DESeq2, edgeR).
#       
##################################################################################################################################################################

HOST_ANNOTATION = config["host_annotation"]
PATHOGEN_ANNOTATION = config["pathogen_annotation"]

# Rule: Quantify genes in the host genome using HTSeq-count
rule host_gene_quantification:
    
    input:
        bam="output/alignment/mapped/host/{sample}_aligned.bam",
        gtf=HOST_ANNOTATION
    
    output:
        count_table="output/gene_count/host/host_gene_count.txt"
    
    params:
        feature_type= config["htseq-count"]["feature_type"]
        gtf_attribute= config["htseq-count"]["gtf_attribute"]
    
    log:
        "logs/gene_quantification/host/host_gene_quantification.log"
    
    conda:
        "environments/htseqcount.yaml"
    
    threads: 16

    benchmark:
        "output/benchmarks/gene_count/host/host.gene.count.benchmark.txt"
    
    resources:
        mem_mb= config["gene_quant_host"]["mem_mb"]
        
    shell:
        """
        mkdir -p output/gene_count/host &&
        htseq-count -f bam -r pos -s no -t {params.feature_type} -i {params.gtf_attribute} --with-header {input.bam} {input.gtf} > {output.count_table} 2> {log}
        """

# Rule: Quantify genes in the pathogen genome using HTSeq-count
rule pathogen_gene_quantification:
    
    input:
        bam="output/alignment/mapped/pathogen/{sample}_aligned.bam",
        gtf=PATHOGEN_ANNOTATION
    
    output:
        count_table="output/gene_count/pathogen/pathogen_gene_count.txt"
    
    params:
        feature_type= config["htseq-count"]["feature_type"]
        gtf_attribute= config["htseq-count"]["gtf_attribute"]
        
    log:
        "logs/gene_quantification/pathogen/pathogen_gene_quantification.log"
    
    threads: 16

    conda:
        "environments/htseqcount.yaml"
    
    benchmark:
        "output/benchmarks/gene_count/pathogen/pathogen.gene.count.benchmark.txt"

    resources:
        mem_mb= config["gene_quant_pathogen"]["mem_mb"]
    
    shell:
        """
        mkdir -p output/gene_count/pathogen &&
        htseq-count -f bam -r pos -s no -t {params.feature_type} -i {params.gtf_attribute} --with-header {input.bam} {input.gtf} > {output.count_table} 2> {log}
        """
