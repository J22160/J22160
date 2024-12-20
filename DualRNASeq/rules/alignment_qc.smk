###########################################   alignment_qc.smk   ############################################################################### 
#
#     Description:
#        This Snakemake rule script performs quality control (QC) checks on aligned BAM files for both host and pathogen genomes. 
#        The script takes BAM files generated from alignment steps and produces detailed QC reports using the QualiMap tool. 
#        These reports are essential for assessing the quality of the alignment and ensuring the data is suitable for downstream 
#        analysis, such as gene quantification or differential expression studies.
#
#     Inputs:
#        - Host aligned BAM files: "output/alignment/mapped/host/{sample}_aligned.bam"
#        - Pathogen aligned BAM files: "output/alignment/mapped/pathogen/{sample}_aligned.bam"
#
#     Outputs:
#        - QC reports in PDF and HTML formats
#
#     Tools used:
#        - QualiMap: Provides comprehensive QC metrics for BAM files.
#
#
#     Parameters provided by QualiMap:
#        - General statistics (e.g., total reads, mapped reads, mapping quality)
#        - Coverage metrics (e.g., mean coverage, coverage uniformity)
#        - GC content analysis
#        - Read distribution along the genome
#        - Alignment quality metrics (e.g., read duplication, mismatch rates)
#
#     Usage:
#        This script is designed to be part of a larger RNA-Seq or dual-RNASeq pipeline, providing essential QC information for BAM files.
#
####################################################################################################################################################


rule qualimap_qc_host:
    
    input:
        bam="output/alignment/mapped/host/{sample}_aligned.bam",
    
    output:
        qc_data_dir= "output/alignment_qc/host/{sample}"
    
    log:
        "log/qualimap/host/{sample}_qualimap.log"
    
    conda:
        "environments/qualimap.yaml"
    
    benchmark:
        "output/benchmarks/alignment_qc/host/{sample}.host.alignmentqc.benchmark.txt"

    shell:
        """
        mkdir -p output/alignment_qc/host/{wildcards.sample} &&
        qualimap bamqc \
            -bam {input.bam} \
            -outdir {output.qc_data_dir} \
            -outformat PDF:HTML > {log} 2>&1
        """

rule qualimap_qc_pathogen:
    
    input:
        bam= "output/alignment/mapped/pathogen/{sample}_aligned.bam",
    
    output:
        qc_data_dir= "output/alignment_qc/pathogen/{sample}"
    
    log:
        "log/qualimap/pathogen/{sample}_qualimap.log"
    
    conda:
        "environments/qualimap.yaml"
    
    benchmark:
        "output/benchmarks/alignment_qc/pathogen/{sample}.pathogen.alignmentqc.benchmark.txt"
    
    shell:
        """
        mkdir -p output/alignment_qc/pathogen/{wildcards.sample} &&
        qualimap bamqc \
            -bam {input.bam} \
            -outdir {output.qc_data_dir} \
            -outformat PDF:HTML > {log} 2>&1
        """

 
