
#################################  get_data.smk  #########################################################################

#         This snakemake script is designed to automatically detect paired-end FASTQ files in a specified directory,
#         organize the data by sample IDs and read directions (forward or 1/ reverse or 2), and export the structured data into
#         a YAML file (`samples.yaml`). The output file is used as input for dual-RNASeq workflows,
#         facilitating downstream analysis such as alignment, quantification, and differential expression studies.

#   Use Case:
#    - Streamlines the preprocessing of paired-end sequencing data.
#    - Provides an organized structure for bioinformatics workflows requiring paired-end file paths.

###########################################################################################################################

import os
import yaml


def get_paired_samples(data_dir):
    
    #      Detect paired-end FASTQ files in the specified directory and organize them into a dictionary.
    #
    #      Parameters:
    #         data_dir (str): Path to the directory containing FASTQ files.
    #
    #      Returns:
    #         dict: A dictionary with sample IDs as keys and paths to R1/R2 reads as nested values.
    #          
    #      Example structure:
    #          {
    #              "sample1": {"1": "path/to/sample1_1.fastq.gz", "2": "path/to/sample1_2.fastq.gz"},
    #             "sample2": {"1": "path/to/sample2_1.fastq.gz", "2": "path/to/sample2_2.fastq.gz"}
    #          }
    
    files = [f for f in os.listdir(data_dir) if f.endswith(".fastq.gz")]
    samples = {}

    for file in files:
        # Parse sample ID and read direction (e.g., 1{forward pair} or 2{reverse pair}) from filename
        parts = file.split("_")
        sample_id = "_".join(parts[:-1])
        read_dir = parts[-1].replace(".fastq.gz", "")

        # Organize files by sample ID
        if sample_id not in samples:
            samples[sample_id] = {}
        samples[sample_id][read_dir] = os.path.join(data_dir, file)

    return samples


# Load configuration settings from the YAML configuration file
config_file = "config.yml"
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

# Specify the directory containing raw fastq data
data_dir = config["data_dir"]

# Detect and organize paired-end FASTQ files
samples = get_paired_samples(data_dir)

# Save the organized sample data to an output YAML file in output folder
samples_output = "output/samples.yaml"
with open(samples_output, 'w') as f:
    yaml.dump(samples, f)

# Output detected samples to the console for verification
print(f"Detected input samples: {samples}")
