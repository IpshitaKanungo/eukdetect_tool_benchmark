configfile: "config.yaml"

import re
import os
from datetime import datetime

# stamp can come from CLI: --config stamp=20251029-1430
stamp = config.get("stamp") or datetime.now().strftime("%Y%m%d-%H%M%S")
result = config["path"]["result_path"].replace("{stamp}", stamp)


# Raw data
raw_data_dir = config["raw_data_dir"]

file_extensions = tuple(config["file_extensions"])

# Create a list of unique sample basenames
sample_basename = list(
    set(
        os.path.splitext(f)[0].rsplit("_", 1)[0]
        for f in os.listdir(raw_data_dir)
        if f.endswith(file_extensions)
    )
)
print(sample_basename)
def get_path(sample):
    # Strip suffixes and extensions from the sample to match `sample_basename`
    base_sample = sample.rsplit("_", 1)[0]
   # print(base_sample)
    # Ensure the sample exists in `sample_basename`
    if base_sample not in sample_basename:
        raise ValueError(f"Sample '{base_sample}' not found in sample_basename")
    
    # Dynamically build paths
    return [
        os.path.join(raw_data_dir, f"{base_sample}_{pair}.fq")
        for pair in ["1", "2"]
    ]

use_bowtie2 = config["rules"]["bowtie2"]
use_eukrep  = config["rules"]["eukrep"]
use_tiara   = config["rules"]["tiara"]
use_kraken2 = config["rules"]["kraken2"]

bowtie2_result = config["folder"]["bowtie2"]
assembly_result = config["folder"]["assembly"]
megahit_result = config["folder"]["megahit"]
binning_result = config["folder"]["binning"]
refinement_result = config["folder"]["refinement"]
eukrep_result = config["folder"]["eukrep"]
tiara_result = config["folder"]["tiara"]
kraken2_result = config["folder"]["kraken2"]

MEGAHIT_BRANCHES = ["bowtie2"] if use_bowtie2 else ["host"]
BRANCHES = (["bowtie2"] if use_bowtie2 else
    ([b for b, flag in [("eukrep", use_eukrep), ("tiara", use_tiara), ("kraken2", use_kraken2)] if flag] or ["host"]))

rule all:
    input:        
        # Conditionally include Bowtie2 logs if enabled in config
        (expand(f"{result}/{bowtie2_result}/logs/{{sample}}_bowtie2.log", sample=sample_basename) if use_bowtie2 else []),
        expand(f"{result}/{megahit_result}/logs/{{sample}}_{{branch}}_megahit.log", sample=sample_basename, branch=MEGAHIT_BRANCHES),
        (expand(f"{result}/{eukrep_result}/{{sample}}/{{sample}}_eukaryotic.fa", sample=sample_basename) if use_eukrep else []),
        (expand(f"{result}/{kraken2_result}/{{sample}}/{{sample}}_taxonomic_report.tsv", sample=sample_basename) if use_kraken2 else []),
        (expand(f"{result}/{tiara_result}/{{sample}}/{{sample}}_tiara.txt", sample=sample_basename) if use_tiara else []),
        expand(f"{result}/{binning_result}/logs/{{sample}}_{{br}}_binning.log", sample=sample_basename, br=BRANCHES),
        expand(f"{result}/{refinement_result}/logs/{{sample}}_{{br}}_refinement.log", sample=sample_basename, br=BRANCHES),


if use_bowtie2:
    include: "rules/bowtie2.smk"
include: "rules/megahit.smk"
if not use_bowtie2:
    if use_eukrep:
        include: "rules/eukrep.smk"
    if use_tiara:  
        include: "rules/tiara.smk"
    if use_kraken2: 
        include: "rules/kraken2.smk"
include: "rules/binning.smk"
include: "rules/refinement.smk"
