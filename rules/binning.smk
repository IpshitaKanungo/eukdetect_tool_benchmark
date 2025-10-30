rule binning:
     input:
        contigs = lambda wildcards: (
            f"{result}/{kraken2_result}/{wildcards.sample}/{wildcards.sample}_noneukaryotic.fa"  if wildcards.br == "kraken2" else
            f"{result}/{tiara_result}/{wildcards.sample}/{wildcards.sample}_noneukaryotic.fa"    if wildcards.br == "tiara"   else
            f"{result}/{eukrep_result}/{wildcards.sample}/{wildcards.sample}_noneukaryotic.fa"   if wildcards.br == "eukrep"  else
            f"{result}/{megahit_result}/{wildcards.br}/{wildcards.sample}/{wildcards.sample}.contigs.fa"
        ),
        r1 = lambda wildcards: get_path(wildcards.sample)[0],
        r2 = lambda wildcards: get_path(wildcards.sample)[1]  
     output:
        bin_r1 = temp(f"{result}/{binning_result}/{{br}}/{{sample}}/{{sample}}_1.fastq"),
        bin_r2 = temp(f"{result}/{binning_result}/{{br}}/{{sample}}/{{sample}}_2.fastq"),
        metabat_folder = directory(f"{result}/{binning_result}/{{br}}/{{sample}}/metabat2_bins"),
        maxbin_folder = directory(f"{result}/{binning_result}/{{br}}/{{sample}}/maxbin2_bins"),
        concoct_folder = directory(f"{result}/{binning_result}/{{br}}/{{sample}}/concoct_bins")
     log:f"{result}/{binning_result}/logs/{{sample}}_{{br}}_binning.log"
     conda: "../envs/metawrap.yaml"

     params:
        threads=config["binning"]["threads"],
        tools = " ".join(filter(None,[
            "--metabat2" if config["binning"]["binning_tools"].get("metabat2", False) else "",
            "--maxbin2" if config["binning"]["binning_tools"].get("maxbin2", False) else "",
            "--concoct" if config["binning"]["binning_tools"].get("concoct", False) else ""
        ]))

     shell:
        """
        mkdir -p $(dirname {output.bin_r1})
        if [[ "{input.r1}" == *.gz ]]; then zcat {input.r1} > {output.bin_r1}; else ln -sf $(readlink -f {input.r1}) {output.bin_r1}; fi
        if [[ "{input.r2}" == *.gz ]]; then zcat {input.r2} > {output.bin_r2}; else ln -sf $(readlink -f {input.r2}) {output.bin_r2}; fi


        mkdir -p {result}/{binning_result}/logs
        mkdir -p {result}/{binning_result}/{wildcards.br}/{wildcards.sample}
        mkdir -p {result}/{binning_result}/{wildcards.br}/{wildcards.sample}/metabat2_bins
        mkdir -p {result}/{binning_result}/{wildcards.br}/{wildcards.sample}/maxbin2_bins
        mkdir -p {result}/{binning_result}/{wildcards.br}/{wildcards.sample}/concoct_bins

        metawrap binning -o {result}/{binning_result}/{wildcards.br}/{wildcards.sample} \
        -t {params.threads} -l 1500 -a {input.contigs} \
        {params.tools} {output.bin_r1} {output.bin_r2}  &> {log}
        """

