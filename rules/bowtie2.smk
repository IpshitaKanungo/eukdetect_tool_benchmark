rule bowtie2:
    input:
        R1=lambda wildcards: get_path(wildcards.sample)[0],
        R2=lambda wildcards: get_path(wildcards.sample)[1]
    output:
        R1_unaligned=f"{result}/{bowtie2_result}/{{sample}}/{{sample}}_unaligned_1.fq",
        R2_unaligned=f"{result}/{bowtie2_result}/{{sample}}/{{sample}}_unaligned_2.fq",
        sam=f"{result}/{bowtie2_result}/{{sample}}_aligned.sam",
    log:log=f"{result}/{bowtie2_result}/logs/{{sample}}_bowtie2.log"
    params:
        threads=config["bowtie2"]["threads"],
        index=config["bowtie2"]["index"]
    shell:
         """
         mkdir -p {result}/{bowtie2_result}
         mkdir -p {result}/{bowtie2_result}/{wildcards.sample}
         mkdir -p {result}/{bowtie2_result}/logs
    
         bowtie2 -x {params.index} \
                -1 {input.R1} \
                -2 {input.R2} \
                -S {output.sam} \
                -p {params.threads} \
                --un-conc {result}/{bowtie2_result}/{wildcards.sample}/unaligned \
                &> {log}
    
         mv {result}/{bowtie2_result}/{wildcards.sample}/unaligned.1 {output.R1_unaligned}
         mv {result}/{bowtie2_result}/{wildcards.sample}/unaligned.2 {output.R2_unaligned}
         """

