rule megahit:
     input:
        lambda wildcards: (
            [
                f"{result}/{bowtie2_result}/{wildcards.sample}/{wildcards.sample}_unaligned_1.fq",
                f"{result}/{bowtie2_result}/{wildcards.sample}/{wildcards.sample}_unaligned_2.fq",
            ]
            if wildcards.branch == "bowtie2"
            else get_path(wildcards.sample)   
        )
     output:
         contigs=f"{result}/{megahit_result}/{{branch}}/{{sample}}/{{sample}}.contigs.fa",
         outdir=directory(f"{result}/{megahit_result}/{{branch}}/{{sample}}")
     log:log=f"{result}/{megahit_result}/logs/{{sample}}_{{branch}}_megahit.log"
     params:
         threads=config["megahit"]["threads"],
     shell:
         """
         mkdir -p {result}/{megahit_result}/logs
         mkdir -p {result}/{megahit_result}/{wildcards.branch}
         tmpdir="{output.outdir}.tmp"
         rm -rf "$tmpdir"
         megahit -1 {input[0]} -2 {input[1]}\
                 -t {params.threads} \
                 --out-prefix {wildcards.sample} \
                 -o "$tmpdir" \
                 > {log} 2>&1
         rm -rf "{output.outdir}"
         mv "$tmpdir" "{output.outdir}"
         """
