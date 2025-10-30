rule eukrep:
     input:
        contig=f"{result}/{megahit_result}/host/{{sample}}/{{sample}}.contigs.fa"
     output:
        refined_euk = f"{result}/{eukrep_result}/{{sample}}/{{sample}}_eukaryotic.fa",
        refined_noneuk = f"{result}/{eukrep_result}/{{sample}}/{{sample}}_noneukaryotic.fa"

     log:f"{result}/{eukrep_result}/logs/{{sample}}_eukrep.log"
     conda: "../envs/eukrep.yaml"
     shell:
        """
        mkdir -p {result}/{eukrep_result}
        mkdir -p {result}/{eukrep_result}/{wildcards.sample}
        mkdir -p {result}/{eukrep_result}/logs
        EukRep -i {input.contig} \
               -o {output.refined_euk} \
               --prok {output.refined_noneuk} &> {log}
        """
