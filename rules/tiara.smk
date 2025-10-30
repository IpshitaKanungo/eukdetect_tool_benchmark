rule tiara:
     input:
        Assembly_tiara=f"{result}/{megahit_result}/host/{{sample}}/{{sample}}.contigs.fa"
     output:
        output_tiara=f"{result}/{tiara_result}/{{sample}}/{{sample}}_tiara.txt",
        log_out = f"{result}/{tiara_result}/{{sample}}/log_out_{{sample}}.txt",
        euk_fa = f"{result}/{tiara_result}/{{sample}}/{{sample}}_eukaryotic.fa",
        noneuk_fa = f"{result}/{tiara_result}/{{sample}}/{{sample}}_noneukaryotic.fa",
     log:f"{result}/{tiara_result}/logs/{{sample}}_tiara.log"
     conda: "../envs/tiara.yaml"
     params:
         # Toggle in config.yaml: tiara: { prefer_second: true/false }
        prefer_second = lambda wildcards: "--prefer-second" if config.get("tiara", {}).get("prefer_second", True) else "",
        
     shell:
        """
        mkdir -p {result}/{tiara_result}/logs
        mkdir -p {result}/{tiara_result}/{wildcards.sample}
        tiara -i {input.Assembly_tiara} -o {output.output_tiara} --min_len 1500 --tf all &> {log}
        mv {result}/{tiara_result}/{wildcards.sample}/log_{wildcards.sample}_tiara.txt {output.log_out}

        python python_scripts/tiara_split.py \
          -a "{input.Assembly_tiara}" \
          -t "{output.output_tiara}" \
          -e "{output.euk_fa}" \
          -n "{output.noneuk_fa}" \
          {params.prefer_second} \
          >> "{log}" 2>&1
        """

