rule kraken2:
     input:
         assembled_read=f"{result}/{megahit_result}/host/{{sample}}/{{sample}}.contigs.fa"
     output:
         kraken2_output=f"{result}/{kraken2_result}/{{sample}}/{{sample}}_kraken2.krak2",
         kraken2_report=f"{result}/{kraken2_result}/{{sample}}/{{sample}}_taxonomic_report.tsv",
         euk_out=f"{result}/{kraken2_result}/{{sample}}/{{sample}}_eukaryotic.fa",
         noneuk_out=f"{result}/{kraken2_result}/{{sample}}/{{sample}}_noneukaryotic.fa"
     log:
         log=f"{result}/{kraken2_result}/logs/{{sample}}_kraken2.log"
     params:
         threads=config["kraken2"]["threads"],
         db=config["kraken2"]["db"]
     shell: 
         """
         mkdir -p {result}/{kraken2_result}/logs
         mkdir -p {result}/{kraken2_result}/{wildcards.sample}
       
         kraken2 --confidence 0.1 --db {params.db} --use-names --threads {params.threads} \
               --report {output.kraken2_report} \
               --output {output.kraken2_output} \
               {input.assembled_read}&> {log} 
         
         python ../euk_non-euk_segregation_scripts/kraken2_split.py \
          --tax-report "{output.kraken2_report}" \
          --kraken2    "{output.kraken2_output}" \
          --assembly   "{input.assembled_read}" \
          --euk-out    "{output.euk_out}" \
          --noneuk-out "{output.noneuk_out}" \
          >> "{log}" 2>&1
         """
