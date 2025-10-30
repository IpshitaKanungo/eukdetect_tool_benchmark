rule refinement:
    input:
        metabat2_bins =f"{result}/{binning_result}/{{br}}/{{sample}}/metabat2_bins",
        maxbin2_bins =f"{result}/{binning_result}/{{br}}/{{sample}}/maxbin2_bins",
        concoct_bins =f"{result}/{binning_result}/{{br}}/{{sample}}/concoct_bins"

    output:
        ref_folder = directory(f"{result}/{refinement_result}/{{br}}/{{sample}}")
    log:f"{result}/{refinement_result}/logs/{{sample}}_{{br}}_refinement.log",
    conda: "../envs/metawrap.yaml"
    params:
        threads = config["refinement"]["threads"],
        completeness = config["refinement"]["completeness"],
        contamination = config["refinement"]["contamination"],
        use_metabat2 = config["binning"]["binning_tools"]["metabat2"],
        use_maxbin2 = config["binning"]["binning_tools"]["maxbin2"],
        use_concoct = config["binning"]["binning_tools"]["concoct"]
    shell:
        """
        mkdir -p {result}/{refinement_result}/logs 
        mkdir -p {result}/{refinement_result}/{wildcards.br}/{wildcards.sample}

         # Constructing the command
        command="metawrap bin_refinement -o {result}/{refinement_result}/{wildcards.br}/{wildcards.sample} -t {params.threads}"

        # Add input directories only if they are valid
        flag_count=0
        if [ "{params.use_metabat2}" == "True" ] && [ -n "$(ls -A {input.metabat2_bins})" ]; then
            command+=" -A {input.metabat2_bins}"
            flag_count=$((flag_count+1))
        elif [ "{params.use_metabat2}" == "True" ]; then
            echo "Directory {input.metabat2_bins} is empty." >> {log}
        fi

        if [ "{params.use_maxbin2}" == "True" ] && [ -n "$(ls -A {input.maxbin2_bins})" ]; then
            if [ $flag_count -eq 0 ]; then
                command+=" -A {input.maxbin2_bins}"
            elif [ $flag_count -eq 1 ]; then
                command+=" -B {input.maxbin2_bins}"
            fi
            flag_count=$((flag_count+1))
        elif [ "{params.use_maxbin2}" == "True" ]; then
            echo "Directory {input.maxbin2_bins} is empty." >> {log}
        fi

        if [ "{params.use_concoct}" == "True" ] && [ -n "$(ls -A {input.concoct_bins})" ]; then
            if [ $flag_count -eq 0 ]; then
                command+=" -A {input.concoct_bins}"
            elif [ $flag_count -eq 1 ]; then
                command+=" -B {input.concoct_bins}"
            elif [ $flag_count -eq 2 ]; then
                command+=" -C {input.concoct_bins}"
            fi
            flag_count=$((flag_count+1))
        elif [ "{params.use_concoct}" == "True" ]; then
            echo "Directory {input.concoct_bins} is empty." >> {log}
        fi


        # Check if command has options
        if [[ $command == *"-A"* ]]; then
            # Execute command if at least one option is present
            command+=" -c {params.completeness} -x {params.contamination}"
            echo "Running: $command" >> {log}
            $command &>> {log} || true
        else
            echo "No non-empty input directories found. Skipping bin refinement." >> {log}
        fi
        """

