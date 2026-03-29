import os
import pandas as pd

rule checkm:
    input:
        bininfo = "iter_{iteration}/binning/iter{iteration}_bins.BinInfo.txt"

    output:
        quality_report = "iter_{iteration}/checkm2/quality_report.tsv"
    threads:
        config["threads"]
    conda:
        "../envs/checkm.yaml"
    shell:
        """
        if ls iter_{wildcards.iteration}/binning/*.fa 1> /dev/null 2>&1; then
            checkm2 predict -t {threads} \
                -x fa \
                -o iter_{wildcards.iteration}/checkm2 \
                -i iter_{wildcards.iteration}/binning/ \
                --database_path {config[checkm_db_path]} \
                --threads {threads} \
                --remove_intermediates
        else
            mkdir -p iter_{wildcards.iteration}/checkm2
            echo -e "Name\tCompleteness\tContamination" > {output.quality_report}
        fi
        """

rule filter_checkm:
    input:
        quality_report = "iter_{iteration}/checkm2/quality_report.tsv"
    output:
        filtered_report = "iter_{iteration}/binning_filtered/quality_report.filtered.tsv"
    run:
        def apply_simlink(bin_name):
            source_path = f"iter_{wildcards.iteration}/binning/{bin_name}.fa"
            dest_path = f"iter_{wildcards.iteration}/binning_filtered/{bin_name}.fa"
            if os.path.exists(dest_path):
                os.remove(dest_path)
            os.symlink(os.path.abspath(source_path), dest_path)
        
        df = pd.read_csv(input.quality_report, sep="\t")
        df_filtered = df[(df["Completeness"] > 70) & (df["Contamination"] < 10)]
        os.makedirs(f"iter_{wildcards.iteration}/binning_filtered", exist_ok=True)
        df_filtered["Name"].apply(apply_simlink)
        df_filtered.to_csv(output.filtered_report, sep="\t", index=False)
