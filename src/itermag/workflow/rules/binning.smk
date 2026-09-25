rule binning:
    input:
        coverage = "iter_{iteration}/assembly_mapping/metabat_coverage.tsv",
        contigs = "iter_{iteration}/assembly/contigs.fa"
    output:
        bins = "iter_{iteration}/binning/iter{iteration}_bins.BinInfo.txt"
    threads:
        config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p iter_{wildcards.iteration}/binning
        rm -f iter_{wildcards.iteration}/binning/iter{wildcards.iteration}_bins.*

        if metabat2 \
            -i {input.contigs} \
            -o iter_{wildcards.iteration}/binning/iter{wildcards.iteration}_bins \
            -m 2500 \
            -a {input.coverage} \
            -t {threads}; then
            : # metabat2 succeeded
        else
            if ls iter_{wildcards.iteration}/binning/iter{wildcards.iteration}_bins.*.fa 1>/dev/null 2>&1; then
                echo "Warning: metabat2 exited non-zero but produced bins."
            else
                echo "metabat2 produced no bins."
                touch {output.bins}
            fi
        fi
        """