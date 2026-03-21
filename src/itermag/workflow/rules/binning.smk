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
        metabat2 -i {input.contigs} \
        -o iter_{wildcards.iteration}/binning/iter{wildcards.iteration}_bins \
        -m 2500 \
        -a {input.coverage} || true

        touch {output.bins}
        """


#metabat2 -i ../assembly/contigs.fa -o mb2 -m 1500 -a metabat_coverage.tsv