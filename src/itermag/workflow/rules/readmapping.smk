rule map_reads:
    input:
        r1 = "iter_{iteration}/reads/reads.1.fq.gz",
        r2 = "iter_{iteration}/reads/reads.2.fq.gz",
        assembly = "iter_{iteration}/assembly/contigs.fa"
    output:
        coverage = "iter_{iteration}/assembly_mapping/metabat_coverage.tsv"
    threads:
        config["threads"]
    shell:
        """
        coverm contig -t {threads} \
            --methods metabat \
            --reference {input.assembly} \
            -1 {input.r1} \
            -2 {input.r2} \
            --min-read-aligned-percent 75 \
            --min-read-percent-identity 95 \
            -o {output.coverage}
        """