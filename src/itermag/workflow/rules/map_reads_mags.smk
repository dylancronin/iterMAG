rule extract_reads:
    input:
        quality_report = "iter_{iteration}/binning_filtered/quality_report.filtered.tsv"
    output:
        bam = "iter_{iteration}/bin_mapping/coverm.filtered.bam",
        reads_unmapped_1 = "iter_{iteration}/bin_mapping/reads.1.fq.gz",
        reads_unmapped_2 = "iter_{iteration}/bin_mapping/reads.2.fq.gz"
    threads:
        config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p iter_{wildcards.iteration}/bin_mapping

        cat iter_{wildcards.iteration}/binning_filtered/*.fa > iter_{wildcards.iteration}/bin_mapping/all_mag_contigs.fna 2>/dev/null || : 

        if [ ! -s iter_{wildcards.iteration}/bin_mapping/all_mag_contigs.fna ]; then
            echo "No high-quality MAGs found. Skipping read mapping."
            touch {output.bam} {output.reads_unmapped_1} {output.reads_unmapped_2}
        else
            rm -rf iter_{wildcards.iteration}/bin_mapping/coverm
            coverm make \
                -t {threads} \
                -r iter_{wildcards.iteration}/bin_mapping/all_mag_contigs.fna \
                -1 {config[forward]} \
                -2 {config[reverse]} \
                -o iter_{wildcards.iteration}/bin_mapping/coverm

            mv iter_{wildcards.iteration}/bin_mapping/coverm/*.bam iter_{wildcards.iteration}/bin_mapping/coverm.bam
            rm -rf iter_{wildcards.iteration}/bin_mapping/coverm

            coverm filter \
                -t {threads} \
                -b iter_{wildcards.iteration}/bin_mapping/coverm.bam \
                -o {output.bam} \
                --inverse \
                --min-read-percent-identity-pair 95 \
                --min-read-aligned-percent-pair 75 \
                --proper-pairs-only

            samtools fastq \
                -@ $(( {threads} > 1 ? {threads} - 1 : 1 )) \
                {output.bam} \
                -1 {output.reads_unmapped_1} \
                -2 {output.reads_unmapped_2} \
                -0 /dev/null \
                -s /dev/null \
                -n

            if [ ! -s {output.reads_unmapped_1} ]; then
                echo "All input reads mapped to MAGs. No unmapped reads remain."
                rm -f {output.reads_unmapped_1} {output.reads_unmapped_2}
                touch {output.reads_unmapped_1} {output.reads_unmapped_2}
            fi
        fi
        """