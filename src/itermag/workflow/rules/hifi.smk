from pathlib import Path

from itermag.utils.fastq import stage_gzipped_reads


rule stage_hifi_reads:
    input:
        reads=config["long_reads"]
    output:
        reads="iter_{iteration}/reads/reads.fq.gz"
    run:
        stage_gzipped_reads(input.reads, output.reads)


if not warm_start:
    rule myloasm_assembly:
        input:
            reads="iter_{iteration}/reads/reads.fq.gz"
        output:
            contigs="iter_{iteration}/assembly/contigs.fa"
        params:
            prefix=config["myloasm_env"],
            assembly_dir="iter_{iteration}/assembly/myloasm"
        log:
            "iter_{iteration}/assembly/myloasm.log"
        threads:
            config["threads"]
        shell:
            r"""
            set -euo pipefail
            mkdir -p iter_{wildcards.iteration}/assembly
            conda run --no-capture-output --prefix {params.prefix:q} \
                myloasm {input.reads:q} --hifi -t {threads} \
                -o {params.assembly_dir:q} > {log:q} 2>&1
            test -f {params.assembly_dir:q}/assembly_primary.fa
            cp {params.assembly_dir:q}/assembly_primary.fa {output.contigs:q}
            """

    rule map_hifi_assembly:
        input:
            reads="iter_{iteration}/reads/reads.fq.gz",
            assembly="iter_{iteration}/assembly/contigs.fa"
        output:
            coverage="iter_{iteration}/assembly_mapping/metabat_coverage.tsv"
        threads:
            config["threads"]
        shell:
            r"""
            set -euo pipefail
            mkdir -p iter_{wildcards.iteration}/assembly_mapping
            if [ ! -s {input.assembly:q} ]; then
                printf 'contigName\tcontigLen\ttotalAvgDepth\treads.bam\treads.bam-var\n' \
                    > {output.coverage:q}
            else
                coverm contig -t {threads} \
                    --methods metabat \
                    --reference {input.assembly:q} \
                    --single {input.reads:q} \
                    --mapper minimap2-hifi \
                    --min-read-aligned-percent 75 \
                    --min-read-percent-identity 95 \
                    -o {output.coverage:q}
            fi
            """


rule extract_hifi_reads:
    input:
        quality_report="iter_{iteration}/binning_filtered/quality_report.filtered.tsv",
        reads="iter_{iteration}/reads/reads.fq.gz"
    output:
        bam="iter_{iteration}/bin_mapping/coverm.filtered.bam",
        reads="iter_{iteration}/bin_mapping/reads.fq.gz"
    params:
        bin_dir="iter_{iteration}/binning_filtered",
        mapping_dir="iter_{iteration}/bin_mapping"
    threads:
        config["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.mapping_dir:q}
        shopt -s nullglob
        bins=({params.bin_dir:q}/*.fa)
        reference={params.mapping_dir:q}/all_mag_contigs.fna
        : > "$reference"
        if (( ${{#bins[@]}} )); then
            cat "${{bins[@]}}" > "$reference"
        fi

        if [ ! -s "$reference" ]; then
            # No MAGs: preserve reads, rather than fabricating an empty residual set.
            cp {input.reads:q} {output.reads:q}
            printf '@HD\tVN:1.6\tSO:unsorted\n' \
                | samtools view -b -o {output.bam:q} -
        else
            rm -rf {params.mapping_dir:q}/coverm
            coverm make -t {threads} \
                -r "$reference" \
                --single {input.reads:q} \
                --mapper minimap2-hifi \
                -o {params.mapping_dir:q}/coverm
            bams=({params.mapping_dir:q}/coverm/*.bam)
            if (( ${{#bams[@]}} != 1 )); then
                echo "Expected exactly one HiFi mapping BAM" >&2
                exit 1
            fi
            mv "${{bams[0]}}" {params.mapping_dir:q}/coverm.bam
            rm -rf {params.mapping_dir:q}/coverm
            coverm filter -t {threads} \
                -b {params.mapping_dir:q}/coverm.bam \
                -o {output.bam:q} \
                --inverse \
                --exclude-supplementary \
                --min-read-percent-identity 95 \
                --min-read-aligned-percent 75
            samtools fastq -@ {threads} -F 2304 -n {output.bam:q} \
                | gzip -c > {output.reads:q}
        fi
        """
