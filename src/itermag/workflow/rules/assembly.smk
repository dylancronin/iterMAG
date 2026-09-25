rule megahit_assembly:
    input:
        r1 = "iter_{iteration}/reads/reads.1.fq.gz",
        r2 = "iter_{iteration}/reads/reads.2.fq.gz"
    output:
        contigs = "iter_{iteration}/assembly/contigs.fa",
    log:
        "iter_{iteration}/assembly/megahit.log"
    threads:
        config["threads"]
    shell:
        """
        set -euo pipefail
        if [ ! -s {input.r1} ] && [ ! -s {input.r2} ]; then
            echo "Input reads are empty; skipping assembly." > {log}
            mkdir -p iter_{wildcards.iteration}/assembly
            touch {output.contigs}
            exit 0
        fi

        rm -rf iter_{wildcards.iteration}/assembly_tmp
        mkdir -p iter_{wildcards.iteration}/assembly
        rm -f {output.contigs}
        megahit -t {threads} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o iter_{wildcards.iteration}/assembly_tmp \
            --presets meta-large \
            > {log} 2>&1

        mkdir -p iter_{wildcards.iteration}/assembly
        mv iter_{wildcards.iteration}/assembly_tmp/final.contigs.fa {output.contigs}

        rm -rf iter_{wildcards.iteration}/assembly_tmp
        """