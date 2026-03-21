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
        megahit -t {threads} \
        -1 {input.r1} \
        -2 {input.r2} \
        -o iter_{wildcards.iteration}/assembly_tmp \
        --presets meta-large \
        > {log} 2>&1

        mv iter_{wildcards.iteration}/assembly_tmp/final.contigs.fa {output.contigs}
        
        rm -rf iter_{wildcards.iteration}/assembly_tmp
        """
        