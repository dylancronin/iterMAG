from itermag.utils.fastq import stage_gzipped_reads


rule symlink_reads:
    input:
        r1=config["forward"],
        r2=config["reverse"]
    output:
        r1="iter_{iteration}/reads/reads.1.fq.gz",
        r2="iter_{iteration}/reads/reads.2.fq.gz"
    run:
        stage_gzipped_reads(input.r1, output.r1)
        stage_gzipped_reads(input.r2, output.r2)
