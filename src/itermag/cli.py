import argparse

from itermag.workflow import DEFAULT_CHECKM_DB, DEFAULT_MYLOASM_ENV, run_workflow
from itermag.utils.config import InputError, validate_read_inputs


def main():
    parser = argparse.ArgumentParser(description="iterMAG: Iterative MAG assembly.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    subparsers.required = True

    run_parser = subparsers.add_parser("run", help="Run iterMAG.")
    run_parser.add_argument("-1", "--forward", help="Forward reads (short mode).")
    run_parser.add_argument("-2", "--reverse", help="Reverse reads (short mode).")
    run_parser.add_argument("-o", "--output", required=True, help="Output directory.")
    run_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads (default: 1).",
    )
    run_parser.add_argument(
        "-g", "--genomes", required=False, help="Genome directory for warm start."
    )
    run_parser.add_argument(
        "-it",
        "--max_iterations",
        type=int,
        default=5,
        help="Maximum number of iterations (default: 5).",
    )
    run_parser.add_argument(
        "--checkm_db_path",
        default=DEFAULT_CHECKM_DB,
        help=(
            "Path to CheckM2 database "
            "(default: {})".format(DEFAULT_CHECKM_DB)
        ),
    )
    run_parser.add_argument(
        "--conda-prefix",
        default=None,
        help="Directory for Snakemake-managed conda envs (default: inside the output directory).",
    )

    run_parser.add_argument("--read-type", choices=("short", "hifi"), default="short",
                            help="Sequencing type (default: short).")
    run_parser.add_argument("--long-reads", help="PacBio HiFi FASTQ, plain or gzipped.")
    run_parser.add_argument("--myloasm-env", default=DEFAULT_MYLOASM_ENV,
                            help="Existing conda prefix containing myloasm (HiFi mode).")

    args = parser.parse_args()
    try:
        validate_read_inputs(args.forward, args.reverse, args.read_type, args.long_reads)
    except InputError as exc:
        run_parser.error(str(exc))

    if args.command == "run":
        print(
            f"Running iterMAG in {args.read_type} mode "
            f"with {args.threads} threads..."
        )
        run_workflow(
            forward=args.forward,
            reverse=args.reverse,
            output=args.output,
            threads=args.threads,
            genomes=args.genomes,
            max_iterations=args.max_iterations,
            checkm_db_path=args.checkm_db_path,
            conda_prefix=args.conda_prefix,
            read_type=args.read_type,
            long_reads=args.long_reads,
            myloasm_env=args.myloasm_env,
        )


if __name__ == "__main__":
    main()
