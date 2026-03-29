import argparse
from itermag.workflow import run_workflow

def main():
    parser = argparse.ArgumentParser(description="iterMAG: Iterative MAG assembly.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    subparsers.required = True
    
    #run command
    run_parser = subparsers.add_parser("run", help="Run iterMAG.")
    run_parser.add_argument("-1", "--forward", required=True, help="Forward reads.")
    run_parser.add_argument("-2", "--reverse", required=True, help="Reverse reads.")
    run_parser.add_argument("-o", "--output", required=True, help="Output directory.")
    run_parser.add_argument("-t", "--threads", required=True, help="Number of threads.")
    run_parser.add_argument("-g", "--genomes", required=False, help="Genome directory.")
    run_parser.add_argument("-it", "--max_iterations", type=int, default=5, help="Maximum number of iterations (default: 5).")
    run_parser.add_argument("--checkm_db_path", required=True, help="Path to CheckM2 database.")
    
    args = parser.parse_args() #parse user input

    if args.command == "run":
        print(f"Running iterMAG on {args.forward} and {args.reverse} with {args.threads} threads...")
        run_workflow(
            forward=args.forward,
            reverse=args.reverse,
            output=args.output,
            threads=args.threads,
            genomes=args.genomes,
            max_iterations=args.max_iterations,
            checkm_db_path=args.checkm_db_path
        )
    

if __name__ == '__main__':
    main()
