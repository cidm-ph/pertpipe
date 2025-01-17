import argparse

__version__ = "1.0.0"

def create_parser():
    parser = argparse.ArgumentParser(description="pertpipe", prog="pertpipe")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--R1", help="Path to R1 file")
    parser.add_argument("--R2", help="Path to R2 file")
    parser.add_argument("--fasta", "-f", help="Path to Fasta file")
    parser.add_argument("--datadir", "-d", help="Path to MLST data directory")
    parser.add_argument(
        "--meta", 
        "-m", 
        action="store_true",
        help="For use in targeted or clinical metagenomics, to skip Spades and takes into consideration a metagenomic style input")
    parser.add_argument(
        "--longread",
        "-l",
        action="store_true",
        help="Genome was assembled by long reads",
    )
    parser.add_argument(
        "--threads",
        "-t",
        nargs="?",
        const=4,
        type=int,
        help="Specify number of threads used",
    )
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        help="get pertpipe version",
        version=f"pertpipe v{__version__}",
    )
    """
    parser.add_argument(
        "--parallel",
        "-p",
        action="store_true",
        help="Parallelise the mapping portion if computer allows",
    )
    """
    return parser