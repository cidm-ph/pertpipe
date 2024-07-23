import datetime
import logging
import os
import sys
import warnings
import glob
from scripts import assists
from scripts import arguments
from scripts import virulence_info
from scripts import mres_blast
from scripts import mres_copy_no

__version__ = "0.0.1"
warnings.simplefilter(action="ignore", category=FutureWarning)
logging.getLogger().setLevel(logging.INFO)
formatter = logging.Formatter(
    "pertpipe:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
)

dependency_list = ["abricate", "spades.py", "mlst", "minimap2", "samtools", "bcftools"]
ref_list = []

def pertpipe(args):
    """
    Running order of pertpipe
    """
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    is_assembly = bool(args.fasta is not None)
    is_reads = bool(args.R1 is not None)

    # set outdir defaults - if no outdir is set, it will default to either the fasta or R1 location
    if args.outdir is None and args.fasta is not None:
        default = os.path.dirname(args.fasta)
        outdir = default
    elif args.outdir is None and args.R1 is not None:
        default = os.path.dirname(args.R1)
        outdir = default
    else:
        outdir = args.outdir
    
    # force creation of new folder within set outdir
    maindir = outdir 
    
    #newdir = maindir + "/bams"
    folder_exists = os.path.exists(maindir)
    if not folder_exists:
        os.makedirs(maindir)
        logging.info("Making output folder")
    else:
        logging.info(f"Folder exists")

    # error log
    errorlog = os.path.join(outdir, "pertpipe_" + date + ".log")

    # Clear existing handlers
    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.handlers.clear()

    stdout_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(errorlog, mode="w+")
    for handler in [stdout_handler, file_handler]:
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # cmd checks
    if is_reads is True:
        if args.R2 is None:
            logging.error("R2 was not provided, please provide the paired reads")
            sys.exit(1)

    # launch line
    logging.info(
        "Launching pertpipe v%s and writing output files to directory %s",
        __version__,
        outdir,
    )
    # checking file integrity and existence of output directory
    if all(item is not None for item in [args.fasta, args.R1, args.R2]):
        assists.check_files(args.R1)
        assists.check_files(args.R2)
        assists.check_files(args.fasta)
        logging.info("Found fasta, R1 and R2, skipping Skesa")

    # checking all the versions and installations of dependencies.
    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)
    if "abricate" in dependency_list:
        assists.check_abricate()
    elif "mlst" in dependency_list:
        assists.check_mlst(args.datadir)
    
    if is_reads and is_assembly is False:
        # assembly
        spades_outdir = maindir + "/spades"
        folder_exists = os.path.exists(spades_outdir)
        if not folder_exists:
            os.makedirs(spades_outdir)
            logging.info("Making spades output folder")
        else:
            logging.info(f"Spades folder exists")

        spades_result = assists.check_spades_finished(spades_outdir)
        if spades_result is False:
            spades = f"spades.py --careful --only-assembler --pe1-1 {args.R1} --pe1-2 {args.R2} -o {maindir}/spades"
            assists.run_cmd(spades)

        else:
            logging.info("Spades has already finished for this sample. Skipping.")
        assembly = spades_outdir + "/contigs.fasta"
        assists.check_files(assembly)
        closed = assists.check_closed_genome(assembly)
    elif is_reads is False and is_assembly:
        assembly = args.fasta
        closed = assists.check_closed_genome(assembly)
        
    prn_outdir = maindir + "/analysis"
    folder_exists = os.path.exists(prn_outdir)
    if not folder_exists:
        os.makedirs(prn_outdir)
        logging.info("Making vaccine antigens output folder")
    else:
        logging.info(f"Spades folder exists")
    final_dict = {
        "Folder": maindir
    }
    res_dict = virulence_info.virlence_analysis(assembly, prn_outdir, closed)
    final_dict.update(res_dict)

    # 23s rRNA for macrolide resistance
    analysis_outdir = maindir + "/analysis"
    mutation_list = mres_blast.mres_detection(assembly, analysis_outdir)
    if mutation_list != [] and is_assembly is False:
        res_dict = mres_copy_no.mres_copy_numbers(args.R1, args.R2, analysis_outdir, mutation_list)
    else:
        res_dict = {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A"
            }
    final_dict.update(res_dict)

    # Extract headers and values
    headers = list(final_dict.keys())
    values = list(final_dict.values())

    tsv_lines = ["\t".join(headers), "\t".join(values)]
    tsv_string = "\n".join(tsv_lines)
    output_path = outdir + "/vir_res.tsv"
    with open(output_path, 'w') as output_file:
        output_file.write(tsv_string)
        logging.info(f"Writing information to {output_path}")
    logging.info(f"Complete!")

if __name__ == "__main__":
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    pertpipe(args)