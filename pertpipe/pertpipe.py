import datetime
import logging
import os
import sys
import warnings
import glob
from pertpipe import assists
from pertpipe import arguments
from pertpipe import prn_info
from pertpipe import mres_blast
from pertpipe import mres_copy_no

__version__ = "0.0.1"

warnings.simplefilter(action="ignore", category=FutureWarning)

dependency_list = ["abricate", "spades.py", "mlst", "minimap2", "samtools", "bcftools"]
ref_list = []

import logging
import datetime
import os
import sys

def pertpipe():
    """
    Running order of pertpipe
    """
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    logger = logging.getLogger()
    is_reads = bool(args.R1 is not None)

    # set outdir defaults - if no outdir is set, it will default to either the fasta or R1 location
    if args.outdir is None and args.R1 is not None:
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

    # Ensure no duplicate handlers are added
    if not logger.hasHandlers():
        formatter = logging.Formatter(
            "pertpipe:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
        )
        stdout_handler = logging.StreamHandler(sys.stdout)
        file_handler = logging.FileHandler(errorlog, mode="w+")
        for handler in [stdout_handler, file_handler]:
            handler.setLevel(logging.INFO)
            handler.setFormatter(formatter)
            logger.addHandler(handler)

    logger.setLevel(logging.INFO)

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

    assists.check_files(args.R1)
    assists.check_files(args.R2)

    # checking all the versions and installations of dependencies.
    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)
    if "abricate" in dependency_list:
        assists.check_abricate()
    
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

    # abricate for vaccine antigens
    assembly = spades_outdir + "/contigs.fasta"
    assists.check_files(assembly)

    prn_outdir = maindir + "/analysis"
    folder_exists = os.path.exists(prn_outdir)
    if not folder_exists:
        os.makedirs(prn_outdir)
        logging.info("Making vaccine antigens output folder")
    else:
        logging.info(f"Spades folder exists")
    prn_info.prn_analysis(assembly, prn_outdir)

    # 23s rRNA for macrolide resistance
    analysis_outdir = maindir + "/analysis"
    folder_exists = os.path.exists(analysis_outdir)
    if not folder_exists:
        os.makedirs(analysis_outdir)
        logging.info("Making analysis output folder")
    else:
        logging.info(f"analysis folder exists")
    mutation_list = mres_blast.mres_detection(assembly, analysis_outdir)
    if mutation_list != []:
        mres_copy_no.mres_copy_numbers(args.R1, args.R2, analysis_outdir, mutation_list)



if __name__ == "__main__":
    pertpipe()