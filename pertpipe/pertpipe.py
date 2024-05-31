import datetime
import logging
import os
import sys
import warnings
import glob
import subprocess
from pertpipe import assists
from pertpipe import arguments
from pertpipe import abricate_info
from pertpipe import mres_blast

__version__ = "0.0.1"
logging.getLogger().setLevel(logging.INFO)
warnings.simplefilter(action="ignore", category=FutureWarning)
formatter = logging.Formatter(
    "pertpipe:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
)

dependency_list = ["abricate", "spades.py", "mlst"]
ref_list = []

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

    # error log
    errorlog = os.path.join(outdir + "/pertpipe_" + date + ".log"
                            )
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
    # force creation of new folder within set outdir
    maindir = outdir 
    
    #newdir = maindir + "/bams"
    folder_exists = os.path.exists(maindir)
    if not folder_exists:
        os.makedirs(maindir)
        logging.info("Making output folder")
    else:
        logging.info(f"Folder exists")

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

    if assists.check_spades_finished is False:
        spades = f"spades.py --careful --only-assembler --pe1-1 {args.R1} --pe1-2 {args.R2} -o {maindir}/spades"
        assists.run_cmd(spades)
    else:
        logging.info("Spades has already finished for this sample. Skipping.")

    # abricate for vaccine antigens
    assembly = spades_outdir + "/contigs.fasta"
    assists.check_files(assembly)

    abricate_outdir = maindir + "/abricate"
    abricate_outfile = abricate_outdir + "/vfdb.txt"
    folder_exists = os.path.exists(abricate_outdir)
    if not folder_exists:
        os.makedirs(abricate_outdir)
        logging.info("Making spades output folder")
    else:
        logging.info(f"Spades folder exists")
    if os.path.exists(abricate_outfile) is True and os.stat(abricate_outfile).st_size != 0:
        logging.info("VFDB results from Abricate already exists. Skipping")
    else:
        abricate = f"abricate --datadir {assists.bor_vfdb_db} --db bp-only_vfdb --quiet {assembly} > {abricate_outfile}"
        assists.run_cmd(abricate)
        
    assists.check_files(abricate_outfile)
    abricate_info.abricate_parse(abricate_outfile)

    # 23s rRNA for macrolide resistance
    analysis_outdir = maindir + "/analysis"
    folder_exists = os.path.exists(analysis_outdir)
    if not folder_exists:
        os.makedirs(analysis_outdir)
        logging.info("Making analysis output folder")
    else:
        logging.info(f"analysis folder exists")
    mres_blast.mres_detection(assembly, analysis_outdir)



if __name__ == "__main__":

    pertpipe()