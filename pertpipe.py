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
from scripts import mres_map

__version__ = "0.0.1"
warnings.simplefilter(action="ignore", category=FutureWarning)
logging.getLogger().setLevel(logging.INFO)
formatter = logging.Formatter(
    "pertpipe:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
)

dependency_list = ["abricate", "spades.py", "mlst", "minimap2", "samtools", "bcftools", "prokka", "kallisto"]
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
        logging.info("Found fasta, R1 and R2, skipping assembly")

    # checking all the versions and installations of dependencies.
    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)
    if "abricate" in dependency_list:
        assists.check_abricate()
    elif "mlst" in dependency_list:
        assists.check_mlst(args.datadir)
    
    if is_reads and is_assembly is False:
        # spades assembly
        spades_outdir = maindir + "/spades"
        folder_exists = os.path.exists(spades_outdir)
        if not folder_exists:
            os.makedirs(spades_outdir)
            logging.info("Making spades output folder")
        else:
            logging.info(f"Spades folder exists")

        spades_result = assists.check_spades_finished(spades_outdir)
        if spades_result is False and args.meta is False:
            spades = f"spades.py --careful --only-assembler --pe1-1 {args.R1} --pe1-2 {args.R2} -o {maindir}/spades"
            assists.run_cmd(spades)
        elif spades_result is False and args.meta is True:
            spades = f"spades.py --meta --only-assembler --pe1-1 {args.R1} --pe1-2 {args.R2} -o {maindir}/spades"
            assists.run_cmd(spades)
        else:
            logging.info("Spades has already finished for this sample. Skipping.")
        assembly = spades_outdir + "/contigs.fasta"
        assists.check_files(assembly)
        closed = assists.check_closed_genome(assembly)
    
    elif is_assembly:
        assembly = args.fasta
        closed = assists.check_closed_genome(assembly)
    
    if is_reads and is_assembly is False and args.meta is True:
        megahit_outdir = maindir + "/megahit"
        folder_exists = os.path.exists(megahit_outdir)

        megahit_result = assists.check_megahit_finished(megahit_outdir)
        if folder_exists and megahit_result is False:
            logging.info(f"Removing existing failed megahit folder")
            os.rmdir(megahit_outdir)
        
        if megahit_result is False:
            megahit = f"megahit -1 {args.R1} -2 {args.R2} -o {megahit_outdir}"
            assists.run_cmd(megahit)
        else:
            logging.info("Megahit has already finished for this sample. Skipping.")

        megahit = megahit_outdir + "/final.contigs.fa"
        #assists.megahit_assembly_graphs(megahit_outdir)
        
        assists.check_files(megahit)
        #closed = assists.check_closed_genome(assembly)

    prokka_outdir = maindir + "/prokka"
    folder_exists = os.path.exists(prokka_outdir)
    name = os.path.basename(maindir)
    if not folder_exists:
        os.makedirs(prokka_outdir)
        logging.info("Making prokka output folder")
    else:
        logging.info(f"Prokka folder exists")
    prokka_result = assists.check_prokka_finished(prokka_outdir, name)
    if prokka_result is False and args.meta is False:
        logging.info(f"Running Prokka")
        prokka = f"prokka --outdir {prokka_outdir} --force --cpus 8 --prefix {name} --locustag {name} --compliant --gcode 11 {assembly}"
        assists.run_cmd(prokka)
    elif prokka_result is False and args.meta is True:
        logging.info(f"Running Prokka in Metagenomics mode")
        prokka = f"prokka --outdir {prokka_outdir} --force --cpus 8 --prefix {name} --locustag {name} --compliant --gcode 11 {assembly} --metagenome"
        assists.run_cmd(prokka)
    else:
        logging.info("Prokka has already finished for this sample. Skipping.")
    

    prn_outdir = maindir + "/analysis"
    folder_exists = os.path.exists(prn_outdir)
    if not folder_exists:
        os.makedirs(prn_outdir)
        logging.info("Making analysis output folder")
    else:
        logging.info(f"Analysis folder exists")
    
    final_dict = {
        "Folder": maindir
    }

    res_dict = virulence_info.virulence_analysis(assembly, prn_outdir, closed, args.datadir, prokka_outdir)
    final_dict.update(res_dict)

    # 23s rRNA for macrolide resistance
    analysis_outdir = maindir + "/analysis"
    if args.meta is False:
        mutation_list, copies, detected = mres_blast.mres_detection(assembly, analysis_outdir, args.meta)
    else:
        mutation_list, copies, detected = mres_blast.mres_detection(megahit, analysis_outdir, args.meta)
    if is_reads is True:
        res_dict = mres_map.mres_map(args.R1, args.R2, analysis_outdir, mutation_list, args.meta)
    elif is_assembly and closed:
        logging.info(f"Assembly only mode")
        if mutation_list != []:
            positions = ",".join(mutation_list)
            logging.info(f"23s mutation occurs as a {positions} in {copies} copies")
            if "A2037G" in mutation_list:
                res_dict = {
                    "Resistance": "Resistant",
                    "Mutation": positions,
                    "Copy No": f"{str(copies)} copies",
                }
            else: 
                res_dict = {
                    "Resistance": "Mutations in 23S rRNA detected",
                    "Mutation": positions,
                    "Copy No": f"{str(copies)} copies",
            }
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