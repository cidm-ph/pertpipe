import os
import logging
import subprocess
import sys
import shutil
from scripts import assists

# copy bpertussis from databases to mlst db folder
logging.getLogger().setLevel(logging.INFO)

def pertpipe_setup():
    dependency_list = ["abricate", "spades.py", "mlst", "minimap2", "samtools", "bcftools"]
    bpertussis_db = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pertpipe/databases/bpertussis")

    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)
    if "abricate" in dependency_list:
        assists.check_abricate()

    result = subprocess.run('which mlst', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    if result.returncode == 0:
        default_mlst_path = result.stdout.strip().replace("/bin/mlst", "/db/pubmlst")
    else:
        logging.error(f"Error finding path with 'which': {result.stderr}")
        sys.exit(1)

    if not os.path.isdir(default_mlst_path):
        logging.error(f"Error finding path with 'which': {default_mlst_path}")
        sys.exit(2)

    new_destination = default_mlst_path + "/bpertussis"
    if not os.path.isdir(new_destination):
        logging.info(f"Creating 'bpertussis' folder in {default_mlst_path}")
        os.makedirs(new_destination)
        try:
            shutil.copytree(bpertussis_db, new_destination,dirs_exist_ok=True)
            logging.info(f"Copied {bpertussis_db} to {new_destination}")
        except FileExistsError:
            logging.error(f"Destination folder {new_destination} already exists")
        except PermissionError:
            logging.error(f"Permission denied while copying {bpertussis_db} to {new_destination}")
        except Exception as e:
            logging.error(f"An error occurred while copying {bpertussis_db} to {new_destination}")
    else:
        logging.info(f"{new_destination} already exists")

# run the makemlstdb thingo
    bpertussis_mlst = None
    if os.path.isdir(new_destination):
        file_list = [
            "23SrRNA.tfa",
            "bpertussis.txt",
            "fhaB24005550.tfa",
            "fim2.tfa",
            "fim3.tfa",
            "prn.tfa",
            "ptxA.tfa",
            "ptxB.tfa",
            "ptxC.tfa",
            "ptxD.tfa",
            "ptxE.tfa",
            "ptxP.tfa",
        ]
        for file in file_list:
            if os.path.exists(os.path.join(new_destination, file)):
                logging.info(f"{file} exists.")
            else:
                logging.error(f"{file} does not exist! check if it is missing in the database folder {bpertussis_db}")
                sys.exit(3)

    cmd_list = [
        f"mlst-make_blast_db",
        f"mlst --longlist | grep bpertussis"
    ]
    expected_result = 'bpertussis\tBPagST\tptxP\tptxA\tptxB\tptxC\tptxD\tptxE\tfhaB24005550\tfim2\tfim3\t23SrRNA\n'
    check_existing = subprocess.run(cmd_list[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    if check_existing.returncode == 0 and check_existing.stdout != expected_result:
        for cmd in cmd_list:
            mlst_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
            if mlst_result.returncode == 0:
                logging.info(f"Successfully run {cmd}")
            else:
                logging.error(f"{mlst_result.stderr}")
            bpertussis_mlst = mlst_result.stdout
        
        if bpertussis_mlst == expected_result:
            logging.info(f"Successfully created 'bpertussis' mlst scheme" )
        else:
            logging.critical(f"Unsuccessful creation of MLST scheme")
    elif check_existing.returncode == 0 and check_existing.stdout == expected_result:
        logging.info("Everything looks good, setup of MLST complete")

    # test MLST cos having issues with ptxA being ~33 not 1.
    pert_test = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/H640.fasta")
    cmd = f"mlst --scheme bpertussis {pert_test}"
    mlst_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    if mlst_result.returncode == 0:
        logging.info(f"Successfully run {cmd}")
    else:
        logging.error(f"{mlst_result.stderr}")
        sys.exit()
    bpertussis_mlst = mlst_result.stdout.split("\t")[5]
    expected_result = "ptxA(1)"
    bad_result = "ptxA(~33)"
    if bpertussis_mlst == expected_result:
        logging.info(f"MLST scheme test success")
    elif bpertussis_mlst == bad_result:
        logging.info(f"MLST scheme test failed, this sample should be ptxA1 not ptxA~33. Need to re-copy/re-download database and re-do setup.")
    else:
        logging.critical(f"MLST scheme test failed. Need to re-copy/re-download database and re-do.")

# set up abricate automatically.

pertpipe_setup()