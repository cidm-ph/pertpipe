import os
import subprocess
import sys
import logging
import shutil
import os.path
import pkg_resources

bor_vfdb_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases")

def run_cmd(command):
    """
    Run commands with error outputs.
    """
    logging.info("Running command: %s", command)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if result.returncode != 0:
        logging.critical("Failed to run command: %s", result.args)
        logging.critical("stdout: %s", result.stdout)
        logging.critical("stderr: %s", result.stderr)
        sys.exit(1)
    return result


def check_files(file):
    """
    Check input files if they exist and have contents
    """

    if os.path.isfile(file) is True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = (
            file
            + " either file is does not exist or is empty, please check files. Exiting."
        )
        logging.critical(msg)
        sys.exit(1)


def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) is True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)


def check_dependencies(cmd_exec):
    cmd_path = shutil.which(cmd_exec)
    vcmd = subprocess.run([cmd_exec, "--version"], capture_output=True, text=True)
    result = vcmd.stdout.splitlines()
    if cmd_exec == "abricate":
        version = " ".join(result).replace("abricate ", "v")
        if pkg_resources.parse_version(version) < pkg_resources.parse_version(
            "0.9.8"
        ):
            logging.critical("Abricate version too old, please upgrade to v1.0.0+")
            sys.exit(1)
    elif cmd_exec == "spades.py":
        version = result[0].replace("SPAdes genome assembler ", "")
    elif cmd_exec == "mlst":
        version = result[0].replace("mlst ", "v")

    if cmd_path is not None:
        msg = "Located " + cmd_exec + " " + version + " in " + cmd_path
        logging.info(msg)
    else:
        msg = cmd_exec + " was not found, please check installation on your device"
        logging.critical(msg)
        sys.exit(1)


def check_abricate():
    result = subprocess.run(
        ["abricate", "--list", "--datadir", bor_vfdb_db ],
        capture_output=True,
        text=True,
    )
    if result.returncode > 0:
        logging.critical("Abricate database is not prepared")
        logging.critical(
            "correct by running:   abricate --setupdb --datadir " + bor_vfdb_db 
        )
        sys.exit(1)
    dbs = [x.split("\t")[0] for x in result.stdout.splitlines()[1:]]
    if any(x not in dbs for x in ["bp-only_vfdb"]):
        logging.critical("unable to find databases")
        sys.exit(1)

def check_spades_finished(spades_outdir):
    result = f"======= SPAdes pipeline finished."
    spades_log = spades_outdir + "/spades.log"
    with open(spades_log) as log:
        if result in log.read():
            return True
        return False

