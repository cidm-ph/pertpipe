import os
import re
import logging
import pandas as pd
from Bio.Blast import NCBIXML
from scripts import assists
from scripts import prn_assists
#from scripts import draw_figure

prn_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/IR_PRN.fasta")
prn_type_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/bpertussis/prn.tfa") # all the prn types
fha_type_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/fhaB.fasta") # all the fhaB types
is_elements = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/IS_elements.fasta") # IS elements.

def virulence_analysis(assembly, prn_outdir, closed, datadir, prokka_outdir):
    # commands needed for prn analysis
    abricate_cmd = f"abricate --datadir {assists.bor_vfdb_db} --db bp-only_vfdb --quiet {assembly} > {prn_outdir}/vfdb.txt"
    mlst_cmd = f"mlst --scheme bpertussis {assembly} > {prn_outdir}/mlst.txt"
    mlst_datadir_cmd = f"mlst --scheme bpertussis --datadir {datadir}/pubmlst --blastdb {datadir}/blast/mlst.fa {assembly} > {prn_outdir}/mlst.txt"
    blast_cmds = [
        f"blastn -task megablast -query {assembly} -subject {prn_seq} -outfmt 6 -out {prn_outdir}/blast_prn.txt",
        f"blastn -task megablast -query {assembly} -subject {prn_seq} -outfmt 5 -out {prn_outdir}/blast_prn.xml",
        f"blastn -task megablast -query {assembly} -subject {prn_type_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.txt",
        f"blastn -task megablast -query {assembly} -subject {prn_type_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.xml",
        f"blastn -task megablast -query {assembly} -subject {fha_type_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_fhaB_type.txt",
        f"blastn -task megablast -query {assembly} -subject {fha_type_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_fhaB_type.xml"
    ]
    #prn_cut_position, prn = None
    #is_string = "?"

    # run the commands
    for command in [abricate_cmd] + blast_cmds:
        assists.run_cmd(command)

    if datadir is not None:
        assists.run_cmd(mlst_datadir_cmd)
    else:
        assists.run_cmd(mlst_cmd)

    mandatory_files = [
        f"{prn_outdir}/vfdb.txt", 
        f"{prn_outdir}/mlst.txt"
    ]

    optional_files = [
        f"{prn_outdir}/blast_prn.txt",
        f"{prn_outdir}/blast_prn.xml",
        f"{prn_outdir}/blast_prn_type.txt",
        f"{prn_outdir}/blast_prn_type.xml",
        f"{prn_outdir}/blast_fhaB_type.txt",
        f"{prn_outdir}/blast_fhaB_type.xml"
    ]
    # check the outputs
    for outfile in mandatory_files:
        assists.check_files(outfile)
    for outfile in optional_files:
        if os.path.isfile(outfile) is True and os.stat(outfile).st_size != 0:
            truemsg = outfile + " exists and not empty, continuing..."
            logging.info(truemsg)
       
    # lets check how many contigs contain prn
    try:
        vfdb_info = pd.read_csv(f"{prn_outdir}/vfdb.txt", sep="\t", header=0)
        prn_vfdb = vfdb_info[vfdb_info['GENE'].str.contains('prn')].reset_index()
        prn_type_info = pd.read_csv(f"{prn_outdir}/blast_prn_type.txt", sep="\t", header=None)
        prn_type_xml = open(f"{prn_outdir}/blast_prn_type.xml")
        prn_promoter = pd.read_csv(f"{prn_outdir}/blast_prn.txt", sep="\t", header=None)
        prn_xml = open(f"{prn_outdir}/blast_prn.xml")
    except Exception as e:
        logging.warning(f"Failed to load one or more files: {e}")
        vfdb_info, prn_vfdb, prn_type_info, prn_type_xml, prn_promoter, prn_xml = None, None, None, None, None, None
    file_vars = [vfdb_info, prn_vfdb, prn_type_info, prn_type_xml, prn_promoter, prn_xml]

    #prn_promoter = pd.read_csv(f"{prn_outdir}/blast_type.txt", sep="\t", header=None)


    # this is the 1 PRN gene checking & pathway
    if not any(file is None for file in file_vars):
        if len(prn_vfdb) == 1 or prn_vfdb['COVERAGE'].any() == '1-2733/2733': # sometimes two prn genes can come up even if there is a full length gene.
            logging.info(f"1 PRN gene detected")
            coverage = prn_vfdb['COVERAGE'][0]
            if coverage == '1-2733/2733': # check if coverage is 100.0 so that we know its full length
                logging.info(f"Full length PRN gene detected")
                prn_row, prn_type = prn_assists.prn_type(prn_type_info, "full") # full PRN typing
                is_prn1or2 = bool(re.findall(r'\bprn[12]\b', prn_type)) # check if its prn1 or prn2, if not find snps
                if closed == True:  # check if closed genome
                    prn_promoter_xml = NCBIXML.parse(prn_xml)
                    prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "full") # extract only PRN region in closed genomes.
                # this command only should be used for short read assembled genomes.
                else:
                    prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
                if is_prn1or2 is True:
                    if prn_row.iloc[0][2] != 100.0: # now check if we need to screen for new mutations!
                        blast_prn_xml = NCBIXML.parse(prn_type_xml)
                        mut_type, mutation = prn_assists.snp_mutations(blast_prn_xml, prn_row, prn_type)
                        prn_type = prn_assists.match_known_prn(mut_type, prn_type, mutation, None)
                else:
                    prn_type = prn_assists.promoter_scan(prn_promoter, prn_row, prn_type)
            else:
                logging.info(f"Truncated PRN gene detected")
                prn_row, prn_type = prn_assists.prn_type(prn_type_info, "partial") # partial PRN typing
                if closed == True:
                    prn_promoter_xml = NCBIXML.parse(prn_xml)
                    prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "partial")
                # this command only should be used for short read assembled genomes.
                else:
                    prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
                blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
                assists.run_cmd(blast_cmd_5)
                if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
                    is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
                    prn_type = prn_assists.dupe_type(prn_promoter, prn_row, is_prn, prn_type)
                else:
                    prn_type= prn_assists.dupe_type(prn_promoter, prn_row, None, prn_type)
                

        # this is the 2 PRN genes checking & pathway
        if len(prn_vfdb) > 1:
            is_prn = pd.DataFrame()
            prn_type = None
            logging.info(f"2 or more PRN genes detected, going down PRN deficient analysis")
            prn_row, prn_type = prn_assists.prn_type(prn_type_info, "dupe") # duplicate PRN typing
            if closed == True:
                prn_promoter_xml = NCBIXML.parse(prn_xml)
                prn_contigs = prn_assists.extract_prn(assembly, prn_promoter_xml, prn_promoter, prn_outdir, "dupe")
            else:
                prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
            blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
            assists.run_cmd(blast_cmd_5)
            if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
                is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
                prn_type = prn_assists.dupe_type(prn_promoter, prn_row, is_prn, prn_type)
            else:
                prn_type = prn_assists.dupe_type(prn_promoter, prn_row, None, prn_type)
    else: 
        prn_row = prn_type = "Not Detected"
    if prn_row is None or prn_type is None:
        prn_row = prn_type = "Not Detected"
    # run fhaB checking now.
    try:
        fhaB_type_info = pd.read_csv(f"{prn_outdir}/blast_fhaB_type.txt", sep="\t", header=None)
        fhaB_type_xml = open(f"{prn_outdir}/blast_fhaB_type.xml")
        fhaB_vfdb = vfdb_info[vfdb_info['GENE'].str.contains('fhaB')].reset_index() 
    except Exception as e:
        logging.warning(f"Failed to load one or more files: {e}")
        vfdb_info, fhaB_type_info, fhaB_type_xml = None, None, None
    file_vars = [vfdb_info, fhaB_type_info, fhaB_type_xml]



    if not any(file is None for file in file_vars):
        if len(fhaB_vfdb) == 1:
            max_value = fhaB_type_info[11].max()
            rows_with_max_value = fhaB_type_info[fhaB_type_info[11] == max_value]
            coverage = fhaB_vfdb['COVERAGE'][0]
            if len(rows_with_max_value) == 1 and coverage == '1-10773/10773':
                max_length = rows_with_max_value[3][0]
                if max_length > 10700:
                    fhab_len = "full"
                    logging.info(f"Full length fhaB gene detected")
                else:
                    fhab_len = "abnormal"
                fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
            else:
                logging.info(f"Abnormal number {len(rows_with_max_value)} of fhaB genes detected please investigate.")
                fhaB_type = "Abnormal"
        if len(fhaB_vfdb) > 1:
            if any(fhaB_vfdb["COVERAGE"] == "1-10773/10773"):
                fhaB_full_length = fhaB_vfdb[fhaB_vfdb["COVERAGE"] == "1-10773/10773"].reset_index()
                fhaB_contig_name = fhaB_full_length['SEQUENCE'][0]
                logging.info(f"Full length fhaB gene detected")
                fhab_len = "full"
                row_with_seq_name = fhaB_type_info[fhaB_type_info[0] == fhaB_contig_name]
                max_value = fhaB_type_info[11].max()
                rows_with_max_value = row_with_seq_name[row_with_seq_name[11] == max_value]
                fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
            else:
                logging.info(f"fhaB gene truncated, getting best hit")
                fhab_len = "truncated"
                max_value = fhaB_type_info[11].max()
                rows_with_max_value = fhaB_type_info[fhaB_type_info[11] == max_value].reset_index()
                fhaB_type = prn_assists.fhaB_type(rows_with_max_value, fhab_len)
    else:
        fhaB_type = "Not Detected"
    
    #name = os.path.basename(os.path.dirname(prokka_outdir))
    #prokka_gbk, contig_prokka_map = assists.contig_prokka_tag(assembly, name, prokka_outdir)
    #draw_figure.draw_clinker(prn_type, prn_outdir, prokka_gbk, contig_prokka_map)

    # filling the remaining information from MLST/Virulence gene types
    mlst_info = pd.read_csv(f"{prn_outdir}/mlst.txt", sep="\t", header=None)

    # cutting the ptx toxin genes out.
    ptx_toxin_cols = mlst_info.iloc[0, 4:10]

    # removing the ptx suffix from BCDE so that its ptxA(X)B(X)C(X)D(X)
    #ptx_toxin, ptxp, fim2, fim3 = None
    ptx_toxin = ptx_toxin_cols.iloc[1] + ''.join([col[3:] for col in ptx_toxin_cols.iloc[2:]])
    ptxp, fim2, fim3 = mlst_info[4][0], mlst_info[11][0], mlst_info[12][0]

    # final dictionary to be returned.
    virulence_info = {
        "ptxP": ptxp,
        "ptx_toxin": ptx_toxin,
        "prn": prn_type,
        "fim2": fim2,
        "fim3": fim3,
        "fhaB": fhaB_type
    }
    return virulence_info
   