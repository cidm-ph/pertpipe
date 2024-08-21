import os
import re
import logging
import pandas as pd
from Bio.Blast import NCBIXML
from scripts import assists
from scripts import prn_assists
#from scripts import draw_figure

prn_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/IR_PRN.fasta")
prn_type_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/bpertussis/prn.tfa") # all the types
is_elements = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/IS_elements.fasta") # IS elements.

def virulence_analysis(assembly, prn_outdir, closed, datadir, prokka_outdir):
    # commands needed for prn analysis
    abricate_cmd = f"abricate --datadir {assists.bor_vfdb_db} --db bp-only_vfdb --quiet {assembly} > {prn_outdir}/vfdb.txt"
    mlst_cmd = f"mlst --scheme bpertussis {assembly} > {prn_outdir}/mlst.txt"
    mlst_datadir_cmd = f"mlst --scheme bpertussis --datadir {datadir}/pubmlst --blastdb {datadir}/blast/mlst.fa {assembly} > {prn_outdir}/mlst.txt"
    blast_cmd = f"blastn -task megablast -query {assembly} -subject {prn_seq} -outfmt 6 -out {prn_outdir}/blast_prn.txt"
    blast_cmd_2 = f"blastn -task megablast -query {assembly} -subject {prn_seq} -outfmt 5 -out {prn_outdir}/blast_prn.xml"
    blast_cmd_3 = f"blastn -task megablast -query {assembly} -subject {prn_type_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.txt"
    blast_cmd_4 = f"blastn -task megablast -query {assembly} -subject {prn_type_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.xml"
    #prn_cut_position, prn = None
    #is_string = "?"

    # run the commands
    for command in [abricate_cmd, blast_cmd, blast_cmd_2, blast_cmd_3, blast_cmd_4]: 
        assists.run_cmd(command)
    if datadir is not None:
        assists.run_cmd(mlst_datadir_cmd)
    else:
        assists.run_cmd(mlst_cmd)

    # check the outputs
    for outfile in [f"{prn_outdir}/vfdb.txt", f"{prn_outdir}/mlst.txt", f"{prn_outdir}/blast_prn.txt"]:
        assists.check_files(outfile)

    # lets check how many contigs contain prn
    vfdb_info = pd.read_csv(f"{prn_outdir}/vfdb.txt", sep="\t", header=0)
    prn_vfdb = vfdb_info[vfdb_info['GENE'].str.contains('prn')].reset_index()
    prn_type_info = pd.read_csv(f"{prn_outdir}/blast_prn_type.txt", sep="\t", header=None)
    #prn_promoter = pd.read_csv(f"{prn_outdir}/blast_type.txt", sep="\t", header=None)
    prn_type_xml = open(f"{prn_outdir}/blast_prn_type.xml")
    prn_promoter = pd.read_csv(f"{prn_outdir}/blast_prn.txt", sep="\t", header=None)
    prn_xml = open(f"{prn_outdir}/blast_prn.xml")

    # this is the 1 PRN gene checking & pathway
    if len(prn_vfdb) == 1:
        logging.info(f"1 PRN gene detected")
        coverage = prn_vfdb['COVERAGE'][0]
        if coverage == '1-2733/2733': # check if coverage is 100.0 so that we know its full length
            logging.info(f"Full length PRN gene detected")
            prn_row, prn_type = prn_assists.prn_type(prn_type_info, "full") # full PRN typing
            is_prn1or2 = bool(re.findall(r'\bprn[12]\b', prn_type)) # check if its prn1 or prn2, if not find snps
            if closed == True:  # check if closed genome
                prn_promoter_xml = NCBIXML.parse(prn_xml)
                prn_contigs = prn_assists.extract_prn(prn_promoter_xml, prn_promoter, prn_outdir) # extract only PRN region in closed genomes.
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
                prn_contigs = prn_assists.extract_prn(prn_promoter_xml, prn_promoter, prn_outdir)
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
            prn_assists.extract_prn(prn_promoter_xml, prn_promoter, prn_outdir)
        else:
            prn_contigs = prn_assists.extract_contigs(assembly, prn_row, prn_outdir) # extracting the contigs matching the prn
            blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
            assists.run_cmd(blast_cmd_5)
            if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
                is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
                prn_type = prn_assists.dupe_type(prn_promoter, prn_row, is_prn, prn_type)
            else:
                prn_type = prn_assists.dupe_type(prn_promoter, prn_row, None, prn_type)

    name = os.path.basename(os.path.dirname(prokka_outdir))
    prokka_gbk, contig_prokka_map = assists.contig_prokka_tag(assembly, name, prokka_outdir)
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
    }
    return virulence_info
   