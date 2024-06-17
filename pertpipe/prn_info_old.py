import os
import logging
import pandas as pd
from Bio.Blast import NCBIXML
from pertpipe import assists

prn_seq = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases/IR_PRN.fasta")
prn_type_seq = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases/bpertussis/prn.tfa") # all the types
is_elements = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases/IS_elements.fasta") # IS elements.

def match_known_prn(is_string, prn_type, prn_cut_position):
    mut_name = None
    known_prn = {
       # "146-2733/2733" : "prn2::del(-1513,145)", #H646 prn2
        ('?', 'prn2', -42) : "prn2::del(-283, -40)", #J078
        ('?', 'prn28', -331) : "prn2::del(666, 667)", #J473
        ('?', 'prn24', -331) :"prn2::Stop-C1273T", # H697
        ('?', 'prn1', 1326): "prn2::del(−292, 1340)", # 24-013-0032 & J625 
        ('IS481', 'prn2', 1613) : "prn2::IS481(1613)", # most common IS481 insertion & 24-013-0032
        ('IS481', 'prn2', 245) : "prn2::IS481(240)", #H378 & 17-0520-4681
        ('IS481', 'prn9', 1613) : "prn9::IS481(1613)", # J038 
        ('?', 'prn2', '-1185C') : "prn2::ins1185G", #H883
        ('?', 'prn22', ''): "prn2::Stop-C233T", #J018
        ('?', 'prn2', -73): "prn::dis(-73)", #H806 <- just says "promoter disruption"
        ('?', 'prn2', 146): "prn2::del(−1513, 145)" # H696
        # : "prn2::Stop(C739T)",
    }
    mut_name = known_prn.get((is_string, prn_type, prn_cut_position), f"Unknown PRN mutation {prn_type}::{is_string}({prn_cut_position}), Please investigate :)")
    return mut_name

def blast_type_filter(blast_type):
    aln_length = None
    blast_type[2] = blast_type[2].astype('float64')
    blast_type_filter = blast_type[(blast_type[2].isin([100, 100.0, 100.00, 100.000]))]
    if blast_type[3].any() > 2700:
        aln_length = "Full"
        blast_type_filter = blast_type[(blast_type[2].isin([100, 100.0, 100.00, 100.000])) & (blast_type[3] > 2700)]
    elif blast_type[3].any() < 2700:
        aln_length = "Partial"
    if len(blast_type_filter) == 0:
        blast_type_filter = blast_type[(blast_type[2] > 90)]
        if blast_type[3].any() > 2700:
            aln_length = "Full"
            blast_type_filter = blast_type[(blast_type[2] > 90) & (blast_type[3] > 2700)]
        elif blast_type[3].any() < 2700:
            aln_length = "Partial"
    #blast_type[(blast_type[2] == 100) | (blast_type[2] == 100.000)]
    unique_types = len(blast_type_filter[1].drop_duplicates())
    if len(blast_type_filter) == 1 and unique_types > 0 and aln_length == "Full":
        prn_type = blast_type_filter.loc[blast_type_filter[11].idxmax()][1]
        prn_row = blast_type[blast_type[1].str.match(prn_type)]
    elif len(blast_type_filter) >= 2 and unique_types < 1 and aln_length == "Partial":
        blast_type_filter = blast_type_filter[blast_type_filter[1].map(blast_type_filter[1].value_counts()) == 2]
        prn_type = blast_type_filter[1]
        prn_row = blast_type[blast_type[1] == prn_type]
    elif len(blast_type_filter) >= 2 and unique_types > 1 and aln_length == "Partial":
        grouped = blast_type_filter.groupby(1)[11].sum().reset_index()
        prn_type = grouped.loc[grouped[11].idxmax()][1]
        prn_row = blast_type[blast_type[1] == prn_type]
    else:
        prn_row = None
        prn_type = None
    return prn_row, prn_type.replace("_", "")
        


def prn_mutations(blast_xml, hit_list):
    for blast_result in blast_xml:
        if blast_result.query in hit_list:
            for alignment in blast_result.alignments:
                for hsp in alignment.hsps:
                    midline = hsp.match
                    match_count = midline.count('|') + midline.count(' ')
                    space_positions = [pos for pos, char in enumerate(midline) if char == ' ']
                    formatted_positions = [
                        f"{hsp.sbjct[pos]}{pos}{hsp.query[pos]}"
                        for pos in space_positions
                    ]
                    return formatted_positions
                
def prn_analysis(assembly, prn_outdir):
    # commands needed for prn analysis
    abricate_cmd = f"abricate --datadir {assists.bor_vfdb_db} --db bp-only_vfdb --quiet {assembly} > {prn_outdir}/vfdb.txt"
    mlst_cmd = f"mlst --scheme bpertussis {assembly} > {prn_outdir}/mlst.txt"
    blast_cmd = f"blastn -task megablast -query {assembly} -subject {prn_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn.txt"
    blast_cmd_2 = f"blastn -task megablast -query {assembly} -subject {prn_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn.xml"
    blast_cmd_3 = f"blastn -task megablast -query {assembly} -subject {prn_type_seq} -outfmt 6 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.txt"
    blast_cmd_4 = f"blastn -task megablast -query {assembly} -subject {prn_type_seq} -outfmt 5 -min_raw_gapped_score 100 -out {prn_outdir}/blast_prn_type.xml"
    prn_cut_position = None
    is_string = "?"
    # run the commands
    for command in [abricate_cmd, mlst_cmd, blast_cmd, blast_cmd_2, blast_cmd_3, blast_cmd_4]: 
        assists.run_cmd(command)

    # check the outputs
    for outfile in [f"{prn_outdir}/vfdb.txt", f"{prn_outdir}/mlst.txt", f"{prn_outdir}/blast_prn.txt"]:
        assists.check_files(outfile)

    # look at vfdb for cuts first. If full do different things.
    vfdb_info = pd.read_csv(f"{prn_outdir}/vfdb.txt", sep="\t", header=0)
    prn_info = vfdb_info[vfdb_info['GENE'].str.contains('prn')].reset_index()
    prn_info['%IDENTITY'] = prn_info['%IDENTITY'].astype('float64')
    blast_prn = pd.read_csv(f"{prn_outdir}/blast_prn.txt", sep="\t", header=None)
    #naming_prn = blast_prn[(blast_prn[1] == "AJ011092.1") & ((blast_prn[8] == 1) | (blast_prn[9] == 1))]
    promoter_prn = blast_prn[(blast_prn[1] == "NC_002929_IR_PRN")] # & ((blast_prn[8] == 1) | (blast_prn[9] == 1))]
    blast_type = pd.read_csv(f"{prn_outdir}/blast_prn_type.txt", sep="\t", header=None)
    prn_row, prn_type = blast_type_filter(blast_type)

    # extracting the prn contigs if more than 2 occur (likely split) and identifying if the split is caused by an IS element.
    if len(prn_row) > 1:
        is_prn = pd.DataFrame()
        #prn_contigs = extract_contigs(assembly, prn_row, prn_outdir)
        blast_cmd_5 = f"blastn -task megablast -query {prn_contigs} -subject {is_elements} -outfmt 6 -out {prn_outdir}/blast_prn_is.txt"
        assists.run_cmd(blast_cmd_4)
        if os.path.isfile(prn_outdir + "/blast_prn_is.txt") is True and os.stat(prn_outdir + "/blast_prn_is.txt").st_size != 0:
            is_prn = pd.read_csv(f"{prn_outdir}/blast_prn_is.txt", sep="\t", header=None)
        if not is_prn.empty:
            pident_gr_90 = is_prn[2] > 90
            alnlen_gr_25 = is_prn[3] > 25

            filt_is_prn = is_prn[pident_gr_90 & alnlen_gr_25]
            if len(filt_is_prn[1]) >= 2:
                contains_irl = filt_is_prn[1].str.contains("_IRL").any()
                contains_irr = filt_is_prn[1].str.contains("_IRR").any()
                if contains_irl and contains_irr:
                    irl_string = filt_is_prn[filt_is_prn[1].str.contains('_IRL')].iloc[0][1].rstrip('_IRL')
                    irr_string = filt_is_prn[filt_is_prn[1].str.contains('_IRR')].iloc[0][1].rstrip('_IRR')
                    if irl_string == irr_string:
                        is_string = irl_string
                    else:
                        logging.info("IS elements found, however they are different IS elements.")

                
            elif len(is_prn[1]) < 1:
                logging.info("Only one section of the IS element detected, further investigation is required.") 

        prn_row = prn_row[((prn_row[8] == 1) | (prn_row[9] == 1))]
    logging.info(f"PRN Type: {prn_type}")

    # full prn genes
    mutation_list = None
    if len(prn_info) == 1:
        if prn_info['COVERAGE'][0] == "1-2733/2733" and prn_row[2][0] == 100.00:
            logging.info(f"One PRN ({prn_type}) section detected. PRN gene intact. Therefore likely to be Pertactin positive. Investigating promoter region mutations...")
        elif prn_info['COVERAGE'][0] == "1-2733/2733" and prn_row[2][0] != 100.00:
            hit_list = prn_info['SEQUENCE'][0]
            xml_handle = open(f"{prn_outdir}/blast_prn_type.xml")
            blast_xml = NCBIXML.parse(xml_handle)
            mutation_list = prn_mutations(blast_xml, hit_list)
            prn_cut_position = ''.join(mutation_list)
        elif prn_info['COVERAGE'][0] != "1-2733/2733":
            logging.info(f"One PRN ({prn_type}) section detected. PRN truncated. Therefore likely to be Pertactin negative. Investigating promoter region mutations...")
        if not promoter_prn.empty and mutation_list == None:
            row = promoter_prn.iloc[0]
            if row[2] == 100 and row[3] == 2733:
                logging.info(f"No mutations in promoter or PRN regions.")
            if row[8] > row[9]:
                prn_cut_position = row[9] - 332
                logging.info(f"promoter_cut_position: {prn_cut_position}")
            elif row[8] < row[9]:
                prn_cut_position = row[8] - 332
                logging.info(f"promoter_cut_position: {prn_cut_position}")
            else:
                logging.info("Something is not right.")

    # cut prn genes
    
    elif len(prn_info) > 1:
        logging.info(f"PRN in 2 or more sections detected.")
        if not prn_row.empty:
            row = prn_row.iloc[0]
            if row[8] > row[9] and row[9] == 1:
                prn_cut_position = row[8]
                logging.info(f"prn_cut_position: {prn_cut_position}")
            elif row[8] < row[9] and row[8 == 1]:
                prn_cut_position = row[9]
                logging.info(f"prn_cut_position: {prn_cut_position}")
            else:
                logging.info("Something is not right.")

    if prn_cut_position != None:
        mutation_name = match_known_prn(is_string, prn_type, prn_cut_position)
        logging.info(f"Detected prn mutation: {mutation_name}")
    return mutation_name



    
    
