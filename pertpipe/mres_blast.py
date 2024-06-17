import os
import pandas as pd
import logging
from Bio.Blast import NCBIXML
from pertpipe import assists


rrna_seq = os.path.join(os.path.dirname(os.path.abspath(__file__)), "databases/23S_rRNA.fasta")
def mres_detection(assembly, outdir):

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    hit_list = []
    blast_cmd = f"blastn -task megablast -query {assembly} -subject {rrna_seq} -outfmt 6 -out {outdir}/blast_23s.txt"
    assists.run_cmd(blast_cmd)
    blast_cmd_2 = f"blastn -task megablast -query {assembly} -subject {rrna_seq} -outfmt 5 -out {outdir}/blast_23s.xml"
    assists.run_cmd(blast_cmd_2)
    blast_df = pd.read_csv(f"{outdir}/blast_23s.txt", sep="\t", header=None)
    result_count = len(blast_df)
    if result_count > 1:
        filt_blast_df = blast_df[blast_df[3] == 2882]
    elif result_count == 1: 
        filt_blast_df = blast_df
    
    filt_count = len(filt_blast_df)
    for hit in filt_blast_df.iterrows():
        aln_len, perc_id, hit_name = hit[1][3], hit[1][2], hit[1][0]
        hit_list.append(hit_name)
        if perc_id == 100.0 and aln_len == 2882:
            logging.info(f'Full length 23s rRNA detected, with no mutations')
            return []
        elif perc_id == 99.965 and aln_len == 2882:
            logging.info(f'23s rRNA detected, with one mutation')
            xml_handle = open(f"{outdir}/blast_23s.xml")
            blast_xml = NCBIXML.parse(xml_handle)
            mutation_list = mres_position(blast_xml, hit_list)
            logging.info(f"{mutation_list}")
            return mutation_list
        else:
            logging.error(f"Encountered issue, potentially truncated 23S rRNA detected, or error in assembly occurred.")
        
    
def mres_position(blast_xml, hit_list):
    for blast_result in blast_xml:
        accession_id = blast_result.query.split(" ")
        if accession_id[0] in hit_list:
            for alignment in blast_result.alignments:
                for hsp in alignment.hsps:
                    midline = hsp.match
                    match_count = midline.count('|') + midline.count(' ')
                    space_positions = [pos for pos, char in enumerate(midline) if char == ' ']
                    formatted_positions = [
                        f"{hsp.sbjct[pos]}{pos + 1}{hsp.query[pos]}"
                        for pos in space_positions
                    ]
                    return formatted_positions
    
