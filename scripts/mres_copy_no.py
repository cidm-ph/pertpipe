import os
import logging
import pandas as pd
from io import StringIO
from scripts import assists

rrna_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23S_rRNA.fasta")
def mres_copy_numbers(reads1, reads2, outdir, mutation_list):
    if os.path.exists(f"{outdir}/23s_aln.sort.bam") is False or os.stat(f"{outdir}/23s_aln.sort.bam").st_size == 0:
        minimap2_cmd = f"minimap2 -ax sr {rrna_seq} {reads1} {reads2} | samtools view -b > {outdir}/23s_aln.bam"
        assists.run_cmd(minimap2_cmd)
        samtools_sort_cmd = f"samtools sort {outdir}/23s_aln.bam > {outdir}/23s_aln.sort.bam"
        samtools_index_cmd = f"samtools index {outdir}/23s_aln.sort.bam"
        assists.run_cmd(samtools_sort_cmd)
        assists.run_cmd(samtools_index_cmd)
    if os.path.exists(f"{outdir}/mres.vcf") is False or os.stat(f"{outdir}/mres.vcf").st_size == 0:
        bcftools_cmd = f"bcftools mpileup -f {rrna_seq} {outdir}/23s_aln.sort.bam | bcftools call -mv -Ov -o {outdir}/mres.vcf"
        assists.run_cmd(bcftools_cmd)

    res_dict = vcf_calc_and_blast_match(f"{outdir}/mres.vcf", mutation_list)

    return res_dict

def vcf_calc_and_blast_match(bcftools_vcf, mutation_list):
    
    with open(bcftools_vcf, "r") as vcf:
        oneline = ""
        lines = vcf.readlines()
        for line in lines:
            if not line.startswith("##"):
                oneline += line
    
    vcf_df = pd.read_csv(StringIO(oneline), sep="\t", header=0)

    info_data = vcf_df["INFO"].str.split(";", expand = True)
    mres_df = pd.DataFrame()
    mres_df['DP4'] = info_data[7].str.lstrip('DP4=')
    split_columns = mres_df['DP4'].str.split(',', expand=True)
    mres_df['sum_all_four'] = split_columns.sum(axis=1).astype(int)
    mres_df['sum_last_two'] = split_columns.iloc[:, -2:].sum(axis=1).astype(int)
    mres_df['FREQ'] = mres_df['sum_last_two'] / mres_df['sum_all_four']
    mres_df["REFPOSALT"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"]
    copy_no = mres_df["FREQ"].apply(determine_copy_number)[0]
    mres_df = mres_df[['REFPOSALT', 'FREQ']]
    mask = mres_df['REFPOSALT'].isin(mutation_list)
    mres_df = mres_df[mask]
    positions = ",".join(mutation_list)
    logging.info(f"23s mutation occurs as a {positions}, very likely in {copy_no}")
    result_dict = {
        "Resistance": "Resistant",
        "Mutation": positions,
        "Copy No": copy_no
    }
    return result_dict

def determine_copy_number(freq):
    if 0.68 <= freq <= 1:
        return "3 copies"
    elif 0.05 <= freq <= 0.316:
        return "1 copy"
    else:
        return "2 copies"
