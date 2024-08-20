import os
import logging
import pandas as pd
from io import StringIO
from scripts import assists
from scripts import mres_copy_no

rrna_seq = os.path.join(os.path.dirname(os.path.dirname(__file__)), "databases/23S_rRNA.fasta")
def mres_mapping_only(reads1, reads2, outdir):
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
    result_dict = map_calculations(f"{outdir}/mres.vcf")
    return result_dict

def map_calculations(bcftools_vcf):
    with open(bcftools_vcf, "r") as vcf:
        oneline = ""
        lines = vcf.readlines()
        for line in lines:
            if not line.startswith("##"):
                oneline += line
    
    vcf_df = pd.read_csv(StringIO(oneline), sep="\t", header=0)
    if vcf_df.empty is False:
        info_data = vcf_df["INFO"].str.split(";", expand = True)
        mres_df = pd.DataFrame()
        mres_df['DP4'] = info_data[7].str.lstrip('DP4=')
        split_columns = mres_df['DP4'].str.split(',', expand=True)
        mres_df['sum_all_four'] = split_columns.sum(axis=1).astype(int)
        mres_df['sum_last_two'] = split_columns.iloc[:, -2:].sum(axis=1).astype(int)
        mres_df['FREQ'] = mres_df['sum_last_two'] / mres_df['sum_all_four']
        mres_df["REFPOSALT"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"]
        copy_no = mres_df["FREQ"].apply(mres_copy_no.determine_copy_number)[0]
        mres_df = mres_df[['REFPOSALT', 'FREQ']]
        positions = mres_df["REFPOSALT"][0]
        logging.info(f"23S mutation occurs as a {positions}, very likely in {copy_no}")
        result_dict = {
            "Resistance": "Resistant",
            "Mutation": positions,
            "Copy No": copy_no
        }
    else:
        logging.info(f"No 23S mutations were detected in this sample")
        result_dict = {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A"
            }
    return result_dict
