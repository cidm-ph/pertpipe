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
        samtools_filt_q50 = f"samtools view -b -q 60 {outdir}/23s_aln.sort.bam > {outdir}/23s_aln.q60.sort.bam"
        samtools_index_cmd = f"samtools index {outdir}/23s_aln.q60.sort.bam"
        assists.run_cmd(samtools_sort_cmd)
        assists.run_cmd(samtools_filt_q50)
        assists.run_cmd(samtools_index_cmd)
    if os.path.exists(f"{outdir}/mres.vcf") is False or os.stat(f"{outdir}/mres.vcf").st_size == 0:
        bcftools_cmd = f"bcftools mpileup -f {rrna_seq} {outdir}/23s_aln.q60.sort.bam | bcftools call -mv -Ov -o {outdir}/mres.vcf"
        assists.run_cmd(bcftools_cmd)

    result_dict, mres_df = map_calculations(f"{outdir}/mres.vcf")
    mres_df.to_csv(f"{outdir}/23s_mres.txt", sep="\t", index=False, header=False) 
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
        mres_df['DP4'] = info_data.apply(lambda row: next((cell for cell in row if "DP4=" in cell), None), axis=1).str.lstrip('DP4=')
        split_columns = mres_df['DP4'].str.split(',', expand=True)
        split_columns = split_columns.apply(pd.to_numeric, errors='coerce')
        mres_df['sum_all_four'] = split_columns.sum(axis=1).astype(int)
        mres_df['sum_last_two'] = split_columns.iloc[:, -2:].sum(axis=1).astype(int)
        mres_df['FREQ'] = mres_df['sum_last_two'] / mres_df['sum_all_four']
        
        # extract only V domain
        vdomain_start, vdomain_end = 1918, 2444
        mres_df["REFPOSALT"], mres_df["POS"] = vcf_df["REF"] + vcf_df["POS"].astype(str) + vcf_df["ALT"], vcf_df["POS"]
        mres_df = mres_df[(mres_df["POS"] >= vdomain_start) & (mres_df["POS"] <= vdomain_end)]
        mut_total = len(mres_df.index)
        copy_no = mres_df["FREQ"].apply(mres_copy_no.determine_copy_number).to_list()
        mres_df = mres_df[['REFPOSALT', 'FREQ']]
        positions = mres_df["REFPOSALT"].to_list()
        logging.info(f"23S mutation occurs as a {positions}, very likely in {copy_no}")
        if mres_df["REFPOSALT"].str.contains("A2037G").any():
            result_dict = {
                "Resistance": "Resistant",
                "Mutation": ",".join(positions),
                "Copy No": ",".join(copy_no)
            }
        else:
            result_dict = {
                "Resistance": "Mutations in 23S rRNA V domain detected",
                "Mutation": ", ".join(positions),
                "Copy No": ", ".join(copy_no)
            }

    else:
        logging.info(f"No 23S mutations were detected in this sample")
        result_dict = {
            "Resistance": "Susceptible",
            "Mutation": "N/A",
            "Copy No": "N/A"
            }
    return result_dict, mres_df
