import pandas as pd
import logging

def abricate_parse(abricate_outfile):
    vfdb_info = pd.read_csv(abricate_outfile, sep="\t", header=0)
    prn_info = prn_analysis(vfdb_info)
    print(prn_info)

def prn_analysis(vfdb_info):
    known_prn = [
        "ins-IS481::245",
        "ins-IS481::1598",
        "ins-IS481::1613",
        "ins-IS1002::1613",
        "ins-IS481::2745",
        "stop::760",
        "stop::1273",
        "ins-G:1185",
    ]
    
    prn_info = vfdb_info[vfdb_info['GENE'].str.contains('prn')]
    if len(prn_info) < 1:
        if prn_info['COVERAGE'][0] == "1-2733/2733":
            logging.info("One PRN section detected. PRN gene intact. Therefore likely to be Pertactin positive.")
    elif len(prn_info) >= 2:
        if prn_info['COVERAGE'].str.contains("1-1598/2733").any():
            logging.info("Two PRN sections detected. PRN gene disrupted. Therefore likely to be Pertactin negative.")
        #elif prn_info['COVERAGE'].str.contains("")
    if prn_info.empty:
        logging.info(f"No PRN detected")
    return prn_info