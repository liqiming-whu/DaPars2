import argparse
import os
import pandas as pd

def load_dapars2_res(chromosome, new_header, dir_pre, file_pre):
    input_file = os.path.join(f"{dir_pre}_{chromosome}", f"{file_pre}_result_temp.{chromosome}.txt")
    dap_res = pd.read_table(input_file, sep="\t", header=0, names=new_header)
    return dap_res
    

def merge_apa_res(dir_prefix, file_prefix, chr_list, sample_list, output_file):
    chrfo = open(chr_list)
    chrs = [i.strip() for i in chrfo.readlines()]
    chrfo.close()
    col_names = ["Gene","fit_value","Predicted_Proximal_APA","Loci",*sample_list]
    dap_res_list = []
    for chr in chrs:
        dap_res = load_dapars2_res(chr, col_names, dir_prefix, file_prefix)
        dap_res_list.append(dap_res)
    res_df = pd.concat(dap_res_list)
    res_df.to_csv(output_file, sep="\t", index=False)
    


def merge_apa_res_main():
    parser = argparse.ArgumentParser(description="merge_apa_res",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dir-prefix", dest="dir_prefix", required=True, type=str, default="Dapars2_out", help="directory prefix")
    parser.add_argument("--file-prefix", dest="file_prefix", required=True, type=str, default="Dapars2", help="file prefix")
    parser.add_argument("--chr-list", dest="chr_list", required=True, type=str, help="chromosome list file")
    parser.add_argument("--sample-list", dest="sample_list", nargs="+", type=str, help="sample list")
    parser.add_argument("--output-file", dest="output_file", required=True, type=str, help="output file")
    args = parser.parse_args()
    merge_apa_res_main(*vars(args).values())
