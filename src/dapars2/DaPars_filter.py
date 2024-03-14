#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import math
import numpy as np
import scipy as sp

# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector


def DaPars_Filtering(input_file, samples_list, treat_index, control_index, output_file,
                     Coverage_cutoff = 30, Num_least_in_treat=1, Num_least_in_control=1,
                     FDR_cutoff = 0.05, PDUI_cutoff = 0.2, Fold_change_cutoff = 0.58):    
    treat_index = [i-1 for i in treat_index]
    control_index = [i-1 for i in control_index]
    try:
        treat_samples = ",".join(samples_list[i] for i in treat_index)
        control_samples = ",".join(samples_list[i] for i in control_index)
    
        print(f"treat samples: {treat_samples}")
        print(f"control samples: {control_samples}")
    except IndexError:
        print("treat or control index out of range, please check.")
        sys.exit(1)

    output_write = open(output_file,'w')
    num_line = 0
    
    result_dict = {}
    All_P_values = []
    Selected_events_id = []
    All_mean_abundance = []
    
    for line in open(input_file,'r'):
        if num_line > 0:
            fields = line.strip('\n').split('\t')
            group1_coverages = np.zeros(2)
            group2_coverages = np.zeros(2)
            num_group1_pass = 0
            group1_PDUIs = 0
            for i in treat_index:
                curr_long = fields[4+i*3]
                curr_short = fields[5+i*3]
                if curr_long != 'NA' and curr_short != 'NA':
                    curr_long  = float(curr_long)
                    curr_short = float(curr_short)
                    if curr_long + curr_short >= Coverage_cutoff:
                        group1_PDUIs = group1_PDUIs + float(fields[6+i*3])
                        num_group1_pass += 1
                        group1_coverages[0] = group1_coverages[0] + curr_long  
                        group1_coverages[1] = group1_coverages[1] + curr_short
                    else:
                        fields[4+i*3] = 'NA'
                        fields[5+i*3] = 'NA'
                        fields[6+i*3] = 'NA'
            
            
            num_group2_pass = 0
            group2_PDUIs = 0
            for i in control_index:
                curr_long = fields[4+i*3]
                curr_short = fields[5+i*3]
                if curr_long != 'NA' and curr_short != 'NA':
                    curr_long  = float(curr_long)
                    curr_short = float(curr_short)
                    if curr_long + curr_short >= Coverage_cutoff:
                        group2_PDUIs = group2_PDUIs + float(fields[6+i*3])
                        num_group2_pass += 1
                        group2_coverages[0] = group2_coverages[0] + curr_long  
                        group2_coverages[1] = group2_coverages[1] + curr_short
                    else:
                        fields[4+i*3] = 'NA'
                        fields[5+i*3] = 'NA'
                        fields[6+i*3] = 'NA'
            
            
            
            if num_group1_pass >= Num_least_in_treat and num_group2_pass >= Num_least_in_control:
                Final_group_diff = str(group1_PDUIs/num_group1_pass - group2_PDUIs/num_group2_pass)
                
                All_mean_abundance.append([group1_PDUIs/num_group1_pass, group2_PDUIs/num_group2_pass])
                
                fields.append(str(Final_group_diff))
                ratio_val,P_val = sp.stats.fisher_exact([group1_coverages/num_group1_pass,group2_coverages/num_group2_pass])
                
                All_P_values.append(P_val)
                Selected_events_id.append(fields[0])
            else:
                fields.append('NA')
            
            
            result_dict[fields[0]] = fields
                    
        else:
            first_line = line.strip('\n').split('\t')      
            
        num_line += 1
    
    
    # stats = importr('stats')
    # All_p_adjust = stats.p_adjust(FloatVector(All_P_values), method = 'BH')
    All_p_adjust = sp.stats.false_discovery_control(All_P_values, method="bh")
    first_line.extend(['Treat_Mean_PDUI', 'Control_Mean_PDUI', 'P_val','adjusted.P_val','Pass_Filter'])
    output_write.writelines('\t'.join(first_line)+'\n')
    for curr_event_id in result_dict:
        mean_PDUI_group1 = 'NA'
        mean_PDUI_group2 = 'NA'
        curr_P_val = 'NA'
        curr_FDR_val = 'NA'
        Pass_filter = 'N'
        curr_fields = result_dict[curr_event_id]
        if curr_event_id in Selected_events_id:
            sel_ind = Selected_events_id.index(curr_event_id)
            curr_P_val = str(All_P_values[sel_ind])
            curr_FDR_val = str(All_p_adjust[sel_ind])
            
            mean_PDUI_group1 = All_mean_abundance[sel_ind][0]
            mean_PDUI_group2 = All_mean_abundance[sel_ind][1]
            
            
            if float(curr_FDR_val) <= FDR_cutoff and abs(float(curr_fields[-1]))>=PDUI_cutoff and abs(math.log((mean_PDUI_group1+1e-5)/(mean_PDUI_group2+1e-5),2))>=Fold_change_cutoff:
                Pass_filter = 'Y'
        
        curr_fields.append(str(mean_PDUI_group1))
        curr_fields.append(str(mean_PDUI_group2))
        curr_fields.append(curr_P_val)
        curr_fields.append(curr_FDR_val)
        curr_fields.append(Pass_filter)
        
        output_write.writelines('\t'.join(curr_fields) +'\n')

    output_write.close()


def DaPars_Filter_main():
    parser = argparse.ArgumentParser(description="DaPars Filter",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_file", required=True, type=str, help="input file")
    parser.add_argument("-s", "--samples-list", dest="samples_list", nargs="+", required=True, type=str, help="samples list")
    parser.add_argument("-t", "--treat-index", dest="treat_index", nargs="+", required=True, type=int, help="treat sample index, 1 based")
    parser.add_argument("-c", "--control-index", dest="control_index", nargs="+", required=True, type=int, help="control sample index, 1 based")
    parser.add_argument("-o", "--output", dest="output_file", required=True, type=str, help="output file")
    parser.add_argument("--3UTR-Coverage-cutoff", dest="Coverage_cutoff", type=float, default=30, help="3UTR coverage cutoff")
    parser.add_argument("--Num-least-in-treat", dest="Num_least_in_treat", type=int, default=1, help="Num sample pass 3UTR coverage filter least in treat")
    parser.add_argument("--Num-least-in-control", dest="Num_least_in_control", type=int, default=1, help="Num sample pass 3UTR coverage filter least in control")
    parser.add_argument("--FDR-cutoff", dest="FDR_cutoff", type=float, default=0.05, help="FDR cutoff")
    parser.add_argument("--PDUI-cutoff", dest="PDUI_cutoff", type=float, default=0.2, help="PDUI cutoff")
    parser.add_argument("--Fold-change-cutoff", dest="Fold_change_cutoff", type=float, default=0.58, help="Fold change cutoff")
    args = parser.parse_args()
    DaPars_Filtering(*vars(args).values())
    

if __name__ == '__main__':
    DaPars_Filter_main()
