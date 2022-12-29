from cmath import nan
from heapq import merge
import pandas as pd
import os 

#############
#User configurable variables
#############
#scripts_dir is the root directory of the analysis folder which contains the scripts_hiti.py script. All functions and classes
#used in the analysis have been stored in this file.
scripts_dir="/media/data/AtteR/projects/hiti/pipeline_output_reorg/Arc-HITI/"
os.chdir(scripts_dir)
from scripts_hiti import *
#sample directory (inside which all the subdirectories exist)
sample_dir=scripts_dir + "/SATI_mCherry_2-/"
os.chdir(sample_dir)
#path to the analysis folder
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SATI_5p"
export_path=sample_dir = "trimmed_data/"

#where the programs bbduk and starcode are found
# program_path="/media/data/AtteR/Attes_bin"

#############


#SATI 5p
############
read_fwd=True
direc='3p'
lliteral = ' literal=CCATGTTATCCTCCTCGCCC'
rliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
filterlitteral = 'CTCCTCGCCCTTGCTCACCATGGTGGCGCGcctgttAACAGGCTAAGAACTCCTCTGA'
target_sequence = "TTGCTCACCATGGTGGCGCGCCTGTT---AACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"
target_sequence = target_sequence.upper()
############
#read preprocessing for each sample: trim, record read counts before and after trimming, cluster the reads 
#using starcode, calculate percentage, merge into the full dataframe containing read count-and percentage for
#each sample.
df_full=import_reads_process_mini(base_path, target_sequence,filterlitteral,lliteral,rliteral,read_fwd, direc)
df_trim_full=calculate_perc_sd(df_full,3)
result="unaligned/SATI_mCherry_5p_2-.fasta"
save_fasta(result, df_trim_full, target_sequence)

csv_file="/".join(result.split("/")[:-1]) +"/"+ result.split("/")[-1].split(".")[0] + ".csv"
full_df_trim=pd.read_csv(csv_file, index_col=[0])

df_trim_full.iloc[:,-2].sum()

#NT
####################
output_path="aligned/NT/"
result=output_path + "SATI_mCherry_5p_2-_prim.fasta"
aligner(full_df_trim, target_sequence, "align_local2", result, output_path, lliteral, rliteral,3,1,"5p" ,"True")
####################


#AA
#484%3
####################
# corr_frame=1
# result="unaligned/SATI_mCherry_5p_1-.fasta"
# output_html="aligned/AA/SATI_mCherry_5p_1-_AA.html"
# out_csv="aligned/AA/SATI_mCherry_5p_1-_AA.csv"
# translate_NT(result, corr_frame,direc, out_csv)
####################
