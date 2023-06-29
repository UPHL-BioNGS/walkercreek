#!/bin/python

##############################################
# Modified from script written by Erin Young #
# for creating summary file for grandeur #
##############################################

import pandas as pd
import json
from os.path import exists

##########################################
# defining files                         #
##########################################

# input files
flu_names                    = 'input_files.csv'
abricate_flu                 = 'abricate_InsaFlu_typing.tsv'
ncbi_sra_human_read_scrubber = 'TOTAL_SPOTS_REMOVED.txt'
bbduk                        = 'phix.stats.tsv'
irma                         = 'irma.typing.tsv'
nextclade                    = 'NEXTCLADE_REPORT.tsv'
qc_report                    = 'qc_report.tsv'
multiqc_json                 = 'multiqc_data.json'
multiqc_stats                = 'multiqc_general_stats.txt'
kraken2                      = 'read_percentages.tsv'

# final files
final          = 'flu_summary'

##########################################
# grouping similar files                 #
##########################################

csv_files = [ ]
tsv_files = [abricate_flu, bbduk, irma, kraken2, nextclade, qc_report]
top_hit   = [ ]

##########################################
# exiting if no input files              #
##########################################

if not exists(flu_names) :
    print("No analyses to report on for this run!")
    with open(final + '.tsv', 'w') as fp:
        pass
    with open(final + '.txt', 'w') as fp:
        pass
    quit()

##########################################
# creating the summary dataframe         #
##########################################

summary_df = pd.read_csv(flu_names, dtype = str)
summary_df['warnings'] = ''
columns = list(summary_df.columns)

# csv files
for file in csv_files :
    if exists(file) :
        print("Adding results for " + file)
        analysis = str(file).split("_")[0]
        new_df = pd.read_csv(file, dtype = str, index_col= False)
        new_df = new_df.add_prefix(analysis + "_")
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
        summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# csv files
for file in tsv_files :
    if exists(file) :
        print("Adding results for " + file)
        analysis = str(file).split("_")[0]
        new_df = pd.read_table(file, dtype = str, index_col= False)
        new_df = new_df.add_prefix(analysis + "_")
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
        summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# extracting top hit for multiline result
for file in top_hit :
    if exists(file) :
        print("Adding results for " + file)
        analysis = str(file).split("_")[0]
        new_df = pd.read_csv(file, dtype = str, index_col= False)
        new_df = new_df.drop_duplicates(subset=['sample'], keep='first')
        new_df = new_df.add_prefix(analysis + "_")
        summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
        summary_df.drop(analysis + "_sample", axis=1, inplace=True)

# multiqc : bbduk

if exists(multiqc_json) :
    file = multiqc_data
    print("Adding analysis parsed via multiqc in " + file)
    with open(file) as multiqc_data:
        data = json.load(multiqc_data)

        # bbduk phix reads
        if "bbmap" in data['report_saved_raw_data'].keys():
            print("Adding in phix reads from bbmap")
            samples = [sample.replace(".phix", "") for sample in data['report_saved_raw_data']['bbmap']['stats'].keys()]
            phix_reads=[]
            for sample in data['report_saved_raw_data']['bbmap']['stats'].keys() :
                phix_reads.append(data['report_saved_raw_data']['bbmap']['stats'][sample]['kv']['Matched'])
            bbduk_phixreads_df = pd.DataFrame(samples, columns=['bbduk_sample'])
            bbduk_phixreads_df['bbduk_phix_reads'] = phix_reads
            summary_df = pd.merge(summary_df, bbduk_phixreads_df, left_on="sample", right_on="bbduk_sample", how = 'left')
            summary_df.drop("bbduk_sample", axis=1, inplace=True)

if exists(multiqc_stats) :
    file = multiqc_general_stats
    print("Adding analysis parsed via multiqc in " + file)
    new_df = pd.read_table(file, dtype = str, index_col= False)
    if "FastQC_mqc-generalstats-fastqc-avg_sequence_length" in new_df.columns :
        tmp_df = new_df[["Sample","FastQC_mqc-generalstats-fastqc-avg_sequence_length"]].copy()
        tmp_df["fastqc_avg_length"] = tmp_df["FastQC_mqc-generalstats-fastqc-avg_sequence_length"]
        tmp_df.drop("FastQC_mqc-generalstats-fastqc-avg_sequence_length", axis=1, inplace=True)

        summary_df["possible_fastqc_name"] = summary_df['file'].str.split(" ").str[0].str.split(".").str[0]
        summary_df = pd.merge(summary_df, tmp_df, left_on="possible_fastqc_name", right_on="Sample", how = 'left')
        summary_df.drop("Sample", axis=1, inplace=True)
        summary_df.drop("possible_fastqc_name", axis=1, inplace=True)



##########################################
# creating files                         #
##########################################

summary_df.columns = summary_df.columns.str.replace(' ', '_')

# reducing to the top 1 or 2 results for each analysis
final_columns = [
# general information
'abricate_Insaflu_type',
'abricate_Insaflu_subtype',
'IRMA_type',
'IRMA_subtype',
'nextclade_clade',
'reads_before_trimming',
'reads_after_trimming',
'paired_reads_after_trimming',
'unpaired_reads_after_trimming',
'GC_before_trimming',
'GC_after_trimming',
'average_Q_score_before_trimming',
'average_Q_score_after_trimming',
'kraken2_%_reads_Homo_sapiens',
'kraken2_%_reads_Influenza_A_virus',
'kraken2_%_reads_Influenza_B_virus',
'human_read_scrubber_spots_masked',]

set_columns = []
for new_column in final_columns :
    if new_column in summary_df.columns :
        set_columns.append(new_column)

summary_df.to_csv(final + '.tsv', columns = ['sample','file','version'] + set_columns, index=False, sep="\t")
summary_df.to_csv(final + '.txt', columns = ['sample','file','version'] + set_columns, index=False, sep=";")
