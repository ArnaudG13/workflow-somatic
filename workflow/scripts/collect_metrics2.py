#!/usr/bin/env python
# -*- coding: utf-8 -*-

#@uthor : Arnaud Guille
#This script parses qualimapReport html files and others to collect metrics

import re
import sys

input_file_bam_qc = sys.argv[1]
input_file_cov = sys.argv[2]
input_file_bait_bias = sys.argv[3]
input_file_pre_adapt_bias = sys.argv[4]
input_file_error_summary = sys.argv[5]
#input_file_cont_table = sys.argv[6]
input_file_fastqc1 = sys.argv[6]
input_file_fastqc2 = sys.argv[7]
output = sys.argv[8]

file_in_bamqc = open(input_file_bam_qc,"r")
file_in_cov = open(input_file_cov,"r")
file_in_bait_bias = open(input_file_bait_bias,"r")
file_in_pre_adapt_bias = open(input_file_pre_adapt_bias,"r")
file_in_error_summary = open(input_file_error_summary,"r")
#file_in_cont_table = open(input_file_cont_table,"r")
file_in_fastqc1 = open(input_file_fastqc1,"r")
file_in_fastqc2 = open(input_file_fastqc2,"r")
output_file = open(output,"w")
file = file_in_bamqc.read()

name = re.search("<td class=column1>(BAM file: )</td>\s<td class=column2>(.*)</td>",file).group(2)
nb_reads = re.search("<td class=column1>(Number of reads)</td>\s<td class=column2>(.*)</td>",file).group(2)
mapped_reads = re.findall("<td class=column1>(Mapped reads)</td>\s<td class=column2>(.*)</td>",file)
mapped_reads_global = mapped_reads[0][1]
mapped_reads_in_region = mapped_reads[1][1]
unmapped_reads = re.search("<td class=column1>(Unmapped reads)</td>\s<td class=column2>(.*)</td>",file).group(2)
read_stats = re.search("<td class=column1>(Read min/max/mean length)</td>\s<td class=column2>(.*)</td>",file).group(2)
region_size = re.search("<td class=column1>(Regions size/percentage of reference)</td>\s<td class=column2>(.*)</td>",file).group(2)
coverage = re.search("<td class=column1>(Mean)</td>\s<td class=column2>(.*)</td>",file).group(2)
mapping_quality = re.search("<td class=column1>(Mean Mapping Quality)</td>\s<td class=column2>(.*)</td>",file).group(2)
error_rate = re.search("<td class=column1>(General error rate)</td>\s<td class=column2>(.*)</td>",file).group(2)
insert_size = re.search("<td class=column1>P25/Median/P75</td>\s<td class=column2>(.*)</td>",file).group(1)
dup_reads = re.search("<td class=column1>Duplicat.*</td>\s<td class=column2>(.*)</td>",file).group(1)

lines_fastqc1 = file_in_fastqc1.read()
nb_reads_tot1 = re.search("<td>Total Sequences</td><td>(\d+).*",lines_fastqc1).group(1)
lines_fastqc2 = file_in_fastqc2.read()
nb_reads_tot2 = re.search("<td>Total Sequences</td><td>(\d+).*",lines_fastqc2).group(1)
nb_reads_tot = float(nb_reads_tot1) + float(nb_reads_tot2)
nb_reads = nb_reads.replace(",","")
nb_reads = nb_reads.replace(" ","")

percent_nb_reads_keep = round(((float(nb_reads)/nb_reads_tot)*100),2)
format_nb_reads = str(nb_reads) + "/" + str(nb_reads_tot) + " (" + str(percent_nb_reads_keep) + "%)"

#compute mean of coverage uniformity

cov_uni = 0
cpt = 0

lines = file_in_cov.readlines()

for line in lines :
	record = line.split("\t")
	val = record[3]
	if val != "NA" and val != "Nan" and val != "":
		cpt = cpt + 1
		cov_uni += float(record[3])

cov_uni = round(cov_uni/cpt,3)


#bait bias
lines_bait_bias = file_in_bait_bias.readlines()
Cref = lines_bait_bias[10].split("\t")[4]
Gref = lines_bait_bias[15].split("\t")[4]

#pre adaptater bias
lines_adapt_bias = file_in_pre_adapt_bias.readlines()
Deamination = lines_bait_bias[12].split("\t")[4]
OxoG = lines_bait_bias[15].split("\t")[4]

#Error summary
lines_error_sum = file_in_error_summary.readlines()
ac = lines_error_sum[7].split("\t")[5].strip(" \t\n\r")
ag = lines_error_sum[8].split("\t")[5].strip(" \t\n\r")
at = lines_error_sum[9].split("\t")[5].strip(" \t\n\r")
ca = lines_error_sum[10].split("\t")[5].strip(" \t\n\r")
cg = lines_error_sum[11].split("\t")[5].strip(" \t\n\r")
ct = lines_error_sum[12].split("\t")[5].strip(" \t\n\r")

#Cont table
#lines_cont = file_in_cont_table.readlines()
#cont = lines_cont[1].split("\t")[1]
#cont_error = lines_cont[1].split("\t")[2].strip(" \t\n\r")

output_file.write(name+"\t"+format_nb_reads+"\t"+mapped_reads_global+"\t"+unmapped_reads+"\t"+mapped_reads_in_region+"\t"+read_stats+"\t"+region_size+"\t"+coverage+"\t"+mapping_quality+"\t"+error_rate+"\t"+insert_size+"\t"+dup_reads+"\t"+str(cov_uni)+"\t"+str(Deamination)+"\t"+str(OxoG)+"\t"+str(Cref)+"\t"+str(Gref)+"\t"+str(ac)+"\t"+str(ag)+"\t"+str(at)+"\t"+str(ca)+"\t"+str(cg)+"\t"+str(ct)+"\n")

file_in_bamqc.close()
file_in_cov.close()
output_file.close()
