#!/usr/bin/env python
# -*- coding: utf-8 -*-
#bedtools coverage -sorted -d -g /mnt/homedata/ngs_om/travail/annot/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/chrLength.txt -a /mnt/homedata/ngs_om/travail/annot/Homo_sapiens/personal/hg19/Annotation/target_files/Haloplex_Illumina_V10_CCP.bed -b projet_seq/projet_CTC_EM2/01001-P1-M1_CTC_COLON_N_01001/NormalSample/CTC_COLON_N_01001_S3.sort.RG.realigned.recall.bam | sort -k4,4 | bedtools groupby -g 4 -c 6 -o sum > Desktop/All_LakL_mapping/macs2_output/sl4-i27/sl3_i27_peaks.sum.cov


import sys
import os.path
import pandas as pd

# strip string
def strip(x) : return x.strip(' \t\n\r')

##################################
#~~~~~~~~~~~~~~MAIN~~~~~~~~~~~~~~#
##################################

if __name__ == "__main__" :

	input_cov = sys.argv[1]
	out_file = sys.argv[2]

	df = pd.read_csv(input_cov,sep="\t",names=["chrom","start","end","symbol","pos","nb_reads"])
	df_by_symbol = df.groupby("symbol")
	mean_by_symbol = df_by_symbol.nb_reads.mean()
	std_by_symbol = df_by_symbol.nb_reads.std()
	uniformity = std_by_symbol/mean_by_symbol
	cov_at_1 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 1).sum()/ts.count())
	cov_at_5 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 5).sum()/ts.count())
	cov_at_10 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 10).sum()/ts.count())
	cov_at_20 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 20).sum()/ts.count())
	cov_at_50 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 50).sum()/ts.count())
	cov_at_100 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 100).sum()/ts.count())
	cov_at_150 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 150).sum()/ts.count())
	cov_at_200 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 200).sum()/ts.count())
	cov_at_300 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 300).sum()/ts.count())
	cov_at_400 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 400).sum()/ts.count())
	cov_at_500 = df_by_symbol.nb_reads.agg(pos=lambda ts: (ts >= 500).sum()/ts.count())
	out = pd.concat([mean_by_symbol,std_by_symbol,uniformity,cov_at_1,cov_at_5,cov_at_10,cov_at_20,cov_at_50,cov_at_100,cov_at_150,cov_at_200,cov_at_300,cov_at_400,cov_at_500], axis=1)
	out.to_csv(out_file,index=True,header=False,sep="\t")
