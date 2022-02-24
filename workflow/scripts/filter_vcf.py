#!/usr/bin/env python
# -*- coding: utf-8 -*-

#@uthor : Arnaud Guille

import sys
import os.path
import re
import pysam
import fisher
import math
import time
import itertools
import strandBias
from multiprocessing import Process, Manager
from optparse import OptionParser
from string import strip,join

# strip string
def strip(x) : return x.strip(' \t\n\r')

# create a igv link for variants report
def create_igv_link(chrom,pos,port=60151):
	link = "http://localhost:%s/goto?locus=%s:%s"%(port,chrom,pos)
	return link

# compute average of x
def average(x):
	if len(x) > 0 :
		return sum(x)*1.0/len(x)

# compute variance of x
def variance(x):
	if len(x) > 0 :
		avg = average(x)
		variance = map(lambda x: (x-avg)**2,x)
		return average(variance)

# compute Fisher's exact test for the contingency table
# ref+ ref-
# alt+ alt-
def fisher_test(a1,a2,b1,b2) :
	p = fisher.pvalue(a1,a2,b1,b2)
	p = p.two_tail
	p = round(-10*math.log10(p),3)
	return p

#get ATGC DEL/INS count for a given position
def get_bases_at_pos(chrom,pos,sam_file):
	bases_strand_plus = []
	bases_strand_minus = []
	for pileupcolumn in sam_file.pileup(chrom,pos-1,pos):
		posi = pileupcolumn.pos
		n = pileupcolumn.n
		if posi==pos-1 :
			for pileupread in pileupcolumn.pileups:
				try :
					base = pileupread.alignment.seq[pileupread.query_position]
					#qual = ord(pileupread.alignment.qual[pileupread.qpos])
					#if qual < 20 :
						#base = None
				except TypeError :
					base = None
				if pileupread.alignment.is_reverse :
					if pileupread.indel <= -1 :
						bases_strand_minus.append("-1")
					elif pileupread.indel >= 1 :
						bases_strand_minus.append("1")
					else :
						bases_strand_minus.append(base)
				else :
					if pileupread.indel <= -1 :
						bases_strand_plus.append("-1")
					elif pileupread.indel >= 1 :
						bases_strand_plus.append("1")
					else :
						bases_strand_plus.append(base)
	return [bases_strand_plus,bases_strand_minus]

# compute variance of indel size for a given position
def varx(chrom,pos,sam_file):
	tab_len_indel = []
	sum_len_indel = 0.0
	size_indel = []
	for alignedread in sam_file.fetch(chrom,pos-1,pos):
		left_pos_read = alignedread.pos
		cigar = alignedread.cigar
		pos_in_read = 0
		#print cigar,"\t",alignedread.qname,"\t",alignedread.rlen,"\t",alignedread.pos
		nb_ins = 0
		for i in range(0,len(cigar)) :
			if cigar[i][0] != 4 and cigar[i][0] != 5 :
				length = cigar[i][1]
				for j in range(pos_in_read,pos_in_read+length) :
					#print cigar[i],"\t",left_pos_read+j-nb_ins,"\t",alignedread.qname,"\t",pos
					if cigar[i][0] == 1 :
						if (left_pos_read+j) == pos + nb_ins :
							size_indel.append(cigar[i][1])
							#print left_pos_read+j,"\t",alignedread.qname,"\t",cigar[i]
					elif cigar[i][0] == 2 :
						if (left_pos_read+j-nb_ins) == pos :
							size_indel.append(-cigar[i][1])
							#print left_pos_read+j,"\t",alignedread.qname,"\t",cigar[i]
				pos_in_read = pos_in_read + length
				if cigar[i][0] == 1 :
					nb_ins = nb_ins + 1
	#print(size_indel)
	return variance(size_indel)

# get COSMIC occurrence
def nb_oc_cosmic(s):
	s = s.split("(")
	nb_cos = 0
	for elt in s :
		if re.search("(\d+)",elt):
			t = re.search("(\d+)",elt).group(1)
			nb_cos = nb_cos + int(t)
	return nb_cos

# get COSMIC overlapping the region
def cosmic_annot_indel(chrom,start,end) :
	chrom = re.sub("chr", "", chrom)
	cmd = "tabix projet_seq/annotation/hg19/cosmic/hg19_cosmic67.txt.gz "+str(chrom)+":"+str(start)+"-"+str(end)+" | grep -"
	res_cmd = os.popen(cmd).readlines()
	sum_oc = 0
	max_oc = 0
	cosmic_max_oc = ""
	cosmic_in_region = ""
	res_cmd = map(strip,res_cmd)
	for line in res_cmd :
		rec = line.split("\t")
		s = re.search("OCCURENCE=(\d+).(.*)",rec[5]).group()
		oc = nb_oc_cosmic(s)
		#print chrom,"\t",start,"\t",end,"\t",nb_oc_cosmic(s),"\t",rec[5]
		sum_oc = sum_oc + oc
		cosmic_in_region = cosmic_in_region + "|" + rec[5]
		if oc > max_oc :
			max_oc = oc
			cosmic_max_oc = rec[5]
	return [chrom,start,end,sum_oc,max_oc,cosmic_in_region,cosmic_max_oc]

# return lenght of flanking regions
def homopolymer_count(chrom,pos_vcf,fasta):
	count_left=0
	count_right=0
	base_first_left = fasta_file.fetch(chrom,pos_vcf-2,pos_vcf-1).upper()
	base_first_right = fasta_file.fetch(chrom,pos_vcf,pos_vcf+1).upper()
	i = 1
	while True :
		#go_left
		pos = pos_vcf - i
		seq_at_pos = fasta_file.fetch(chrom,pos-1,pos).upper()
		if seq_at_pos != base_first_left :
			break
		else :
			count_left = count_left + 1
			i = i + 1
	i = 1
	while True :
		#go_right
		pos = pos_vcf + i
		seq_at_pos = fasta_file.fetch(chrom,pos-1,pos).upper()
		if seq_at_pos != base_first_right :
			break
		else :
			count_right = count_right + 1
			i = i + 1

	return(count_left,count_right)


#We Reverse complement the G's and A's
def normalize_trinuc(trinuc):
	comp = {'A':'T','T':'A','C':'G','G':'C'}
	if not trinuc[1] in ['C','T']:
		return join(map(lambda x: comp[x], trinuc)[::-1], '')
	else:
		return trinuc

#complement the G's and A's
def normalize_mutcat(x):
	comp = {'A':'T','T':'A','C':'G','G':'C'}
	if not x[0] in ['C','T']:
		return comp[x[0]]+">"+comp[x[2]]
	else:
		return x

#Get Tri-nucleotide context normalized to start from C or T
def make_trinuc(chrom,pos,ref,alt):
	base_left = fasta_file.fetch(chrom,pos-2,pos-1).upper()
	base_right = fasta_file.fetch(chrom,pos,pos+1).upper()
	if base_left == "N" or base_right == "N" :
		return None
	context = base_left+ref+base_right
	context_norm = normalize_trinuc(context)
	mutcat = ref+">"+alt
	mutcat_norm = normalize_mutcat(mutcat)
	tricontext = context_norm[0] + "[" + mutcat_norm + "]" + context_norm[2]
	return tricontext


def process_line(in_queue, out_list, out_list2) :
	while True :
		##read bam file
		sam_file = pysam.Samfile(input_bam,"rb",check_sq=False,check_header=False)
		item = in_queue.get()
		line_no, line = item

		find_AO = False
		find_QA = False
		compute_QAAG = False

		dict_record = {}
		dict_homemade_stats = {}

		# exit signal
		if line == None :
			return

		#process vcf line
		if line[0] != "#" :
			record = line.split("\t")
			record = map(strip,record)
			chrom = record[0]
			pos = int(record[1])
			ref = record[3]
			list_alt = record[4].split(",")
			for alt in list_alt :
				quality = record[5]
				filt = record[6]
				infos = record[7]
				format_fields_id = record[8]
				format_record = record[9]


				str_value = ""
				str_format = ""

				#tri nucleotide context
				if len(ref) > 1 or len(alt) > 1 or ref=="*" or alt=="*" or ref=="N" or alt=="N":
					tricontext = None
				else :
					tricontext = make_trinuc(chrom,pos,ref,alt)
				#print(chrom+"\t"+str(pos)+"\t"+ref+"\t"+alt+"\t"+tricontext)

				for elt in tab_info_id :
					reg = '('+elt+'=(-?\d+\.?\d*)'+'|'+elt+'=(\w+)'+')'
					value = re.search(reg,infos).group() if re.search(reg,infos) else None
					if value != None :
						value = value.split("=")[1]

					dict_record[elt]=value
					str_value = str_value + "\t" + str(value)

				format_fields_id = format_fields_id.split(':')
				format_record = format_record.split(':')


				for elt in tab_format_id :
					try :
						index = format_fields_id.index(elt)
						format = format_record[index]
						#get AO and QA if they are present in the vcf
						if elt == "AO" :
							AO = float(format)
							find_AO = True

						if elt == "QA" :
							QA = float(format)
							find_QA = True
						if elt == "BQ" and stringency_parameter == "mutect" :
							if format == "." :
								format = 0


					except ValueError :
						format = None
					str_format = str_format + "\t" + str(format)
					dict_record[elt]=format

				#variant ID
				ID = chrom+str(pos)+ref+alt

				#compute Quality_AG if AO and QA are present in the vcf
				QAAG = QA/AO if find_AO and find_QA else None

				#count bases at each variant positions
				bases = get_bases_at_pos(chrom,pos,sam_file)

				a_p = bases[0].count("A")
				t_p = bases[0].count("T")
				g_p = bases[0].count("G")
				c_p = bases[0].count("C")
				ins_p = bases[0].count("1")
				del_p = bases[0].count("-1")

				ref_m = bases[1].count(ref)
				alt_m = bases[1].count(alt)
				a_m = bases[1].count("A")
				t_m = bases[1].count("T")
				g_m = bases[1].count("G")
				c_m = bases[1].count("C")
				ins_m = bases[1].count("1")
				del_m = bases[1].count("-1")

				#DEL
				if len(ref) > len(alt) :
					ref_p = bases[0].count(ref[0])
					alt_p = del_p
					ref_m = bases[1].count(ref[0])
					alt_m = del_m
				#INS
				elif len(ref) < len(alt) :
					ref_p = bases[0].count(ref) + ins_p
					alt_p = ins_p
					ref_m = bases[1].count(ref) + ins_m
					alt_m = ins_m
				#SNV
				else :
					ref_p = bases[0].count(ref) if stringency_parameter != "VarScan2" else int(dict_record["RDF"])
					alt_p = bases[0].count(alt) if stringency_parameter != "VarScan2" else int(dict_record["ADF"])
					ref_m = bases[1].count(ref) if stringency_parameter != "VarScan2" else int(dict_record["RDR"])
					alt_m = bases[1].count(alt) if stringency_parameter != "VarScan2" else int(dict_record["ADR"])


				if stringency_parameter == "pindel" :
					total_count = int(format.split(",")[0]) + int(format.split(",")[1])
					ref_count = int(format.split(",")[0])
					alt_count = int(format.split(",")[1])
				elif stringency_parameter == "CombinedPindelHC" :
					total_count = int(dict_record["AD"].split(",")[0]) + int(dict_record["AD"].split(",")[1])
					ref_count = int(dict_record["AD"].split(",")[0])
					alt_count = int(dict_record["AD"].split(",")[1])
				else :
					total_count = len(bases[0])+len(bases[1])
					ref_count = ref_p + ref_m
					alt_count = alt_p + alt_m

				#strand bias with fisher's exact test
				p_fisher = fisher_test(ref_p,ref_m,alt_p,alt_m)
				#strand bias score see https://github.com/tamsen/CallSomaticVariants/blob/master/CallSomaticVariants/StrandBias.cs
				strand_bias_score = strandBias.computeStrandBias(alt_p,alt_m,ref_p,ref_m,0.01)
				#homopolymer_count
				#h_count = homopolymer_count(chrom,pos,fasta_file)
				#h_count_max = max([h_count[0],h_count[1]])
				#h_count_format = fasta_file.fetch(chrom,pos-h_count[0]-1,pos-1).upper()+"|"+ref+"/"+alt+"|"+fasta_file.fetch(chrom,pos,pos+h_count[1]).upper()

				#variant frequency
				if total_count > 0 :
					if stringency_parameter == "VarScan2" or stringency_parameter == "VarScan2-FFPE" :
						variant_frequency = float(dict_record["FREQ"].strip("%"))
					elif stringency_parameter == "none" and which_vaf=="vaf_tumor":
							variant_frequency = float(dict_record["VAF_TUMOR"])*100
					elif stringency_parameter == "none" and which_vaf=="vaf":
							variant_frequency = float(dict_record["VAF"])*100
					else :
						variant_frequency = round(float(alt_count)/float(total_count),3)
				else :
					variant_frequency = 0
				#cov on each strand or not
				if ref_p > 0 and ref_m > 0 :
					cov_ref_each_strand = True
				else :
					cov_ref_each_strand = False
				if alt_p > 0 and alt_m > 0 :
					cov_alt_each_strand = True
				else :
					cov_alt_each_strand = False

				#store stats
				dict_homemade_stats["ID"]=ID
				dict_homemade_stats["chrom"]=chrom
				dict_homemade_stats["pos"]=pos
				dict_homemade_stats["ref"]=ref
				dict_homemade_stats["alt"]=alt
				dict_homemade_stats["quality"]=quality
				dict_homemade_stats["filt"]=filt
				dict_homemade_stats["a_p"]=a_p
				dict_homemade_stats["t_p"]=t_p
				dict_homemade_stats["g_p"]=g_p
				dict_homemade_stats["c_p"]=c_p
				dict_homemade_stats["a_m"]=a_m
				dict_homemade_stats["t_m"]=t_m
				dict_homemade_stats["g_m"]=g_m
				dict_homemade_stats["c_m"]=c_m
				dict_homemade_stats["ins_p"]=ins_p
				dict_homemade_stats["ins_m"]=ins_m
				dict_homemade_stats["del_p"]=del_p
				dict_homemade_stats["del_m"]=del_m
				dict_homemade_stats["ref_count"]=ref_count
				dict_homemade_stats["ref_p"]=ref_p
				dict_homemade_stats["ref_m"]=ref_m
				dict_homemade_stats["alt_count"]=alt_count
				dict_homemade_stats["alt_p"]=alt_p
				dict_homemade_stats["alt_m"]=alt_m
				dict_homemade_stats["total_count"]=total_count
				dict_homemade_stats["p_fisher"]=p_fisher
				dict_homemade_stats["SB"]=strand_bias_score
				dict_homemade_stats["variant_frequency"]=variant_frequency * 100 if stringency_parameter == "mutect2" or stringency_parameter ==  "SomaticSeq" else variant_frequency
				dict_homemade_stats["QAAG"]=QAAG
				dict_homemade_stats["ref_both_strand"]=cov_ref_each_strand
				dict_homemade_stats["alt_both_strand"]=cov_alt_each_strand
				#dict_homemade_stats["h_count_left"]=h_count[0]
				#dict_homemade_stats["h_count_right"]=h_count[1]
				#dict_homemade_stats["h_count_max"]=h_count_max
				#dict_homemade_stats["h_format"]=h_count_format
				dict_homemade_stats["tricontext"]=tricontext

			########################
			#  Variants Filtering  #
			########################

				#Default Filtering For Freebayes SNV
				if stringency_parameter == "freebayes-snp-default" :
					fields_default_filtering = ["DP","AO","p_fisher","variant_frequency","QAAG","MQM"]
					for field in fields_default_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if int(dict_record["DP"]) >= 30 and int(dict_record["AO"]) >= 5 and float(dict_homemade_stats["p_fisher"]) < 60 and float(dict_homemade_stats["variant_frequency"]) > 0.01 and float(dict_homemade_stats["QAAG"]) > 22 :
						if float(dict_homemade_stats["variant_frequency"]) < 0.20 :
							if float(dict_record["MQM"]) > 60 :
								filt = "LOW"
							else :
								filt = "ARTEFACT"
						else :
							if float(dict_homemade_stats["p_fisher"]) < 30 :
								filt = "HIGH"
							else :
								filt = "LOW"
					else :
						filt = "ARTEFACT"

				#SNV Freebayes Filtering based on benchmark with exome NA12878 and NIST
				elif stringency_parameter == "freebayes-snp-high" :
					fields_exome_filtering = ["p_fisher","variant_frequency","ref_both_strand","alt_both_strand","MQM","EPP"]
					for field in fields_exome_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) > 0.3 and float(quality) > 24 and float(dict_homemade_stats["p_fisher"]) < 30 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					elif float(dict_homemade_stats["variant_frequency"]) <= 0.3 and float(dict_record["MQM"]) > 22 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					else :
						filt = "ARTEFACT"

				#SNV Freebayes Filtering based on benchmark with exome NA12878 and NIST, add quality parameter for MiSEQ
				elif stringency_parameter == "freebayes-snp-high-miseq" :
					fields_exome_filtering = ["p_fisher","variant_frequency","ref_both_strand","alt_both_strand","MQM","EPP"]
					for field in fields_exome_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) > 0.3 and float(quality) > 24 and float(dict_homemade_stats["p_fisher"]) < 30 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					elif float(dict_homemade_stats["variant_frequency"]) <= 0.3 and float(dict_record["MQM"]) > 22 and float(dict_homemade_stats["quality"]) > 1 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					else :
						filt = "ARTEFACT"

				#SNV Freebayes Filtering Version 03 23 2015
				elif stringency_parameter == "dev" :

					if ref_p > 1 and ref_m > 1 :
						cov_ref_each_strand = True
					else :
						cov_ref_each_strand = False
					if alt_p > 1 and alt_m > 1 :
						cov_alt_each_strand = True
					else :
						cov_alt_each_strand = False

					dict_homemade_stats["ref_both_strand"]=cov_ref_each_strand
					dict_homemade_stats["alt_both_strand"]=cov_alt_each_strand

					fields_exome_filtering = ["p_fisher","variant_frequency","ref_both_strand","alt_both_strand","MQM","EPP"]

					for field in fields_exome_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) > 0.3 and float(quality) > 1 and float(dict_homemade_stats["p_fisher"]) < 30 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					elif float(dict_homemade_stats["variant_frequency"]) <= 0.3 and float(dict_record["MQM"]) > 1 and float(dict_homemade_stats["quality"]) > 1 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					elif float(dict_homemade_stats["quality"]) <= 1 :
						if dict_homemade_stats["alt_p"] >= 10 and dict_homemade_stats["alt_m"] >= 10 and float(dict_record["MQM"]) > 1  and float(dict_homemade_stats["p_fisher"]) < 30:
							filt = "PROBABLY ARTEFACT"
						else :
							filt = "ARTEFACT"
					else :
						filt = "ARTEFACT"

				#use for benchmark
				elif stringency_parameter == "bench" :

					if ref_p >= 0 and ref_m >= 0 :
						cov_ref_each_strand = True
					else :
						cov_ref_each_strand = False
					if alt_p >= 0 and alt_m >= 0 :
						cov_alt_each_strand = True
					else :
						cov_alt_each_strand = False

					dict_homemade_stats["ref_both_strand"]=cov_ref_each_strand
					dict_homemade_stats["alt_both_strand"]=cov_alt_each_strand

					fields_exome_filtering = ["p_fisher","variant_frequency","ref_both_strand","alt_both_strand","MQM","EPP"]

					for field in fields_exome_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) > 0.3 and float(quality) > 1 and float(dict_homemade_stats["p_fisher"]) < 60 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					elif float(dict_homemade_stats["variant_frequency"]) <= 0.3 and float(dict_record["MQM"]) > 1 and float(dict_homemade_stats["quality"]) > 1 :
						if dict_homemade_stats["ref_both_strand"] :
							if dict_homemade_stats["alt_both_strand"] :
								filt = "HIGH"
							else :
								filt = "ARTEFACT"
						else :
							filt = "LOW"
					elif float(dict_homemade_stats["quality"]) <= 1 :
						if dict_homemade_stats["alt_p"] >= 10 and dict_homemade_stats["alt_m"] >= 10 and float(dict_record["MQM"]) > 1  and float(dict_homemade_stats["p_fisher"]) < 60:
							filt = "PROBABLY ARTEFACT"
						else :
							filt = "ARTEFACT"
					else :
						filt = "ARTEFACT"


				#use for VarScan2
				elif stringency_parameter == "VarScan2" or stringency_parameter == "VarScan2-FFPE":

					cov_p = alt_p + ref_p
					cov_m = alt_m + ref_m

					fields_filtering = ["SB","FREQ","ABQ"]

					filt=""

					quality = str(-10*math.log10(max(float(dict_record["PVAL"]),2.2250738585072014e-308)))
					record[5] = quality

					for field in fields_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) < 5 :
						filt = filt + "LowVariantFreq"
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "|SB"
						if int(dict_record["ABQ"]) < 30 :
							filt = filt + "|LowABQ"
					else :
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "SB"
							else :
								filt = "PASS"
						else :
							filt = "PASS"

					if stringency_parameter == "VarScan2-FFPE" :
						tri = dict_homemade_stats["tricontext"]
						if tri == "A[C>T]G" or tri == "C[C>T]G" or tri == "G[C>T]G" or tri == "T[C>T]G" :
							filt = filt + "|FFPE"

				#use for VarScan2-somatic
				elif stringency_parameter == "VarScan2-somatic" :

					cov_p = alt_p + ref_p
					cov_m = alt_m + ref_m

					fields_filtering = ["SB","FREQ","SS","SSC"]

					filt="."

					quality = str(float(dict_record["SSC"]))
					record[5] = quality

					for field in fields_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if dict_record["SS"] == '2' or dict_record["SS"] == '3' :
						filt = "PASS"
					elif dict_record["SS"] == '1' :
						filt = "Germline"
					else :
						filt = "."

					if float(dict_record["SSC"]) < 20 :
						filt = filt + "|LowSSC"


					if float(dict_homemade_stats["variant_frequency"]) < 0.05 :
						filt = filt + "|LowVariantFreq"
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "|SB"
						#if int(dict_record["ABQ"]) < 30 :
						#	filt = filt + "|LowABQ"
					else :
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "|SB"


				#use for Mutect
				elif stringency_parameter == "mutect" :

					cov_p = alt_p + ref_p
					cov_m = alt_m + ref_m

					fields_filtering = ["SB","variant_frequency","BQ"]

					filt=""

					#print chrom,pos,ref

					for field in fields_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) < 0.05 :
						filt = filt + "LowVariantFreq"
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "|SB"
						if int(dict_record["BQ"]) < 30 :
							filt = filt + "|LowABQ"
					else :
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "SB"
							else :
								filt = "PASS"
						else :
							filt = "PASS"

				#use for Mutect2
				elif stringency_parameter == "mutect2" :

					cov_p = alt_p + ref_p
					cov_m = alt_m + ref_m

					fields_filtering = ["SB","variant_frequency","AD","TLOD","QSS"]

					filt=record[6]

					qss = dict_record["QSS"].split(",")
					ad = dict_record["AD"].split(",")

					try :
						if float(ad[1]) != 0 :
							abq = float(qss[1])/float(ad[1])
						else :
							abq = 0
					except IndexError :
						abq = 0

					#print chrom,pos,ref
					for field in fields_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) < 5 :
						filt = filt + "|LowVariantFreq"
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "|SB"
						if abq < 30 :
							filt = filt + "|LowABQ"
						if float(dict_record["TLOD"]) < 10 :
							filt = filt + "|LowTLOD"
					else :
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["SB"] > -3.01029995664 :
								filt = filt + "|SB"
							if abq < 30 :
								filt = filt + "|LowABQ"
							if float(dict_record["TLOD"]) < 10 :
								filt = filt + "|LowTLOD"
						else :
							if abq < 30 :
								filt = filt + "|LowABQ"
							if float(dict_record["TLOD"]) < 10 :
								filt = filt + "|LowTLOD"




				#SNV Varscan Filtering based on benchmark with exome NA12878 and NIST
				elif stringency_parameter == "varscan-snp-high" :
					fields_exome_filtering = ["p_fisher","variant_frequency","ref_both_strand","alt_both_strand","ABQ","GQ","alt_count","DP"]
					for field in fields_exome_filtering :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_record["GQ"]) > 12 :
						if float(dict_homemade_stats["variant_frequency"]) > 0.2 :
							if float(dict_record["ABQ"]) > 30 :
								if float(dict_homemade_stats["p_fisher"]) < 30 :
									filt="HIGH"
								else :
									filt="ARTEFACT"
							else :
								filt="ARTEFACT"
						else :
							if float(dict_homemade_stats["alt_count"]) > 8 :
								if float(dict_record["ABQ"]) > 30 :
									if dict_homemade_stats["ref_both_strand"] :
										if dict_homemade_stats["alt_both_strand"] :
											filt = "LOW"
										else :
											filt = "ARTEFACT"
									else :
										filt = "LOW"
								else :
									filt = "ARTEFACT"
							else :
								filt = "ARTEFACT"



					else:
						if float(dict_homemade_stats["variant_frequency"]) > 0.3  and float(dict_record["DP"]) > 8:
							if float(dict_record["GQ"]) > 6 :
								if float(dict_record["ABQ"]) > 30 :
									filt="LOW"
								else :
									filt="ARTEFACT"
							else :
								if dict_homemade_stats["ref_both_strand"] :
									if dict_homemade_stats["alt_both_strand"] :
										filt = "HIGH"
									else :
										filt = "ARTEFACT"
								else :
									filt = "LOW"
						else :
							filt = "ARTEFACT"


				#INDEL Filtering PINDEL
				elif stringency_parameter == "pindel" :
					if float(dict_homemade_stats["variant_frequency"]) >= 0.05 and int(dict_homemade_stats["alt_count"]) >= 10 :
						filt="PASS"
					else :
						filt="."

				elif stringency_parameter == "CombinedPindelHC" :
					if float(dict_homemade_stats["variant_frequency"] >= 0.05 and int(dict_homemade_stats["alt_count"]) >= 10) :
						filt="PASS"
					else :
						filt="."

				elif stringency_parameter == "IndelLocator" :
					if re.search('SOMATIC',line) :
						filt = "SOMATIC"
					else :
						filt = "."

				elif stringency_parameter == "scalpel":
					filt = record[6]

				elif stringency_parameter == "none" :
					filt = record[6]



				#INDEL Filtering UnifiedGenotyper
				elif stringency_parameter == "unifiedGenotyper-indel" :
					fields_unifiedGenotyper = ["variant_frequency","ins_p","ins_m","del_p","del_m"]
					for field in fields_unifiedGenotyper :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if len(alt) > len(ref) :
						if float(dict_homemade_stats["variant_frequency"]) > 0.2 :
							if int(dict_homemade_stats["ins_p"]) > 1 and int(dict_homemade_stats["ins_m"]) > 1 :
								filt = "LOW"
							else :
								filt = "ARTEFACT"
						else :
							filt="ARTEFACT"

					else :
						if float(dict_homemade_stats["variant_frequency"]) > 0.8 :
							if float(quality) > 1000 :
								if int(dict_homemade_stats["del_p"]) > 5 and int(dict_homemade_stats["del_m"]) > 5 :
									filt = "LOW"
								else :
									filt = "ARTEFACT"
							else :
								filt = "ARTEFACT"
						else :
							filt = "ARTEFACT"


				#SomaticSeq
				elif stringency_parameter == "SomaticSeq" or stringency_parameter == "SomaticSeq-FFPE" :
					fields_somaticSeq = ["p_fisher","variant_frequency"]

					cov_p = alt_p + ref_p
					cov_m = alt_m + ref_m

					for field in fields_somaticSeq :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_homemade_stats["variant_frequency"]) < 5 :
						filt = filt + "|LowVariantFreq"
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["p_fisher"] >= 60 :
								filt = filt + "|SB"
					else :
						if cov_p > 0 and cov_m > 0 :
							if dict_homemade_stats["p_fisher"] >= 60 :
								filt = filt + "|SB"

					if stringency_parameter == "SomaticSeq-FFPE" :
						tri = dict_homemade_stats["tricontext"]
						if tri == "A[C>T]G" or tri == "C[C>T]G" or tri == "G[C>T]G" or tri == "T[C>T]G" :
							filt = filt + "|FFPE"

				#Strelka
				elif stringency_parameter == "Strelka" :
					fields_Strelka = ["SomaticEVS"]


					for field in fields_Strelka :
						if field not in dict_record.keys() and field not in dict_homemade_stats.keys() :
							print "Error : Could not find field %s for filtering"%field
							sys.exit(1)

					if float(dict_record["SomaticEVS"]) >= 4 and filt == "LowEVS" :
						filt="PASS"

					quality = dict_record["SomaticEVS"]

				record[5] = quality
				record[6] = filt


				line_format = ID + "\t" + sample_name + "\t" + chrom + "\t" + str(pos) + "\t" + ref + "\t" + alt + "\t" + str(quality) + "\t" + str(filt) + "\t" + str(dict_homemade_stats["tricontext"]) + "\t" + str(a_p) + "\t" + str(t_p) + "\t" + str(g_p) + "\t" + str(c_p) + "\t" + str(a_m) + "\t" + str(t_m) + "\t" + str(g_m) + "\t" + str(c_m) + "\t" + str(ins_p) + "\t" + str(ins_m) + "\t" + str(del_p) + "\t" + str(del_m) + "\t" + str(ref_count) + "\t" + str(ref_p) + "\t" + str(ref_m) + "\t" + str(alt_count) + "\t" + str(alt_p) + "\t" + str(alt_m) + "\t" + str(total_count) + "\t" + str(p_fisher) + "\t" + str(strand_bias_score) + "\t" + str(dict_homemade_stats["variant_frequency"]) + "\t" + str(QAAG) + str_value + str_format

				res_stats = (line_no,line_format)
				out_list.append(res_stats)

			vcf_line = (line_no,record)
			out_list2.append(vcf_line)

		else :
			record = line.split("\t")
			record = map(strip,record)
			vcf_line = (line_no,record)
			out_list2.append(vcf_line)

##################################
#~~~~~~~~~~~~~~MAIN~~~~~~~~~~~~~~#
##################################

if __name__ == "__main__" :

	usage = "usage: %prog [options] <input.vcf> <input.bam> <ref.fa> <output>"
	desc = "Filter vcf file and output two files. The filtered vcf file with FILTER field updated and another file with all the statistics contains in the original vcf file plus some home-made statistics"
	ver = "1.0"
	parser = OptionParser(usage=usage,description=desc,version=ver)
	parser.add_option("-p", help="number of thread [default: %default]", dest="thread", default=1)
	parser.add_option("-f","--filter", help="filtering mode [default: %default]", default="default", dest="mode")

	(options, args) = parser.parse_args()

	nb_args = len(sys.argv)

	if nb_args < 4 :
		print "bad usage type --help"
		sys.exit(1)

	#mandatory args
	input_vcf = sys.argv[nb_args-4]
	input_bam = sys.argv[nb_args-3]
	input_fasta = sys.argv[nb_args-2]
	output_path = sys.argv[nb_args-1]

	try :
		num_workers = int(options.thread)
	except ValueError :
		print "-p must be an integer"
		sys.exit(1)

	stringency_parameter = options.mode
	if stringency_parameter not in ["freebayes-snp-default","freebayes-snp-high","freebayes-snp-high-miseq","varscan-snp-high","freebayes-indel","unifiedGenotyper-snp","unifiedGenotyper-indel","pindel","dev","VarScan2","bench","CombinedPindelHC","IndelLocator","scalpel","mutect","mutect2","VarScan2-somatic","SomaticSeq","SomaticSeq-FFPE","VarScan2-FFPE","Strelka","none"] :
		print "-f : accepted values are freebayes-snp-default, freebayes-snp-high,freebayes-snp-high-miseq, freebayes-indel, unifiedGenotyper-snp, unifiedGenotyper-indel, IndelLocator, mutect, scalpel, none"
		sys.exit(1)

	#output stats files
	output_stats = open(output_path+".stats","wb")
	output_vcf_filtered = open(output_path+".filtered.vcf","wb")

	#read vcf file
	f = open(input_vcf,"rb")
	lines = f.readlines()

	#read bam file
	sam_file = pysam.Samfile(input_bam,"rb",check_sq=False,check_header=False)

	#read fasta file
	fasta_file = pysam.Fastafile(input_fasta)

	#format fields ID
	tab_format_id = []

	#info fields ID
	tab_info_id = []

	#sample name
	sample_name = os.path.basename(input_vcf).split(".")[0]

	#parse vcf file
	for line in lines :
		if line[0] == "#" :
			if re.search('[#]+INFO=<ID=([a-zA-Z0-9_-]+),.*',line) :
				info_id = re.search('[#]+INFO=<ID=([a-zA-Z0-9_-]+),.*',line).group(1)
				tab_info_id.append(info_id)

			elif re.search('[#]+FORMAT=<ID=([a-zA-Z0-9_-]+),.*',line) :
				format_id = re.search('[#]+FORMAT=<ID=([a-zA-Z0-9_-]+),.*',line).group(1)
				tab_format_id.append(format_id)
		else :
			break

	if "VAF" in tab_info_id :
		which_vaf="vaf"
	elif "VAF_TUMOR" in tab_info_id :
		which_vaf="vaf_tumor"
	else :
		which_vaf="pysam"

	str_col1 = ""
	str_col2 = ""

	for elt in tab_info_id :
		str_col1 = str_col1 + "\t" + elt
	for elt in tab_format_id :
		str_col2 = str_col2 + "\t" + elt

	col = "ID\tSAMPLE\tCHROM\tPOS\tREF\tALT\tQUALITY\tFILTER\tTRI_CONTEXT\tA+\tT+\tG+\tC+\tA-\tT-\tG-\tC-\tINS+\tINS-\tDEL+\tDEL-\tREF_COUNT\tREF+\tREF-\tALT_COUNT\tALT+\tALT-\tTOTAL_COUNT\tSTRAND_BIAS\tSB\tVARIANT_FREQUENCY\tQUALITY_AG" + str_col1 + str_col2

	try :
		with open(input_vcf) as f :
			#multithreading to parse vcf file
			manager = Manager()
			results = manager.list()
			results2 = manager.list()
			work = manager.Queue(num_workers)

			# start for workers
			pool = []
			for i in xrange(num_workers):
				p = Process(target=process_line, args=(work,results,results2))
				p.start()
				pool.append(p)

			iters = itertools.chain(f, (None, )*num_workers)

			for num_and_line in enumerate(iters):
				work.put(num_and_line)

		for p in pool:
			p.join()

		output_stats.write("%s"%col)
		output_stats.write("\n")

		for line in sorted(results) :
			output_stats.write("%s"%line[1])
			output_stats.write("\n")

		for line in sorted(results2) :
			output_vcf_filtered.write("\t".join(line[1]))
			output_vcf_filtered.write("\n")


		output_stats.close()
		output_vcf_filtered.close()

	except IOError :
   			print "Error !!! Cant open vcf %s"%input_vcf
   			sys.exit(1)
