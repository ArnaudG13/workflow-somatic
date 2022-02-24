#!/usr/bin/env python
# -*- coding: utf-8 -*-

#@uthor : Arnaud Guille

import sys
import os.path
import csv
import re
import math
from optparse import OptionParser

csv.field_size_limit(sys.maxsize)

# strip string
def strip(x) : return x.strip(' \t\n\r')

#get COSMIC occurrence
def nb_oc_cosmic(s):
	s = s.split("(")
	nb_cos = 0
	for elt in s :
		if re.search("(\d+)",elt):
			t = re.search("(\d+)",elt).group(1)
			nb_cos = nb_cos + int(t)
	return nb_cos

#get consensus vote
def vote(row):

	score = 0
	vote = 0

	if row[sift_key] != "." and  row[sift_key] != "" and row[sift_key] != " " :
		vote = vote + 1
	if row[sift_key] == "D" :
		score = score + 1

	if row[polyphen_key] != "." and row[polyphen_key] != "" and  row[polyphen_key] != " " :
		vote = vote + 1
	if row[polyphen_key] == "D" :
		score = score + 1

	if lrt_found :
		if row[lrt_key] != "." and row[lrt_key] != "" and row[lrt_key] != " " :
			vote = vote + 1
		if row[lrt_key] == "D" :
			score = score + 1

	if mut_taster_found :
		if row[mut_taster_key] != "." and row[mut_taster_key] != "" and row[mut_taster_key] != " " :
			vote = vote + 1
		if row[mut_taster_key] == "D" or row[mut_taster_key] == "A" :
			score = score + 1

	if mut_asses_found :
		if row[mut_asses_key] != "." and row[mut_asses_key] != "" and row[mut_asses_key] != " " :
			vote = vote + 1
		if row[mut_asses_key] == "M" or row[mut_asses_key] == "H" :
			score = score + 1

	if fat_hmm_found :
		if row[fat_hmm_key] != "." and row[fat_hmm_key] != "" and row[fat_hmm_key] != " " :
			vote = vote + 1
		if row[fat_hmm_key] == "D" :
			score = score + 1

	if provean_found :
		if row[provean_key] != "." and row[provean_key] != "" and row[provean_key] != " " :
			vote = vote + 1
		if row[provean_key] == "D" :
			score = score + 1

	if metasvm_found :
		if row[metasvm_key] != "." and row[metasvm_key] != "" and row[metasvm_key] != " " :
			vote = vote + 1
		if row[metasvm_key] == "D" :
			score = score + 1

	if metaLR_found :
		if row[metaLR_key] != "." and row[metaLR_key] != "" and row[metaLR_key] != " " :
			vote = vote + 1
		if row[metaLR_key] == "D" :
			score = score + 1

	if mcap_found :
		if row[mcap_key] != "." and row[mcap_key] != "" and row[mcap_key] != " " :
			vote = vote + 1
		if row[metasvm_key] == "D" :
			score = score + 1

	if fatMKL_found :
		if row[fatMKL_key] != "." and row[fatMKL_key] != "" and row[fatMKL_key] != " " :
			vote = vote + 1
		if row[metasvm_key] == "D" :
			score = score + 1

	return (score,vote)


##################################
#~~~~~~~~~~~~~~MAIN~~~~~~~~~~~~~~#
##################################

if __name__ == "__main__" :

	usage = "usage: %prog [options] <input.vcf> <output>"
	desc = "Prioritize variants according to COSMIC, prediction tools, clinvar, dbsnp"
	ver = "1.0"
	parser = OptionParser(usage=usage,description=desc,version=ver)
	parser.add_option("-m", help="mode [default: %default]", dest="mode", default="snp")

	(options, args) = parser.parse_args()

	nb_args = len(sys.argv)

	if nb_args < 2 :
		print "bad usage type --help"
		sys.exit(1)

	mode = options.mode
	if mode not in ["snp","indel"] :
		print "-m : accepted values are snp, indel"
		sys.exit(1)

	#mandatory args
	input_variants = sys.argv[nb_args-2]
	output_path = sys.argv[nb_args-1]

	output = open(output_path,"wb")
	variants = open(input_variants,"rb")
	csv_variants = csv.DictReader(variants,delimiter="\t")

	dict_mut = {}
	csv_variants = csv.DictReader(variants,delimiter="\t",skipinitialspace=True)

	header = csv_variants.fieldnames

	#mandatory
	cosmic_key = [x for x in header if re.search(r'^cosmic[0-9]+[_coding]?', x)]
	try :
		cosmic_key = header[header.index(cosmic_key[0])]
	except IndexError :
		print("COSMIC FIELD MISSING")
		sys.exit(0)

	#mandatory 1
	sift_key = [x for x in header if re.search(r'SIFT_pred', x)]
	try :
		sift_key = header[header.index(sift_key[0])]
	except IndexError :
		print("SIFT FIELD MISSING")
		sys.exit(0)

	#optional 2
	polyphen_key = [x for x in header if re.search(r'Polyphen2_HDIV_pred', x)]
	try :
		polyphen_key = header[header.index(polyphen_key[0])]
		polyphen_found=True
	except IndexError :
		polyphen_found=False
		print("Polyphen2 FIELD MISSING")

	#optional 3
	lrt_key = [x for x in header if re.search(r'LRT_pred', x)]
	try :
		lrt_key = header[header.index(lrt_key[0])]
		lrt_found=True
	except IndexError :
		lrt_found=False
		print("LRT FIELD MISSING")

	#opional 4
	mut_taster_key = [x for x in header if re.search(r'MutationTaster_pred', x)]
	try :
		mut_taster_key = header[header.index(mut_taster_key[0])]
		mut_taster_found=True
	except IndexError :
		mut_taster_found=False
		print("MutationTaster FIELD MISSING")

	#optional 5
	mut_asses_key = [x for x in header if re.search(r'MutationAssessor_pred', x)]
	try :
		mut_asses_key = header[header.index(mut_asses_key[0])]
		mut_asses_found=True
	except IndexError :
		mut_asses_found=False
		print("MutationAssessor FIELD MISSING")

	#opional 6
	fat_hmm_key = [x for x in header if re.search(r'FATHMM_pred', x)]
	try :
		fat_hmm_key = header[header.index(fat_hmm_key[0])]
		fat_hmm_found=True
	except IndexError :
		fat_hmm_found=False
		print("FATHMM FIELD MISSING")

	#opional 7
	provean_key = [x for x in header if re.search(r'PROVEAN_pred', x)]
	try :
		provean_key = header[header.index(provean_key[0])]
		provean_found = True
	except IndexError :
		provean_found = False
		print("PROVEAN FIELD MISSING")

	#optional 8
	metasvm_key = [x for x in header if re.search(r'MetaSVM_pred', x)]
	try :
		metasvm_key = header[header.index(metasvm_key[0])]
		metasvm_found = True
	except IndexError :
		metasvm_found = False
		print("METASVM FIELD MISSING")

	#optional 9
	metaLR_key = [x for x in header if re.search(r'MetaLR_pred', x)]
	try :
		metaLR_key = header[header.index(metaLR_key[0])]
		metaLR_found = True
	except IndexError :
		metaLR_found = False
		print("METALR FIELD MISSING")

	#optional 10
	mcap_key = [x for x in header if re.search(r'M-CAP_pred', x)]
	try :
		mcap_key = header[header.index(mcap_key[0])]
		mcap_found = True
	except IndexError :
		mcap_found = False
		print("M-CAP FIELD MISSING")

	#optional 11
	fatMKL_key = [x for x in header if re.search(r'fathmm-MKL_coding_pred', x)]
	try :
		fatMKL_key = header[header.index(fatMKL_key[0])]
		fatMKL_found = True
	except IndexError :
		fatMKL_found = False
		print("fathmm-MKL FIELD MISSING")

	#optional 12
	metaRNN_key = [x for x in header if re.search(r'MetaRNN_pred', x)]
	try :
		metaRNN_key = header[header.index(metaRNN_key[0])]
		metaRNN_found = True
	except IndexError :
		metaRNN_found = False
		print("MetaRNN FIELD MISSING")

	################################
	#pop allele frequency database #
	################################

	#exac non TCGA ALL
	exac_key = [x for x in header if re.search(r'ExAC_nontcga_ALL', x)]
	try :
		exac_key = header[header.index(exac_key[0])]
		exac_found = True
	except IndexError :
		exac_found = False
		print("Exac_nontcga_ALL FIELD MISSING")

	#Kaviar AF
	kaviar_key = [x for x in header if re.search(r'Kaviar_AF', x)]
	try :
		kaviar_key = header[header.index(kaviar_key[0])]
		kaviar_found = True
	except IndexError :
		kaviar_found = False
		print("Kaviar_AF FIELD MISSING")

	#HRC AF
	hrc_key = [x for x in header if re.search(r'HRC_AF', x)]
	try :
		hrc_key = header[header.index(hrc_key[0])]
		hrc_found = True
	except IndexError :
		hrc_found = False
		print("HRC_AF FIELD MISSING")

	#ESP 6500
	esp_key = [x for x in header if re.search(r'esp6500siv2_all', x)]
	try :
		esp_key = header[header.index(esp_key[0])]
		esp_found = True
	except IndexError :
		esp_found = False
		print("esp6500siv2_all FIELD MISSING")

	#1000g mandatory
	g1000_key = [x for x in header if re.search(r'1000g', x)]
	try :
		g1000_key = header[header.index(g1000_key[0])]
	except IndexError :
		print("1000G FIELD MISSING")
		sys.exit(0)

	#avsnp
	avsnp_key = [x for x in header if re.search(r'avsnp150', x)]
	try :
		avsnp_key = header[header.index(avsnp_key[0])]
		avsnp_found = True
	except IndexError :
		avsnp_found = False
		print("avsnp150 FIELD MISSING")

	#avsnp AF
	af_key = [x for x in header if re.search(r'^AF$', x)]
	try :
		af_key = header[header.index(af_key[0])]
		af_found = True
	except IndexError :
		af_found = False
		print("AF FIELD MISSING")

	#snp129
	snp129_key = [x for x in header if re.search(r'snp129', x)]
	try :
		snp129_key = header[header.index(snp129_key[0])]
		snp129_found = True
	except IndexError :
		snp129_found = False
		print("snp129 FIELD MISSING")

	#snp138NonFlagged
	snp138_key = [x for x in header if re.search(r'snp138NonFlagged', x)]
	try :
		snp138_key = header[header.index(snp138_key[0])]
		snp138_found = True
	except IndexError :
		snp138_found = False
		print("snp138NonFlagged FIELD MISSING")


	for record in csv_variants :
		try :
			if record["ID"] in dict_mut.keys():
				dict_mut[record["ID"]]=dict_mut[record["ID"]] + 1
			else :
				dict_mut[record["ID"]]=1
		except KeyError :
			chrid = record["Chr"] + record["position"] + record["ref_allele"] + record["alt_allele"]
			if chrid in dict_mut.keys():
				dict_mut[chrid]=dict_mut[chrid] + 1
			else :
				dict_mut[chrid]=1

	variants.close

	variants = open(input_variants,"rb")
	csv_variants = csv.DictReader(variants,delimiter="\t")
	try :
		nb_sample = len(set(record["SAMPLE"] for record in csv_variants))
	except KeyError :
		nb_sample = len(set(record["tumor_name"] for record in csv_variants))
	variants.close()

	variants = open(input_variants,"rb")
	csv_variants = csv.DictReader(variants,delimiter="\t")
	fieldnames = csv_variants.fieldnames

	if mode == "snp" or mode == "indel" :

		for field in fieldnames :
			output.write("%s\t"%field)
		output.write("COMMENT\tVOTE_DAMAGING\tNB_VOTE\tEFFECT\tFREQ")
		output.write("\n")

		for record in csv_variants :

			func=""
			v=""
			comment = ""
			effect = ""

			#OCCURRENCE IN DATASET
			try :
				freq_in_dataset = dict_mut[record["ID"]]
			except KeyError :
				chrid = record["Chr"] + record["position"] + record["ref_allele"] + record["alt_allele"]
				freq_in_dataset = dict_mut[chrid]

			#COSMIC OCCURRENCE
			if record[cosmic_key] is not None :
				if  re.search('COS[MV]',record[cosmic_key]) :
					s = re.search("OCCURENCE=(\d+).(.*)",record[cosmic_key]).group()
					nb_cos=nb_oc_cosmic(s)
				else :
					nb_cos=0
			else :
				nb_cos=0

			#rs dbSNP129
			if snp129_found is True :
				if record[snp129_key] is not None :
					if re.search('rs',record[snp129_key]):
						rs129 = record[snp129_key]
					else :
						rs129 = None
				else :
					rs129 = None
			else :
				rs129 = None

			#rs dbsnp138NonFlagged
			if snp138_found is True :
				if record[snp138_key] is not None :
					if re.search('rs',record[snp138_key]):
						rs138 = record[snp138_key]
					else :
						rs138 = None
				else : 
					rs138 = None
			else :
				rs138 = None

			if avsnp_found is True :
				if record[avsnp_key] is not None :
					if re.search('rs',record[avsnp_key]):
						rs = record[avsnp_key]
					else :
						rs = None
				else :
					rs = None
			else :
				rs = None

			#If the variants are very frequent in our dataset and absent from public database, it's a new hotspot, but in fact an artefact
			if nb_cos < 10 and rs129 is None and rs138 is None and rs is None :
				if freq_in_dataset >= math.ceil(nb_sample*0.5) :
					comment = comment + "NEW HOTSPOT" + " | "
				else :
					comment = "NEW" + " | "

			af  = []

			#AF 1000g
			try :
				if record[g1000_key] is not None :
					if re.search('-?\d+\.?\d*',record[g1000_key]) :
						g1000 = float(record[g1000_key])
						af.append(g1000)
					else :
						g1000 = None
				else :
					g1000 = None
			except KeyError :
				if record["1000g2015aug_all"] is not None :
					if re.search('-?\d+\.?\d*',record["1000g2015aug_all"]) :
						g1000 = float(record["1000g2015aug_all"])
						af.append(g1000)
					else :
						g1000 = None
				else :
					g1000 = None

			#AF ESP6500
			try :
				if record["esp6500siv2_all"] is not None :
					if re.search('-?\d+\.?\d*',record["esp6500siv2_all"]) :
						esp6500 = float(record["esp6500siv2_all"])
						af.append(esp6500)
					else :
						esp6500 = None
				else :
					esp6500 = None
			except KeyError :
				if record["esp6500_all"] is not None :
					if re.search('-?\d+\.?\d*',record["esp6500_all"]) :
						esp6500 = float(record["esp6500_all"])
						af.append(esp6500)
					else :
						esp6500 = None
				else :
					esp6500 = None

			#AF exac
			try :
				if record["ExAC_nontcga_ALL"] is not None :
					if re.search('-?\d+\.?\d*',record["ExAC_nontcga_ALL"]) :
						exac = float(record["ExAC_nontcga_ALL"])
						af.append(exac)
					else :
						exac = None
				else :
					exac = None
			except KeyError :
				exac = None

			#AF Kaviar
			try :
				if record["Kaviar_AF"] is not None :
					if re.search('-?\d+\.?\d*',record["Kaviar_AF"]) :
						kaviar = float(record["Kaviar_AF"])
						af.append(kaviar)
					else :
						kaviar = None
				else :
					kaviar = None
			except KeyError :
				kaviar = None

			#AF HRC
			try :
				if record["HRC_AF"] is not None :
					if re.search('-?\d+\.?\d*',record["HRC_AF"]) :
						hrc = float(record["HRC_AF"])
						af.append(hrc)
					else :
						hrc = None
				else :
					hrc = None
			except KeyError :
				hrc = None

			#AF avsnp
			try :
				if record["AF"] is not None :
					if re.search('-?\d+\.?\d*',record["AF"]) :
						af_avsnp = float(record["AF"])
						af.append(af_avsnp)
					else :
						af_avsnp = None
				else :
					af_avsnp = None
			except KeyError :
				af_avsnp = None

			try :
				af_max = max(af)
			except ValueError :
				af_max = None

			#prediction of functional effects of human nsSNPs
			v = [".","."]
			if metaRNN_found is True :
				if record[metaRNN_key] == "D" :
					majority = True
				else :
					majority = False
			else :
				v = vote(record)
				try :
					majority = float(v[0])/float(v[1]) > 0.5
				except ZeroDivisionError :
					majority = False

			#DBSNP
			if rs129 is not None or rs138 is not None or rs is not None :
				if af_max is None :
					comment = comment + "SNP" + " | "
				elif af_max > 0.01 :
					comment = comment + "COMMON SNP" + " | "
				elif af_max <= 0.01 :
					comment = comment + "RARE SNP" + " | "

			#COSMIC
			if nb_cos > 0 and nb_cos < 10 :
				comment = comment + "COSMIC" + " | "
			elif nb_cos >= 10 :
				comment = comment + "COSMIC HOTSPOT" + " | "

			#PREDICTOR MAJORITY RULE
			if majority :
				effect = "DAMAGING"
			else :
				effect = "NEUTRAL"

			for field in fieldnames :
				output.write("%s\t"%record[field])
			output.write(comment+"\t"+str(v[0])+"\t"+str(v[1])+"\t"+effect+"\t"+str(freq_in_dataset))
			output.write("\n")

	output.close()
	variants.close()
