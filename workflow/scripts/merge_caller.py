#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from collections import OrderedDict
import numpy
from datetime import datetime
import re
import math
import argparse

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def parse_pindel(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			

	return {'indels':indels}


def parse_HaplotypeCaller(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			

	return {'indels':indels}

def parse_HaplotypeCallerSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample
			snvs[chrid]['qual']=qual	
			snvs[chrid]['filter']=filt		

	return {'snvs':snvs}

def parse_FreeBayesSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[2]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_VarScan2SNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			info_8 = info[8].split(":")
			info_9 = info[9].split(":")
			if "AD" in info_8 and "RD" in info_8 :
				AD=info_8.index("AD")
				RD=info_8.index("RD")
				ad = info_9[AD]
				rd = info_9[RD]
				ad_sample = rd + "," + ad
			else :
				ad_sample = ".,."
			qual = info[5]
			filt = info[6]
			quality = str(-10*math.log10(max(float(info[9].split(":")[8]),2.2250738585072014e-308)))
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=quality
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_VarDictSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[3]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_PlatypusSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			info_8 = info[8].split(":")
			info_9 = info[9].split(":")
			if "NR" in info_8 and "NV" in info_8 :
				NR=info_8.index("NR")
				NV=info_8.index("NV")
				try :
					depth = int(info_9[NR])
					alt_count = int(info_9[NV])
				except ValueError :
					depth = int(info_9[NR].split(",")[0])
					alt_count = int(info_9[NV].split(",")[0])
				ref_count = depth - alt_count
				ad_sample = str(ref_count) + "," + str(alt_count)
			else :
				ad_sample = ".,."
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_PiscesSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[2]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_Mutect2SNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			val = info[7]
			reg="TLOD=(-?\d+\.?\d*)"
			tlod = re.search(reg,val).group(1) if re.search(reg,val) else None
			qual = tlod
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_LoFreqSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			qual = info[5]
			filt = info[6]
			val = info[7]
			DP4 = val.split(";")[3]
			DP4 = DP4.split("=")[1]
			DP4 = DP4.split(",")
			ref_count = int(DP4[0]) + int(DP4[1])
			alt_count = int(DP4[2]) + int(DP4[3])
			ad_sample = str(ref_count) + "," + str(alt_count)
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def get_af(ad):
	try :
		ref_count = float(ad.split(",")[0])
		alt_count = float(ad.split(",")[1])
		af = alt_count/(alt_count+ref_count)
	except ValueError :
		af = numpy.nan
	except ZeroDivisionError :
		af = numpy.nan
	return af

def mergeSNV(freebayes_snv, hc_snv, lofreq_snv, mutect2_snv, pisces_snv, platypus_snv, vardict_snv, varscan2_snv, output, n_concordant):
	all_snvs = list()
	sf = open(output,"w")
	sf.write("%s\n" %("##fileformat=VCFv4.2"))
	sf.write("%s%s\n" %("##date=",str(datetime.now())))
	sf.write("%s\n" %("##source=MergeCaller"))
	if freebayes_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=FreeBayes,Description=\"Called by FreeBayes\""))
		all_snvs = all_snvs + list(freebayes_snv['snvs'].keys())
	if hc_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=HC,Description=\"Called by HaplotypeCaller\""))
		all_snvs = all_snvs + list(hc_snv['snvs'].keys())
	if lofreq_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=Lofreq,Description=\"Called by LoFreq\""))
		all_snvs = all_snvs + list(lofreq_snv['snvs'].keys())
	if mutect2_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=Mutect2,Description=\"Called by Mutect2\""))
		all_snvs = all_snvs + list(mutect2_snv['snvs'].keys())
	if pisces_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=Pisces,Description=\"Called by Pisces\""))
		all_snvs = all_snvs + list(pisces_snv['snvs'].keys())
	if platypus_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=Platypus,Description=\"Called by Platypus\""))
		all_snvs = all_snvs + list(platypus_snv['snvs'].keys())
	if vardict_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=Vardict,Description=\"Called by Vardict\""))
		all_snvs = all_snvs + list(vardict_snv['snvs'].keys())
	if varscan2_snv is not None :
		sf.write("%s\n" %("##FILTER=<ID=Varscan2,Description=\"Called by Varscan2\""))
		all_snvs = all_snvs + list(varscan2_snv['snvs'].keys())
	sf.write("%s\n" %("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Median vaf between callers\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADP1,Number=R,Type=Integer,Description=\"Allelic depths reported by FreeBayes for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADHC,Number=R,Type=Integer,Description=\"Allelic depths reported by HaplotypeCaller for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADLF,Number=R,Type=Integer,Description=\"Allelic depths reported by LoFreq for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADM2,Number=R,Type=Integer,Description=\"Allelic depths reported by Mutect2 for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPS,Number=R,Type=Integer,Description=\"Allelic depths reported by Pisces for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPY,Number=R,Type=Integer,Description=\"Allelic depths reported by Platypus for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADVC,Number=R,Type=Integer,Description=\"Allelic depths reported by Vardict for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADVS2,Number=R,Type=Integer,Description=\"Allelic depths reported by Varscan2 for the ref and alt alleles in the order listed\""))
	sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('#CHROM', 'POS','ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO','FORMAT', "SAMPLE"))
	all_snvs = sorted(set(all_snvs))
	for snv in all_snvs :
		vcfinfo = {}
		if freebayes_snv is not None :
			if snv in freebayes_snv['snvs'] :
				vcfinfo['freebayes']=snv
		if hc_snv is not None :
			if snv in hc_snv['snvs'] :
				vcfinfo['HC']=snv
		if lofreq_snv is not None :
			if snv in lofreq_snv['snvs'] :
				vcfinfo['lofreq']=snv
		if mutect2_snv is not None :
			if snv in mutect2_snv['snvs'] :
				vcfinfo['mutect2']=snv
		if pisces_snv is not None :
			if snv in pisces_snv['snvs'] :
				vcfinfo['pisces']=snv
		if platypus_snv is not None :
			if snv in platypus_snv['snvs'] :
				vcfinfo['platypus']=snv
		if vardict_snv is not None :
			if snv in vardict_snv['snvs'] :
				vcfinfo['vardict']=snv
		if varscan2_snv is not None :
			if snv in varscan2_snv['snvs'] :
				vcfinfo['varscan2']=snv
		called_by = list(vcfinfo.keys())
		if all(value == vcfinfo[called_by[0]] for value in vcfinfo.values()):
			format=''
			gf_sample=''
			callers=''
			nb_callers_pass=0
			af = []
			for c in called_by :
				if c=='freebayes':
					format=format+'ADP1:'
					gf_sample=gf_sample+freebayes_snv['snvs'][snv]['ad']+':'
					af.append(get_af(freebayes_snv['snvs'][snv]['ad']))
					if float(freebayes_snv['snvs'][snv]['qual']) > 20 :
						nb_callers_pass += 1
						callers=callers+'FreeBayes|'
				elif c=='HC':
					format=format+'ADHC:'
					gf_sample=gf_sample+hc_snv['snvs'][snv]['ad']+':'
					af.append(get_af(hc_snv['snvs'][snv]['ad']))
					if float(hc_snv['snvs'][snv]['qual']) > 100 :
						nb_callers_pass +=1
						callers=callers+'HC|'
				elif c=='lofreq':
					filter1=re.compile('min_dp_10')
					filter2=re.compile('sb_fdr')
					filter3=re.compile('min_snvqual_76')
					filter4=re.compile('min_indelqual_54')
					f1=filter1.search(lofreq_snv['snvs'][snv]['filter'])
					f2=filter2.search(lofreq_snv['snvs'][snv]['filter'])
					f3=filter3.search(lofreq_snv['snvs'][snv]['filter'])
					f4=filter4.search(lofreq_snv['snvs'][snv]['filter'])
					format=format+'ADLF:'
					gf_sample=gf_sample+lofreq_snv['snvs'][snv]['ad']+':'
					af.append(get_af(lofreq_snv['snvs'][snv]['ad']))
					if not (f1 or f2 or f3 or f4) :
						nb_callers_pass +=1
						callers=callers+'Lofreq|'
				elif c=='mutect2':
					filter1=re.compile('alt_allele_in_normal')
					filter2=re.compile('clustered_events')
					filter3=re.compile('germline_risk')
					filter4=re.compile('homologous_mapping_event')
					filter5=re.compile('multi_event_alt_allele_in_normal')
					filter6=re.compile('panel_of_normals')
					filter7=re.compile('str_contraction')
					filter8=re.compile('t_lod_fstar')
					filter9=re.compile('triallelic_site')
					f1=filter1.search(mutect2_snv['snvs'][snv]['filter'])
					f2=filter2.search(mutect2_snv['snvs'][snv]['filter'])
					f3=filter3.search(mutect2_snv['snvs'][snv]['filter'])
					f4=filter4.search(mutect2_snv['snvs'][snv]['filter'])
					f5=filter5.search(mutect2_snv['snvs'][snv]['filter'])
					f6=filter6.search(mutect2_snv['snvs'][snv]['filter'])
					f7=filter7.search(mutect2_snv['snvs'][snv]['filter'])
					f8=filter8.search(mutect2_snv['snvs'][snv]['filter'])
					f9=filter9.search(mutect2_snv['snvs'][snv]['filter'])
					format=format+'ADM2:'
					gf_sample=gf_sample+mutect2_snv['snvs'][snv]['ad']+':'
					af.append(get_af(mutect2_snv['snvs'][snv]['ad']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9) :
						nb_callers_pass += 1
						callers=callers+'Mutect2|'
				elif c=='pisces':
					filter1=re.compile('q30')
					filter2=re.compile('SB')
					filter3=re.compile('R5x9')
					filter4=re.compile('NC')
					f1=filter1.search(pisces_snv['snvs'][snv]['filter'])
					f2=filter2.search(pisces_snv['snvs'][snv]['filter'])
					f3=filter3.search(pisces_snv['snvs'][snv]['filter'])
					f4=filter4.search(pisces_snv['snvs'][snv]['filter'])
					format=format+'ADPS:'
					gf_sample=gf_sample+pisces_snv['snvs'][snv]['ad']+':'
					af.append(get_af(pisces_snv['snvs'][snv]['ad']))
					if not (f1 or f2 or f3 or f4) :
						nb_callers_pass += 1
						callers=callers+'Pisces|'
				elif c=='platypus':
					filter1=re.compile('GOF')
					filter2=re.compile('hp10')
					filter3=re.compile('Q20')
					filter4=re.compile('HapScore')
					filter5=re.compile('MQ')
					filter6=re.compile('strandBias')
					filter7=re.compile('SC')
					filter8=re.compile('QualDepth')
					filter9=re.compile('REFCALL')
					filter10=re.compile('QD')
					f1=filter1.search(platypus_snv['snvs'][snv]['filter'])
					f2=filter2.search(platypus_snv['snvs'][snv]['filter'])
					f3=filter3.search(platypus_snv['snvs'][snv]['filter'])
					f4=filter4.search(platypus_snv['snvs'][snv]['filter'])
					f5=filter5.search(platypus_snv['snvs'][snv]['filter'])
					f6=filter6.search(platypus_snv['snvs'][snv]['filter'])
					f7=filter7.search(platypus_snv['snvs'][snv]['filter'])
					f8=filter8.search(platypus_snv['snvs'][snv]['filter'])
					f9=filter9.search(platypus_snv['snvs'][snv]['filter'])
					f10=filter10.search(platypus_snv['snvs'][snv]['filter'])
					format=format+'ADPY:'
					gf_sample=gf_sample+platypus_snv['snvs'][snv]['ad']+':'
					af.append(get_af(platypus_snv['snvs'][snv]['ad']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10) :
						nb_callers_pass += 1
						callers=callers+'Platypus|'
				elif c=='vardict':
					filter1=re.compile('AMPBIAS')
					filter2=re.compile('Bias')
					filter3=re.compile('Cluster0bp')
					filter4=re.compile('InGap')
					filter5=re.compile('InIns')
					filter6=re.compile('LongMSI')
					filter7=re.compile('MSI12')
					filter8=re.compile('NM5.25')
					filter9=re.compile('Q10')
					filter10=re.compile('SN1.5')
					filter11=re.compile('d3')
					filter12=re.compile('f0.02')
					filter13=re.compile('p8')
					filter14=re.compile('pSTD')
					filter15=re.compile('q22.5')
					filter16=re.compile('v2')
					f1=filter1.search(vardict_snv['snvs'][snv]['filter'])
					f2=filter2.search(vardict_snv['snvs'][snv]['filter'])
					f3=filter3.search(vardict_snv['snvs'][snv]['filter'])
					f4=filter4.search(vardict_snv['snvs'][snv]['filter'])
					f5=filter5.search(vardict_snv['snvs'][snv]['filter'])
					f6=filter6.search(vardict_snv['snvs'][snv]['filter'])
					f7=filter7.search(vardict_snv['snvs'][snv]['filter'])
					f8=filter8.search(vardict_snv['snvs'][snv]['filter'])
					f9=filter9.search(vardict_snv['snvs'][snv]['filter'])
					f10=filter10.search(vardict_snv['snvs'][snv]['filter'])
					f11=filter11.search(vardict_snv['snvs'][snv]['filter'])
					f12=filter12.search(vardict_snv['snvs'][snv]['filter'])
					f13=filter13.search(vardict_snv['snvs'][snv]['filter'])
					f14=filter14.search(vardict_snv['snvs'][snv]['filter'])
					f15=filter15.search(vardict_snv['snvs'][snv]['filter'])
					f16=filter16.search(vardict_snv['snvs'][snv]['filter'])
					format=format+'ADVC:'
					gf_sample=gf_sample+vardict_snv['snvs'][snv]['ad']+':'
					af.append(get_af(vardict_snv['snvs'][snv]['ad']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10 or f11 or f12 or f13 or f14 or f15 or f16) :
						callers=callers+'Vardict|'
						nb_callers_pass += 1
				elif c=='varscan2':
					format=format+'ADVS2:'
					gf_sample=gf_sample+varscan2_snv['snvs'][snv]['ad']+':'
					af.append(get_af(varscan2_snv['snvs'][snv]['ad']))
					if float(varscan2_snv['snvs'][snv]['qual']) > 30 :
						nb_callers_pass += 1
						callers=callers+'Varscan2|'

			if nb_callers_pass > 0 :
				vaf = round(numpy.nanmedian(af),4) * 100
				callers = callers[:-1]
				if vaf < 5 :
					filt = "LowVariantFreq|"+callers
				else :
					filt = callers
				if nb_callers_pass > n_concordant :
					filt =  "CONCORDANT|"+filt
				else :
					filt = "DISCORDANT|"+filt
				info = "VAF="+str(vaf)
				format = format[:-1]
				gf_sample = gf_sample[:-1]
				qual=nb_callers_pass
				vcfinfolist=vcfinfo[called_by[0]].split("\t")
				vcfinfolist=vcfinfo[called_by[0]].split('\t')
				baseinfo=vcfinfolist[0]+'\t'+vcfinfolist[1]+'\t.\t'+vcfinfolist[2]+'\t'+vcfinfolist[3]
				sf.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(baseinfo,qual, filt, info, format, gf_sample))
		else :
			print("Conflict in ref and alt alleles between callers at pos "+indel)

parser = argparse.ArgumentParser(description="Merge snv vcf files from different variants callers")
parser.add_argument('--FreeBayes', type=str, required=False)
parser.add_argument('--HaplotypeCaller', type=str, required=False)
parser.add_argument('--LoFreq', type=str, required=False)
parser.add_argument('--Mutect2', type=str, required=False)
parser.add_argument('--Pisces', type=str, required=False)
parser.add_argument('--Platypus', type=str, required=False)
parser.add_argument('--VarDict', type=str, required=False)
parser.add_argument('--VarScan2', type=str, required=False)
parser.add_argument('-N',type=int, required=True, help="Number of vote to be concordant")
parser.add_argument('output', type=str)
args = parser.parse_args()

n_concordant = args.N
n_vc = 0

if args.FreeBayes is not None :
	freebayes_snv = parse_FreeBayesSNV(args.FreeBayes)
	n_vc = n_vc + 1
else :
	freebayes_snv = None

if args.HaplotypeCaller is not None :
	hc_snv = parse_HaplotypeCallerSNV(args.HaplotypeCaller)
	n_vc = n_vc + 1
else :
	HC_snv = None

if args.LoFreq is not None :
	lofreq_snv = parse_LoFreqSNV(args.LoFreq)
	n_vc = n_vc + 1
else :
	lofreq_snv = None

if args.Mutect2 is not None :
	mutect2_snv = parse_Mutect2SNV(args.Mutect2)
	n_vc = n_vc + 1
else :
	mutect2_snv = None

if args.Pisces is not None :
	pisces_snv = parse_PiscesSNV(args.Pisces)
	n_vc = n_vc + 1
else :
	pisces_snv = None

if args.Platypus is not None :
	platypus_snv = parse_PlatypusSNV(args.Platypus)
	n_vc = n_vc + 1
else :
	platypus_snv = None

if args.VarDict is not None :
	vardict_snv = parse_VarDictSNV(args.VarDict)
	n_vc = n_vc + 1
else :
	vardict_snv = None

if args.VarScan2 is not None :
	varscan2_snv = parse_VarScan2SNV(args.VarScan2)
	n_vc = n_vc + 1
else :
	varscan2_snv = None

output = args.output

if n_concordant > n_vc :
	sys.exit("N concordant cannot be greater than the number of variant caller")

mergeSNV(freebayes_snv, hc_snv, lofreq_snv, mutect2_snv, pisces_snv, platypus_snv, vardict_snv, varscan2_snv, output, n_concordant)
