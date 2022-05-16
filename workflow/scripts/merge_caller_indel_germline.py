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


def parse_HaplotypeCallerindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample
			indels[chrid]['qual']=qual	
			indels[chrid]['filter']=filt		

	return {'indels':indels}

def parse_FreeBayesindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_VarScan2indels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			rd = info[9].split(":")[4]
			ad = info[9].split(":")[5]
			ad_sample = rd + "," + ad
			qual = info[5]
			filt = info[6]
			quality = str(-10*math.log10(max(float(info[9].split(":")[8]),2.2250738585072014e-308)))
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=quality
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Platypusindels(vcf):
	indels = {}
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
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Piscesindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Pindelindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Scalpelindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Strelkaindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample		
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

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

def mergeindels(freebayes_indels, hc_indels, pindel_indels, pisces_indels, platypus_indels, scalpel_indels, strelka_indels, varscan2_indels, output, n_concordant):
	all_indels = list()
	sf = open(output,"w")
	sf.write("%s\n" %("##fileformat=VCFv4.2"))
	sf.write("%s%s\n" %("##date=",str(datetime.now())))
	sf.write("%s\n" %("##source=MergeCaller"))
	if freebayes_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=FreeBayes,Description=\"Called by FreeBayes\""))
		all_indels = all_indels + list(freebayes_indels['indels'].keys())
	if hc_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=HC,Description=\"Called by HaplotypeCaller\""))
		all_indels = all_indels + list(hc_indels['indels'].keys())
	if pindel_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Pindel,Description=\"Called by Pindel\""))
		all_indels = all_indels + list(pindel_indels['indels'].keys())
	if pisces_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Pisces,Description=\"Called by Pisces\""))
		all_indels = all_indels + list(pisces_indels['indels'].keys())
	if platypus_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Platypus,Description=\"Called by Platypus\""))
		all_indels = all_indels + list(platypus_indels['indels'].keys())
	if strelka_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Strelka,Description=\"Called by Strelka\""))
		all_indels = all_indels + list(strelka_indels['indels'].keys())
	if scalpel_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Scalpel,Description=\"Called by Scalpel\""))
		all_indels = all_indels + list(scalpel_indels['indels'].keys())
	if varscan2_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=VarScan2,Description=\"Called by VarScan2\""))
		all_indels = all_indels + list(varscan2_indels['indels'].keys())
	sf.write("%s\n" %("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Median vaf between callers\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADP1,Number=R,Type=Integer,Description=\"Allelic depths reported by FreeBayes for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADHC,Number=R,Type=Integer,Description=\"Allelic depths reported by HaplotypeCaller for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPI,Number=R,Type=Integer,Description=\"Allelic depths reported by Pindel for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPS,Number=R,Type=Integer,Description=\"Allelic depths reported by Pisces for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPY,Number=R,Type=Integer,Description=\"Allelic depths reported by Platypus for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADSC,Number=R,Type=Integer,Description=\"Allelic depths reported by Scalpel for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADSK,Number=R,Type=Integer,Description=\"Allelic depths reported by Strelka for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADVS2,Number=R,Type=Integer,Description=\"Allelic depths reported by Varscan2 for the ref and alt alleles in the order listed\""))
	sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('#CHROM', 'POS','ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO','FORMAT', "SAMPLE"))
	all_indels = sorted(set(all_indels))
	for indels in all_indels :
		vcfinfo = {}
		if freebayes_indels is not None :
			if indels in freebayes_indels['indels'] :
				vcfinfo['freebayes']=indels
		if hc_indels is not None :
			if indels in hc_indels['indels'] :
				vcfinfo['HC']=indels
		if pindel_indels is not None :
			if indels in pindel_indels['indels'] :
				vcfinfo['pindel']=indels
		if pisces_indels is not None :
			if indels in pisces_indels['indels'] :
				vcfinfo['pisces']=indels
		if platypus_indels is not None :
			if indels in platypus_indels['indels'] :
				vcfinfo['platypus']=indels
		if scalpel_indels is not None :
			if indels in scalpel_indels['indels'] :
				vcfinfo['scalpel']=indels
		if strelka_indels is not None :
			if indels in strelka_indels['indels'] :
				vcfinfo['strelka']=indels
		if varscan2_indels is not None :
			if indels in varscan2_indels['indels'] :
				vcfinfo['varscan2']=indels
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
					gf_sample=gf_sample+freebayes_indels['indels'][indels]['ad']+':'
					af.append(get_af(freebayes_indels['indels'][indels]['ad']))
					if float(freebayes_indels['indels'][indels]['qual']) > 20 :
						nb_callers_pass += 1
						callers=callers+'FreeBayes|'
				elif c=='HC':
					format=format+'ADHC:'
					gf_sample=gf_sample+hc_indels['indels'][indels]['ad']+':'
					af.append(get_af(hc_indels['indels'][indels]['ad']))
					if float(hc_indels['indels'][indels]['qual']) > 100 :
						nb_callers_pass +=1
						callers=callers+'HC|'
				elif c=='pindel':
					format=format+'ADVPI:'
					gf_sample=gf_sample+pindel_indels['indels'][indels]['ad']+':'
					af.append(get_af(pindel_indels['indels'][indels]['ad']))
					nb_callers_pass += 1
					callers=callers+'Pindel|'
				elif c=='pisces':
					filter1=re.compile('q30')
					filter2=re.compile('SB')
					filter3=re.compile('R5x9')
					filter4=re.compile('NC')
					f1=filter1.search(pisces_indels['indels'][indels]['filter'])
					f2=filter2.search(pisces_indels['indels'][indels]['filter'])
					f3=filter3.search(pisces_indels['indels'][indels]['filter'])
					f4=filter4.search(pisces_indels['indels'][indels]['filter'])
					format=format+'ADPS:'
					gf_sample=gf_sample+pisces_indels['indels'][indels]['ad']+':'
					af.append(get_af(pisces_indels['indels'][indels]['ad']))
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
					f1=filter1.search(platypus_indels['indels'][indels]['filter'])
					f2=filter2.search(platypus_indels['indels'][indels]['filter'])
					f3=filter3.search(platypus_indels['indels'][indels]['filter'])
					f4=filter4.search(platypus_indels['indels'][indels]['filter'])
					f5=filter5.search(platypus_indels['indels'][indels]['filter'])
					f6=filter6.search(platypus_indels['indels'][indels]['filter'])
					f7=filter7.search(platypus_indels['indels'][indels]['filter'])
					f8=filter8.search(platypus_indels['indels'][indels]['filter'])
					f9=filter9.search(platypus_indels['indels'][indels]['filter'])
					f10=filter10.search(platypus_indels['indels'][indels]['filter'])
					format=format+'ADPY:'
					gf_sample=gf_sample+platypus_indels['indels'][indels]['ad']+':'
					af.append(get_af(platypus_indels['indels'][indels]['ad']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10) :
						nb_callers_pass += 1
						callers=callers+'Platypus|'
				elif c=='scalpel':
					filter1=re.compile('LowAltCnt')
					filter2=re.compile('LowCoverage')
					f1=filter1.search(scalpel_indels['indels'][indels]['filter'])
					f2=filter2.search(scalpel_indels['indels'][indels]['filter'])
					format=format+'ADSC:'
					gf_sample=gf_sample+scalpel_indels['indels'][indels]['ad']+':'
					af.append(get_af(scalpel_indels['indels'][indels]['ad']))
					if not (f1 or f2) :
						nb_callers_pass += 1
						callers=callers+'Scalpel|'
				elif c=='strelka':
					format=format+'ADSK:'
					gf_sample=gf_sample+strelka_indels['indels'][indels]['ad']+':'
					af.append(get_af(strelka_indels['indels'][indels]['ad']))
					if strelka_indels['indels'][indels]['filter'] == 'PASS' :
						nb_callers_pass += 1
						callers=callers+'Strelka|'
				elif c=='varscan2':
					format=format+'ADVS2:'
					gf_sample=gf_sample+varscan2_indels['indels'][indels]['ad']+':'
					af.append(get_af(varscan2_indels['indels'][indels]['ad']))
					if float(varscan2_indels['indels'][indels]['qual']) > 30 :
						nb_callers_pass += 1
						callers=callers+'Varscan2|'
			
			if nb_callers_pass > 0 :	
				vaf = round(numpy.nanmedian(af),4) * 100
				callers = callers[:-1]
				if vaf < 5 :
					filt = "LowVariantFreq|"+callers
				else :
					filt = callers
				if nb_callers_pass >= n_concordant :
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


parser = argparse.ArgumentParser(description="Merge germline indel vcf files from different variants callers")
parser.add_argument('--FreeBayes', type=str, required=False)
parser.add_argument('--HaplotypeCaller', type=str, required=False)
parser.add_argument('--pindel', type=str, required=False)
parser.add_argument('--Pisces', type=str, required=False)
parser.add_argument('--Platypus', type=str, required=False)
parser.add_argument('--Scalpel', type=str, required=False)
parser.add_argument('--Strelka', type=str, required=False)
parser.add_argument('--VarScan2', type=str, required=False)
parser.add_argument('-N',type=int, required=True, help="Number of vote to be concordant")
parser.add_argument('output', type=str)
args = parser.parse_args()

n_concordant = args.N
n_vc = 0

if args.FreeBayes is not None :
	freebayes_indels = parse_FreeBayesindels(args.FreeBayes)
	n_vc = n_vc + 1
else :
	freebayes_indels = None

if args.HaplotypeCaller is not None :
	hc_indels = parse_HaplotypeCallerindels(args.HaplotypeCaller)
	n_vc = n_vc + 1
else :
	hc_indels = None

if args.pindel is not None :
	pindel_indels = parse_Pindelindels(args.pindel)
	n_vc = n_vc + 1
else :
	hc_indels = None

if args.Pisces is not None :
	pisces_indels = parse_Piscesindels(args.Pisces)
	n_vc = n_vc + 1
else :
	pisces_indels = None

if args.Platypus is not None :
	platypus_indels = parse_Platypusindels(args.Platypus)
	n_vc = n_vc + 1
else :
	platypus_indels = None

if args.Scalpel is not None :
	scalpel_indels = parse_Scalpelindels(args.Scalpel)
	n_vc = n_vc + 1
else :
	scalpel_indels = None

if args.Strelka is not None :
	strelka_indels = parse_Strelkaindels(args.Strelka)
	n_vc = n_vc + 1
else :
	strelka_indels = None

if args.VarScan2 is not None :
	varscan2_indels = parse_VarScan2indels(args.VarScan2)
	n_vc = n_vc + 1
else :
	varscan2_indels = None

output = args.output

if n_concordant > n_vc :
	sys.exit("N concordant cannot be greater than the number of variant caller")


mergeindels(freebayes_indels, hc_indels, pindel_indels ,pisces_indels, platypus_indels, scalpel_indels, strelka_indels ,varscan2_indels, output, n_concordant)
