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
			ad_sample = info[9].split(":")[1]
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
			rd = info[9].split(":")[4]
			ad = info[9].split(":")[5]
			ad_sample = rd + "," + ad
			qual = info[5]
			filt = info[6]
			quality = str(-10*math.log10(max(float(info[9].split(":")[8]),2.2250738585072014e-308)))
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=quality
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
			ad_sample = info[9].split(":")[1]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad']=ad_sample			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_StrelkaSNV(vcf):
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

def mergeSNV(freebayes_snv, hc_snv, pisces_snv, platypus_snv, strelka_snv, varscan2_snv, output):
	SAMPLE = os.path.basename(sys.argv[1])
	filter1=re.compile('(.*).SNP.*')
	SAMPLE=filter1.search(sys.argv[1]).group(1)
	sf = open(output,"w")
	sf.write("%s\n" %("##fileformat=VCFv4.2"))
	sf.write("%s%s\n" %("##date=",str(datetime.now())))
	sf.write("%s\n" %("##source=MergeCaller"))
	sf.write("%s\n" %("##FILTER=<ID=FreeBayes,Description=\"Called by FreeBayes\""))
	sf.write("%s\n" %("##FILTER=<ID=HC,Description=\"Called by HaplotypeCaller\""))
	sf.write("%s\n" %("##FILTER=<ID=Pisces,Description=\"Called by Pisces\""))
	sf.write("%s\n" %("##FILTER=<ID=Platypus,Description=\"Called by Platypus\""))
	sf.write("%s\n" %("##FILTER=<ID=Strelka,Description=\"Called by Strelka\""))
	sf.write("%s\n" %("##FILTER=<ID=Varscan2,Description=\"Called by Varscan2\""))
	sf.write("%s\n" %("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Median vaf between callers\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADP1,Number=R,Type=Integer,Description=\"Allelic depths reported by FreeBayes for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADHC,Number=R,Type=Integer,Description=\"Allelic depths reported by HaplotypeCaller for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPS,Number=R,Type=Integer,Description=\"Allelic depths reported by Pisces for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADPY,Number=R,Type=Integer,Description=\"Allelic depths reported by Platypus for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADSK,Number=R,Type=Integer,Description=\"Allelic depths reported by Strelka for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADVS2,Number=R,Type=Integer,Description=\"Allelic depths reported by Varscan2 for the ref and alt alleles in the order listed\""))
	sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('#CHROM', 'POS','ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO','FORMAT', "SAMPLE"))
	all_snvs = sorted(set( list(freebayes_snv['snvs'].keys()) + list(hc_snv['snvs'].keys()) + list(pisces_snv['snvs'].keys()) + list(platypus_snv['snvs'].keys()) + list(strelka_snv['snvs'].keys()) + list(varscan2_snv['snvs'].keys()) ))
	for snv in all_snvs :
		vcfinfo = {}
		if snv in freebayes_snv['snvs'] :
			vcfinfo['freebayes']=snv
		if snv in hc_snv['snvs'] :
			vcfinfo['HC']=snv
		if snv in pisces_snv['snvs'] :
			vcfinfo['pisces']=snv
		if snv in platypus_snv['snvs'] :
			vcfinfo['platypus']=snv
		if snv in strelka_snv['snvs'] :
			vcfinfo['strelka']=snv
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
				elif c=='strelka':
					format=format+'ADSK:'
					gf_sample=gf_sample+strelka_snv['snvs'][snv]['ad']+':'
					af.append(get_af(strelka_snv['snvs'][snv]['ad']))
					if strelka_snv['snvs'][snv]['filter'] == 'PASS' :
						nb_callers_pass += 1
						callers=callers+'Strelka|'
				elif c=='varscan2':
					format=format+'ADVS2:'
					gf_sample=gf_sample+varscan2_snv['snvs'][snv]['ad']+':'
					af.append(get_af(varscan2_snv['snvs'][snv]['ad']))
					if float(varscan2_snv['snvs'][snv]['qual']) > 30 :
						nb_callers_pass += 1
						callers=callers+'Varscan2|'

			if nb_callers_pass > 0 :
				vaf = round(numpy.nanmedian(af),4)
				callers = callers[:-1]
				if vaf < 0.05 :
					filt = "LowVariantFreq|"+callers
				else :
					filt = callers
				if nb_callers_pass > 3 :
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

freebayes_snv = parse_FreeBayesSNV(sys.argv[1])
hc_snv = parse_HaplotypeCallerSNV(sys.argv[2])
pisces_snv = parse_PiscesSNV(sys.argv[3])
platypus_snv = parse_PlatypusSNV(sys.argv[4])
strelka_snv = parse_StrelkaSNV(sys.argv[5])
varscan2_snv = parse_VarScan2SNV(sys.argv[6])
output = sys.argv[7]

mergeSNV(freebayes_snv, hc_snv, pisces_snv, platypus_snv, strelka_snv, varscan2_snv, output)
