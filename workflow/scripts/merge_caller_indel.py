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

def parse_VarDictindels(vcf):
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

def parse_Mutect2indels(vcf):
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
			val = info[7]
			reg="TLOD=(-?\d+\.?\d*)"
			tlod = re.search(reg,val).group(1) if re.search(reg,val) else None
			qual = tlod
			indels[chrid] = {}
			indels[chrid]['ad']=ad_sample			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_LoFreqindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			qual = info[5]
			filt = info[6]
			val = info[7]
			DP4 = val.split(";")[2]
			DP4 = DP4.split("=")[1]
			DP4 = DP4.split(",")
			ref_count = int(DP4[0]) + int(DP4[1])
			alt_count = int(DP4[2]) + int(DP4[3])
			ad_sample = str(ref_count) + "," + str(alt_count)
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

def mergeindels(freebayes_indels, hc_indels, lofreq_indels, mutect2_indels, pisces_indels, platypus_indels, vardict_indels, varscan2_indels, pindel_indels, scalpel_indels, output):
	SAMPLE = os.path.basename(sys.argv[1])
	filter1=re.compile('(.*).INDEL.*')
	SAMPLE=filter1.search(sys.argv[1]).group(1)
	sf = open(output,"w")
	sf.write("%s\n" %("##fileformat=VCFv4.2"))
	sf.write("%s%s\n" %("##date=",str(datetime.now())))
	sf.write("%s\n" %("##source=MergeCaller"))
	sf.write("%s\n" %("##FILTER=<ID=FreeBayes,Description=\"Called by FreeBayes\">"))
	sf.write("%s\n" %("##FILTER=<ID=HC,Description=\"Called by HaplotypeCaller\">"))
	sf.write("%s\n" %("##FILTER=<ID=Lofreq,Description=\"Called by LoFreq\">"))
	sf.write("%s\n" %("##FILTER=<ID=Mutect2,Description=\"Called by Mutect2\">"))
	sf.write("%s\n" %("##FILTER=<ID=Pisces,Description=\"Called by Pisces\">"))
	sf.write("%s\n" %("##FILTER=<ID=Platypus,Description=\"Called by Platypus\">"))
	sf.write("%s\n" %("##FILTER=<ID=Vardict,Description=\"Called by Vardict\">"))
	sf.write("%s\n" %("##FILTER=<ID=Varscan2,Description=\"Called by Varscan2\">"))
	sf.write("%s\n" %("##FILTER=<ID=Pindel,Description=\"Called by Pindel\">"))
	sf.write("%s\n" %("##FILTER=<ID=Scalpel,Description=\"Called by Scalpel\">"))
	sf.write("%s\n" %("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Median vaf between callers\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADP1,Number=R,Type=Integer,Description=\"Allelic depths reported by FreeBayes for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADHC,Number=R,Type=Integer,Description=\"Allelic depths reported by HaplotypeCaller for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADLF,Number=R,Type=Integer,Description=\"Allelic depths reported by LoFreq for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADM2,Number=R,Type=Integer,Description=\"Allelic depths reported by Mutect2 for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADPS,Number=R,Type=Integer,Description=\"Allelic depths reported by Pisces for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADPY,Number=R,Type=Integer,Description=\"Allelic depths reported by Platypus for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADVC,Number=R,Type=Integer,Description=\"Allelic depths reported by Vardict for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADVS2,Number=R,Type=Integer,Description=\"Allelic depths reported by Varscan2 for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADPI,Number=R,Type=Integer,Description=\"Allelic depths reported by Pindel for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADSC,Number=R,Type=Integer,Description=\"Allelic depths reported by Scalpel for the ref and alt alleles in the order listed\">"))
	sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('#CHROM', 'POS','ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO','FORMAT', "SAMPLE"))
	all_indels = sorted(set( list(freebayes_indels['indels'].keys()) + list(hc_indels['indels'].keys()) +  list(lofreq_indels['indels'].keys()) + list(mutect2_indels['indels'].keys()) + list(pisces_indels['indels'].keys()) + list(platypus_indels['indels'].keys()) + list(pindel_indels['indels'].keys()) + list(scalpel_indels['indels'].keys()) + list(vardict_indels['indels'].keys()) + list(varscan2_indels['indels'].keys()) ))
	for indels in all_indels :
		vcfinfo = {}
		if indels in freebayes_indels['indels'] :
			vcfinfo['freebayes']=indels
		if indels in hc_indels['indels'] :
			vcfinfo['HC']=indels
		if indels in lofreq_indels['indels'] :
			vcfinfo['lofreq']=indels
		if indels in mutect2_indels['indels'] :
			vcfinfo['mutect2']=indels
		if indels in pisces_indels['indels'] :
			vcfinfo['pisces']=indels
		if indels in platypus_indels['indels'] :
			vcfinfo['platypus']=indels
		if indels in vardict_indels['indels'] :
			vcfinfo['vardict']=indels
		if indels in varscan2_indels['indels'] :
			vcfinfo['varscan2']=indels
		if indels in pindel_indels['indels'] :
			vcfinfo['pindel']=indels
		if indels in scalpel_indels['indels'] :
			vcfinfo['scalpel']=indels
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
				elif c=='lofreq':
					filter1=re.compile('min_dp_10')
					filter2=re.compile('sb_fdr')
					filter3=re.compile('min_indelsqual_76')
					filter4=re.compile('min_indelqual_54')
					f1=filter1.search(lofreq_indels['indels'][indels]['filter'])
					f2=filter2.search(lofreq_indels['indels'][indels]['filter'])
					f3=filter3.search(lofreq_indels['indels'][indels]['filter'])
					f4=filter4.search(lofreq_indels['indels'][indels]['filter'])
					format=format+'ADLF:'
					gf_sample=gf_sample+lofreq_indels['indels'][indels]['ad']+':'
					af.append(get_af(lofreq_indels['indels'][indels]['ad']))
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
					f1=filter1.search(mutect2_indels['indels'][indels]['filter'])
					f2=filter2.search(mutect2_indels['indels'][indels]['filter'])
					f3=filter3.search(mutect2_indels['indels'][indels]['filter'])
					f4=filter4.search(mutect2_indels['indels'][indels]['filter'])
					f5=filter5.search(mutect2_indels['indels'][indels]['filter'])
					f6=filter6.search(mutect2_indels['indels'][indels]['filter'])
					f7=filter7.search(mutect2_indels['indels'][indels]['filter'])
					f8=filter8.search(mutect2_indels['indels'][indels]['filter'])
					f9=filter9.search(mutect2_indels['indels'][indels]['filter'])
					format=format+'ADM2:'
					gf_sample=gf_sample+mutect2_indels['indels'][indels]['ad']+':'
					af.append(get_af(mutect2_indels['indels'][indels]['ad']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9) :
						nb_callers_pass += 1
						callers=callers+'Mutect2|'
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
					f1=filter1.search(vardict_indels['indels'][indels]['filter'])
					f2=filter2.search(vardict_indels['indels'][indels]['filter'])
					f3=filter3.search(vardict_indels['indels'][indels]['filter'])
					f4=filter4.search(vardict_indels['indels'][indels]['filter'])
					f5=filter5.search(vardict_indels['indels'][indels]['filter'])
					f6=filter6.search(vardict_indels['indels'][indels]['filter'])
					f7=filter7.search(vardict_indels['indels'][indels]['filter'])
					f8=filter8.search(vardict_indels['indels'][indels]['filter'])
					f9=filter9.search(vardict_indels['indels'][indels]['filter'])
					f10=filter10.search(vardict_indels['indels'][indels]['filter'])
					f11=filter11.search(vardict_indels['indels'][indels]['filter'])
					f12=filter12.search(vardict_indels['indels'][indels]['filter'])
					f13=filter13.search(vardict_indels['indels'][indels]['filter'])
					f14=filter14.search(vardict_indels['indels'][indels]['filter'])
					f15=filter15.search(vardict_indels['indels'][indels]['filter'])
					f16=filter16.search(vardict_indels['indels'][indels]['filter'])
					format=format+'ADVC:'
					gf_sample=gf_sample+vardict_indels['indels'][indels]['ad']+':'
					af.append(get_af(vardict_indels['indels'][indels]['ad']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10 or f11 or f12 or f13 or f14 or f15 or f16) :
						callers=callers+'Vardict|'
						nb_callers_pass += 1
				elif c=='varscan2':
					format=format+'ADVS2:'
					gf_sample=gf_sample+varscan2_indels['indels'][indels]['ad']+':'
					af.append(get_af(varscan2_indels['indels'][indels]['ad']))
					if float(varscan2_indels['indels'][indels]['qual']) > 30 :
						nb_callers_pass += 1
						callers=callers+'Varscan2|'
				elif c=='pindel':
					format=format+'ADVPI:'
					gf_sample=gf_sample+pindel_indels['indels'][indels]['ad']+':'
					af.append(get_af(pindel_indels['indels'][indels]['ad']))
					nb_callers_pass += 1
					callers=callers+'Pindel|'
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

			if nb_callers_pass > 0 :
				vaf = vaf = round(numpy.nanmedian(af),4)
				callers = callers[:-1]
				if vaf < 0.05 :
					filt = "LowVariantFreq|"+callers
				else :
					filt = callers
				if nb_callers_pass >= 5 :
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

freebayes_indels = parse_FreeBayesindels(sys.argv[1])
hc_indels = parse_HaplotypeCallerindels(sys.argv[2])
lofreq_indels = parse_LoFreqindels(sys.argv[3])
mutect2_indels = parse_Mutect2indels(sys.argv[4])
pisces_indels = parse_Piscesindels(sys.argv[5])
platypus_indels = parse_Platypusindels(sys.argv[6])
vardict_indels = parse_VarDictindels(sys.argv[7])
varscan2_indels = parse_VarScan2indels(sys.argv[8])
pindel_indels = parse_Pindelindels(sys.argv[9])
scalpel_indels = parse_Scalpelindels(sys.argv[10])
output = sys.argv[11]

mergeindels(freebayes_indels, hc_indels, lofreq_indels, mutect2_indels, pisces_indels, platypus_indels, vardict_indels, varscan2_indels, pindel_indels, scalpel_indels, output)
