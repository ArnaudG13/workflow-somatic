#!/usr/bin/env python

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

def parse_FreeBayesindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_sample_normal
			indels[chrid]['ad_tumor']=ad_sample_tumor			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_VarScan2indels(vcf):
	indels = {}
	datacolumn = {}
	ss_dict = {"0":"Reference", "1":"Germline", "2":"Somatic", "3":"LOH", "5":"Unknown"}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			rd_normal = info[9].split(":")[3]
			ad_normal = info[9].split(":")[4]
			ad_sample_normal = rd_normal + "," + ad_normal
			rd_tumor = info[10].split(":")[3]
			ad_tumor = info[10].split(":")[4]
			ad_sample_tumor = rd_tumor + "," + ad_tumor
			somatic_score = info[7].split(";")[4]
			somatic_score = somatic_score.split("=")[1]
			statut = info[7].split(";")[3]
			statut = statut.split("=")[1]
			statut = ss_dict[statut]
			qual = somatic_score
			filt = statut
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_sample_normal
			indels[chrid]['ad_tumor']=ad_sample_tumor	
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_VarDictindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			if info[4] != "<DEL>" and info[4] != "<INV>" and info[4] != "<DUP>" :
				chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
				ad_sample_normal = info[9].split(":")[1]
				ad_sample_tumor = info[10].split(":")[1]
				status = info[7]
				status = status.split(";")[10]
				status = status.split("=")[1]
				qual = info[5]
				filt = status
				indels[chrid] = {}
				indels[chrid]['ad_normal']=ad_sample_normal
				indels[chrid]['ad_tumor']=ad_sample_tumor			
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
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			val = info[7]
			tlod = val.split(";")[-1]
			tlod = tlod.split("=")[1]
			qual = tlod
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_sample_normal
			indels[chrid]['ad_tumor']=ad_sample_tumor			
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
			indels[chrid]['ad_tumor']=ad_sample
			indels[chrid]['ad_normal']="."			
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
			ad_tumor = info[10].split(":")[4]
			ad_normal = info[9].split(":")[4]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_normal
			indels[chrid]['ad_tumor']=ad_tumor			
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
			ad_normal = info[9].split(":")[1]
			ad_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_normal
			indels[chrid]['ad_tumor']=ad_tumor			
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
			ad_sample_normal = "."
			ad_sample_tumor = "."
			qual = info[5]
			filt = info[6]
			somatic_evs = info[7]
			somatic_evs = somatic_evs.split(";")[12]
			somatic_evs = somatic_evs.split("=")[1]
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_sample_normal
			indels[chrid]['ad_tumor']=ad_sample_tumor			
			indels[chrid]['qual']=somatic_evs
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Seuratindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			val = info[7]
			reg_ar1="AR1=(-?\d+\.?\d*)"
			reg_ar2="AR2=(-?\d+\.?\d*)"
			reg_dp1="DP1=(-?\d+)"
			reg_dp2="DP2=(-?\d+)"
			ar1 = re.search(reg_ar1,val).group(1) if re.search(reg_ar1,val) else None
			ar2 = re.search(reg_ar2,val).group(1) if re.search(reg_ar2,val) else None
			dp1 = re.search(reg_dp1,val).group(1) if re.search(reg_dp1,val) else None
			dp2 = re.search(reg_dp2,val).group(1) if re.search(reg_dp2,val) else None
			if ar1 is not None and dp1 is not None and ar2 is not None and dp2 is not None :
				alt_in_normal = int(float(ar1)*int(dp1))
				alt_in_tumor = int(float(ar2)*int(dp2))
				ref_in_normal = int(dp1) - alt_in_normal
				ref_in_tumor = int(dp2) - alt_in_tumor
				ad_sample_normal = str(ref_in_normal) + "," + str(alt_in_normal)
				ad_sample_tumor = str(ref_in_tumor) + "," + str(alt_in_tumor)
			else :
				ad_sample_normal = "."
				ad_sample_tumor = "."
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_sample_normal
			indels[chrid]['ad_tumor']=ad_sample_tumor			
			indels[chrid]['qual']=qual
			indels[chrid]['filter']=filt

	return {'indels':indels}

def parse_Lancetindels(vcf):
	indels = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			indels[chrid] = {}
			indels[chrid]['ad_normal']=ad_sample_normal
			indels[chrid]['ad_tumor']=ad_sample_tumor			
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

def mergeindels(freebayes_indels, lofreq_indels, mutect2_indels, pindel_indels, scalpel_indels, seurat_indels ,strelka_indels, vardict_indels, varscan2_indels, lancet_indels, output, n_concordant):
	all_indels = list()
	sf = open(output,"w")
	sf.write("%s\n" %("##fileformat=VCFv4.2"))
	sf.write("%s%s\n" %("##date=",str(datetime.now())))
	sf.write("%s\n" %("##source=MergeCaller"))
	if freebayes_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=FreeBayes,Description=\"Called by FreeBayes\">"))
		all_indels = all_indels + list(freebayes_indels['indels'].keys())
	if lancet_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Lancet,Description=\"Called by Lancet\">"))
		all_indels = all_indels + list(lancet_indels['indels'].keys())
	if lofreq_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Lofreq,Description=\"Called by LoFreq\">"))
		all_indels = all_indels + list(lofreq_indels['indels'].keys())
	if mutect2_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Mutect2,Description=\"Called by Mutect2\">"))
		all_indels = all_indels + list(mutect2_indels['indels'].keys())
	if pindel_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Pindel,Description=\"Called by Pindel\">"))
		all_indels = all_indels + list(pindel_indels['indels'].keys())
	if scalpel_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Scalpel,Description=\"Called by Scalpel\">"))
		all_indels = all_indels + list(scalpel_indels['indels'].keys())
	if seurat_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Seurat,Description=\"Called by Seurat\">"))
		all_indels = all_indels + list(seurat_indels['indels'].keys())
	if strelka_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Strelka,Description=\"Called by Strelka\">"))
		all_indels = all_indels + list(strelka_indels['indels'].keys())
	if vardict_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Vardict,Description=\"Called by VarDict\">"))
		all_indels = all_indels + list(vardict_indels['indels'].keys())
	if varscan2_indels is not None :
		sf.write("%s\n" %("##FILTER=<ID=Varscan2,Description=\"Called by Varscan2\">"))
	sf.write("%s\n" %("##INFO=<ID=VAF_NORMAL,Number=1,Type=Float,Description=\"Median vaf between callers in normal\">"))
	sf.write("%s\n" %("##INFO=<ID=VAF_TUMOR,Number=1,Type=Float,Description=\"Median vaf between callers in tumor\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADP1,Number=R,Type=Integer,Description=\"Allelic depths reported by FreeBayes for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADLF,Number=R,Type=Integer,Description=\"Allelic depths reported by LoFreq for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADM2,Number=R,Type=Integer,Description=\"Allelic depths reported by Mutect2 for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADPI,Number=R,Type=Integer,Description=\"Allelic depths reported by Pindel for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADSC,Number=R,Type=Integer,Description=\"Allelic depths reported by Scalpel for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADSE,Number=R,Type=Integer,Description=\"Allelic depths reported by Seurat for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADSK,Number=R,Type=Integer,Description=\"Allelic depths reported by Strelka for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADVC,Number=R,Type=Integer,Description=\"Allelic depths reported by Vardict for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADVS2,Number=R,Type=Integer,Description=\"Allelic depths reported by Varscan2 for the ref and alt alleles in the order listed\">"))
	sf.write("%s\n" %("##FORMAT=<ID=ADLA,Number=R,Type=Integer,Description=\"Allelic depths reported by Lancet for the ref and alt alleles in the order listed\">"))
	sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('#CHROM', 'POS','ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO','FORMAT', "NORMAL", "TUMOR"))
	all_indels = sorted(set(all_indels))
	for indels in all_indels :
		vcfinfo = {}
		if indels in freebayes_indels['indels'] :
			vcfinfo['freebayes']=indels
		if indels in lofreq_indels['indels'] :
			vcfinfo['lofreq']=indels
		if indels in mutect2_indels['indels'] :
			vcfinfo['mutect2']=indels
		if indels in pindel_indels['indels'] :
			vcfinfo['pindel']=indels
		if indels in scalpel_indels['indels'] :
			vcfinfo['scalpel']=indels
		if indels in seurat_indels['indels'] :
			vcfinfo['seurat']=indels
		if indels in strelka_indels['indels'] :
			vcfinfo['strelka']=indels
		if indels in vardict_indels['indels'] :
			vcfinfo['vardict']=indels
		if indels in varscan2_indels['indels'] :
			vcfinfo['varscan2']=indels
		if indels in lancet_indels['indels'] :
			vcfinfo['lancet']=indels
		called_by = list(vcfinfo.keys())
		if all(value == vcfinfo[called_by[0]] for value in vcfinfo.values()):
			format=''
			gf_normal=''
			gf_tumor=''
			callers=''
			nb_callers_pass=0
			af_normal = []
			af_tumor = []
			for c in called_by :
				if c=='freebayes':
					format=format+'ADP1:'
					gf_normal=gf_normal+freebayes_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(freebayes_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+freebayes_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(freebayes_indels['indels'][indels]['ad_tumor']))
					if freebayes_indels['indels'][indels]['filter'] != "REJECT" :
						nb_callers_pass += 1
						callers=callers+'FreeBayes|'
				elif c=='lofreq':
					#filter1=re.compile('min_dp_10')
					#filter2=re.compile('sb_fdr')
					#filter3=re.compile('min_snvqual_76')
					#filter4=re.compile('min_indelqual_54')
					#f1=filter1.search(lofreq_snv['snvs'][snv]['filter'])
					#f2=filter2.search(lofreq_snv['snvs'][snv]['filter'])
					#f3=filter3.search(lofreq_snv['snvs'][snv]['filter'])
					#f4=filter4.search(lofreq_snv['snvs'][snv]['filter'])
					format=format+'ADLF:'
					gf_normal=gf_normal+lofreq_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(lofreq_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+lofreq_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(lofreq_indels['indels'][indels]['ad_tumor']))
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
					filter10=re.compile('germline')
					f1=filter1.search(mutect2_indels['indels'][indels]['filter'])
					f2=filter2.search(mutect2_indels['indels'][indels]['filter'])
					f3=filter3.search(mutect2_indels['indels'][indels]['filter'])
					f4=filter4.search(mutect2_indels['indels'][indels]['filter'])
					f5=filter5.search(mutect2_indels['indels'][indels]['filter'])
					f6=filter6.search(mutect2_indels['indels'][indels]['filter'])
					f7=filter7.search(mutect2_indels['indels'][indels]['filter'])
					f8=filter8.search(mutect2_indels['indels'][indels]['filter'])
					f9=filter9.search(mutect2_indels['indels'][indels]['filter'])
					f10=filter10.search(mutect2_indels['indels'][indels]['filter'])
					format=format+'ADM2:'
					gf_normal=gf_normal+mutect2_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(mutect2_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+mutect2_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(mutect2_indels['indels'][indels]['ad_tumor']))
					#if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10) :
					if mutect2_indels['indels'][indels]['filter']=="PASS":
						nb_callers_pass += 1
						callers=callers+'Mutect2|'
				elif c=='pindel':
					format=format+'ADPI:'
					gf_normal=gf_normal+pindel_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(pindel_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+pindel_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(pindel_indels['indels'][indels]['ad_tumor']))
					nb_callers_pass += 1
					callers=callers+'Pindel|'
				elif c=='scalpel':
					format=format+'ADSC:'
					gf_normal=gf_normal+scalpel_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(scalpel_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+scalpel_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(scalpel_indels['indels'][indels]['ad_tumor']))
					if scalpel_indels['indels'][indels]['filter'] == "PASS" :
						nb_callers_pass += 1
						callers=callers+'Scalpel|'
				elif c=='seurat' :
					format=format+'ADSE:'
					gf_normal=gf_normal+seurat_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(seurat_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+seurat_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(seurat_indels['indels'][indels]['ad_tumor']))
					qual = float(seurat_indels['indels'][indels]['qual'])
					if qual > 10 :
						nb_callers_pass += 1
						callers=callers+'Seurat|'
				elif c=='strelka':
					format=format+'ADSK:'
					gf_normal=gf_normal+strelka_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(strelka_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+strelka_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(strelka_indels['indels'][indels]['ad_tumor']))
					if strelka_indels['indels'][indels]['filter'] == "PASS" :
						nb_callers_pass += 1
						callers=callers+'Strelka|'
				elif c=='vardict':
					#filter1=re.compile('AMPBIAS')
					#filter2=re.compile('Bias')
					#filter3=re.compile('Cluster0bp')
					#filter4=re.compile('InGap')
					#filter5=re.compile('InIns')
					#filter6=re.compile('LongMSI')
					#filter7=re.compile('MSI12')
					#filter8=re.compile('NM5.25')
					#filter9=re.compile('Q10')
					#filter10=re.compile('SN1.5')
					#filter11=re.compile('d3')
					#filter12=re.compile('f0.02')
					#filter13=re.compile('p8')
					#filter14=re.compile('pSTD')
					#filter15=re.compile('q22.5')
					#filter16=re.compile('v2')
					#f1=filter1.search(vardict_indels['indels'][indels]['filter'])
					#f2=filter2.search(vardict_indels['indels'][indels]['filter'])
					#f3=filter3.search(vardict_indels['indels'][indels]['filter'])
					#f4=filter4.search(vardict_indels['indels'][indels]['filter'])
					#f5=filter5.search(vardict_indels['indels'][indels]['filter'])
					#f6=filter6.search(vardict_indels['indels'][indels]['filter'])
					#f7=filter7.search(vardict_indels['indels'][indels]['filter'])
					#f8=filter8.search(vardict_indels['indels'][indels]['filter'])
					#f9=filter9.search(vardict_indels['indels'][indels]['filter'])
					#f10=filter10.search(vardict_indels['indels'][indels]['filter'])
					#f11=filter11.search(vardict_indels['indels'][indels]['filter'])
					#f12=filter12.search(vardict_indels['indels'][indels]['filter'])
					#f13=filter13.search(vardict_indels['indels'][indels]['filter'])
					#f14=filter14.search(vardict_indels['indels'][indels]['filter'])
					#f15=filter15.search(vardict_indels['indels'][indels]['filter'])
					#f16=filter16.search(vardict_indels['indels'][indels]['filter'])
					format=format+'ADVC:'
					gf_normal=gf_normal+vardict_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(vardict_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+vardict_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(vardict_indels['indels'][indels]['ad_tumor']))
					#if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10 or f11 or f12 or f13 or f14 or f15 or f16) :
					#if vardict_indels['indels'][indels]['filter'] == 'StrongSomatic' or vardict_indels['indels'][indels]['filter'] == 'LikelySomatic' :
					if vardict_indels['indels'][indels]['filter'] == 'StrongSomatic' :
						callers=callers+'Vardict|'
						nb_callers_pass += 1
				elif c=='varscan2':
					format=format+'ADVS2:'
					gf_normal=gf_normal+varscan2_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(varscan2_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+varscan2_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(varscan2_indels['indels'][indels]['ad_tumor']))
					#if float(varscan2_indels['indels'][indels]['qual']) > 30 and varscan2_indels['indels'][indels]['filter'] == "Somatic" :
					if float(varscan2_indels['indels'][indels]['qual']) > 20 and varscan2_indels['indels'][indels]['filter'] == "Somatic" :
						nb_callers_pass += 1
						callers=callers+'Varscan2|'
				elif c=='lancet':
					format=format+'ADLA:'
					gf_normal=gf_normal+lancet_indels['indels'][indels]['ad_normal']+':'
					af_normal.append(get_af(lancet_indels['indels'][indels]['ad_normal']))
					gf_tumor=gf_tumor+lancet_indels['indels'][indels]['ad_tumor']+':'
					af_tumor.append(get_af(lancet_indels['indels'][indels]['ad_tumor']))
					if lancet_indels['indels'][indels]['filter'] == "PASS" :
						nb_callers_pass += 1
						callers=callers+'Lancet|'

			if nb_callers_pass > 0 :
				vaf_tumor = round(numpy.nanmedian(af_tumor),4)
				vaf_normal = round(numpy.nanmedian(af_normal),4)
				callers = callers[:-1]
				if vaf_tumor < 0.05 :
					filt = "LowVariantFreq|"+callers
				else :
					filt = callers
				if nb_callers_pass >= n_concordant :
					filt =  "CONCORDANT|"+filt
				else :
					filt = "DISCORDANT|"+filt
				info = "VAF_NORMAL="+str(vaf_normal)+";"+"VAF_TUMOR="+str(vaf_tumor)
				format = format[:-1]
				gf_normal = gf_normal[:-1]
				gf_tumor = gf_tumor[:-1]
				qual=nb_callers_pass
				vcfinfolist=vcfinfo[called_by[0]].split("\t")
				baseinfo=vcfinfolist[0]+'\t'+vcfinfolist[1]+'\t.\t'+vcfinfolist[2]+'\t'+vcfinfolist[3]
				sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(baseinfo,qual, filt, info, format, gf_normal, gf_tumor))
		else :
			print("Conflict in ref and alt alleles between callers at pos "+indels)

parser = argparse.ArgumentParser(description="Merge somatic indel vcf files from different variants callers")
parser.add_argument('--FreeBayes', type=str, required=False)
parser.add_argument('--Lancet', type=str, required=False)
parser.add_argument('--LoFreq', type=str, required=False)
parser.add_argument('--Mutect2', type=str, required=False)
parser.add_argument('--pindel', type=str, required=False)
parser.add_argument('--Scalpel', type=str, required=False)
parser.add_argument('--Seurat', type=str, required=False)
parser.add_argument('--Strelka', type=str, required=False)
parser.add_argument('--VarDict', type=str, required=False)
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

if args.Lancet is not None :
	lancet_indels = parse_Lancetindels(args.Lancet)
	n_vc = n_vc + 1
else :
	lancet_indels = None

if args.LoFreq is not None :	
	lofreq_indels = parse_LoFreqindels(args.LoFreq)
	n_vc = n_vc + 1
else :
	lofreq_indels = None

if args.Mutect2 is not None :
	mutect2_indels = parse_Mutect2indels(args.Mutect2)
	n_vc = n_vc + 1
else :
	mutect2_indels = None

if args.pindel is not None :
	pindel_indels = parse_Pindelindels(args.pindel)
	n_vc = n_vc + 1
else :
	pindel_indels = None

if args.Scalpel is not None :
	scalpel_indels = parse_Scalpelindels(args.Scalpel)
	n_vc = n_vc + 1
else :
	scalpel_indels = None

if args.Seurat is not None :
	seurat_indels = parse_Seuratindels(args.Seurat)
	n_vc = n_vc + 1
else :
	seurat_indels = None

if args.Strelka is not None :
	strelka_indels = parse_Strelkaindels(args.Strelka)
	n_vc = n_vc + 1
else :
	strelka_indels = None

if args.VarDict is not None :
	vardict_indels = parse_VarDictindels(args.VarDict)
	n_vc = n_vc + 1
else :
	vardict_indels = None

if args.VarScan2 is not None :
	varscan2_indels = parse_VarScan2indels(args.VarScan2)
	n_vc = n_vc + 1
else :
	varscan2_indels = None

output = args.output

if n_concordant > n_vc :
	sys.exit("N concordant cannot be greater than the number of variant caller")

mergeindels(freebayes_indels, lofreq_indels, mutect2_indels, pindel_indels, scalpel_indels, seurat_indels ,strelka_indels, vardict_indels, varscan2_indels, lancet_indels, output, n_concordant)
