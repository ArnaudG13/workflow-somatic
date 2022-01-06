#!/bin/sh

input=$1
output=$2

echo -e "##fileformat=VCFv4.1
##INFO=<ID=AAAB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and AB in Tumor\">
##INFO=<ID=AABB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and BB in Tumor\">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR" > $output

awk -F "\t" 'NR!=1 {print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\t.\tAAAB=" $10 ";AABB=" $11 "\tRD:AD\t" $5 ":" $6 "\t" $7 ":" $8}' $input >> $output
