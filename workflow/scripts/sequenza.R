#!/usr/bin/env Rscript
#Script used to analyse allelic copy number with the package sequenza

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*10000)
library(sequenza)
library("scarHRD")

args = commandArgs(TRUE)
in_file = args[1]
out_dir = args[2]

sample_name = gsub("(.*).bin50.seqz.gz","\\1",basename(in_file))

dat = sequenza.extract(in_file, verbose=FALSE)
fit = sequenza.fit(dat)
sequenza.results(sequenza.extract = dat, cp.table=fit, sample.id=sample_name, out.dir=out_dir)
scar_score(in_file,reference="grch37", seqz=TRUE, output_dir=out_dir)
