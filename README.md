# Title

A Snakemake workflow for the analysis of somatic variant (SNVs and indels) and copy number. (Tumor-normal mode).

# Description

* Trimming (option)
* Quality control
* Alignment (BWA mem)
* Deduplication (option)
* BQSR (option)
* somatic SNVs, up to 16 callers [Muse, Mutect2, Lofreq, Virmid, NeuSomatic, DeepSomatic, VarNet, Mutect, Strelka, Varscan2, Seurat, FreeBayes, Vardict, Lancet, SomaticSniper, Shimmer]
* somatic indels, up to 13 callers [NeuSomatic, DeepSomatic, Mutect2, Strelka, Lancet, VarNet, Lofreq, VarScan2, pindel, Vardict, Seurat, FreeBayes, Scalpel]
* germline SNVs, up to 6 callers [FreeBayes, HaplotypeCaller, Pisces, Platypus, Strelka, VarScan2]
* germline indels, up to 8 callers [FreeBayes, HaplotypeCaller, pindel, Pisces, Platypus, Scalpel, Strelka, VarScan2]
* VCF merging and majority rule with the scripts workflow/scripts/merge_caller_somatic.py and workflow/scripts/merge_caller_somatic_indel.py
* Variants annotations (Annovar)

### Dependencies

* singularity >= 3.7.1 https://github.com/sylabs/singularity
* snakemake >= 7.25.0 https://snakemake.readthedocs.io/en/stable/

### Config and log files

###Log files directory

 ```
/home/user/logs/cluster/snakemake/
```

###Cluster config file for Slurm

 ```
config/config.yaml
```

###config file for the somatic analysis (Annototation files etc ..)

 ```
config/config_somatic_hg19.yml
```


### Usage

```
snakemake \
--configfile /scratch/user/my_project/config_somatic_hg19.yml \
--use-singularity \
--singularity-prefix /scratch/user/singularity_cache \
--singularity-args "--bind /scratch/user/" \
--profile slurm
```


## Develop with

* [Python](https://www.python.org/) 
* [bash](http://git.savannah.gnu.org/cgit/bash.git)
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home)

## Author

* **Arnaud Guille**


