# Title

A Snakemake workflow for the analysis of somatic variant (SNVs and indels) and copy number. (Tumor-normal mode).

# Description

* Trimming (option)
* Quality control
* Alignment (BWA mem)
* Deduplication (option)
* BQSR (option)
* somatic SNVs (13 callers) and indels (10 callers)
* germline SNVs (6 callers) and indels (8 callers)
* VCF merging and majority rule with the scripts workflow/scripts/merge_caller_somatic.py and workflow/scripts/merge_caller_somatic_indel.py
* Variants annotations (Annovar)

### Dependencies

* singularity >= 3.5.3 https://github.com/sylabs/singularity
* snakemake >= 6.5.0 https://snakemake.readthedocs.io/en/stable/

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

* [Python](https://www.python.org/) - Language de programmation Python
* [bash](http://git.savannah.gnu.org/cgit/bash.git) - Language de programmation Bash
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) - Gestionnaire de workflow snakemake

## Author

* **Arnaud Guille**


