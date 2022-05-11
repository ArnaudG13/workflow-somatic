# Titre du projet

Workflow snakemake pour l'analyse de variants somatiques et constitutionnels (SNP et INDEL), lorsque la contre-partie normale est disponible.

# Description du workflow

Le workflow intègre toutes les étapes essentielles d'une analyse de variants de type tumeur-normal.

* Trimming des fichiers fastq (optionnel)
* Controle qualité des fichiers fastq
* Alignement des reads
* Déduplication (optionnel)
* Recalibration
* Détection des SNP (13 outils au choix) et INDEL (10 outils au choix) somatiques
* Détection des SNP (6 outils au choix) et INDEL (8 outils au choix) constitutionnels
* Merge des fichiers VCF issus des différents outils et vote majoritaire
* Annotation des variants

### Prérequis outils

* singularity >= 3.5.3 https://github.com/sylabs/singularity
* snakemake >= 6.5.0 https://snakemake.readthedocs.io/en/stable/

### Fichiers de configuration

Le fichier ```config/config_somatic_hg19.yml``` contient tous les chemins vers les fichiers d'annotation, ainsi que les différentes options d'analyse.

Le fichier ```config/cluster2.yml``` contient les paramètres liés au cluster de calcul

###Fichier de log

Par défaut les fichiers de log du cluster de calcul sont générés dans

 ```
/home/guille/logs/cluster/snakemake/
```

Vous pouvez modifier le chemin de destination de ces logs dans le fichier config/cluster2.yml et créer le dossier s'il n'existe pas.


### Utilisation

Exemple de ligne de commande sur le cluster de l'ifb (slurm) à partir du dossier principale du workflow ```workflow-somatic```

```
snakemake --cores 400 -j 100 -T 3 \
--configfile /shared/projects/pmngs/projet_seq/benchmark/config_somatic_hg19.yml \
--use-singularity \
--singularity-prefix /shared/projects/pmngs/singularity_cache \
--cluster-config config/cluster2.yml \
--cluster "sbatch -N {cluster.host} -n {cluster.core} -t {cluster.time} --mem {cluster.mem} -o {cluster.stdout} -e {cluster.stderr}"
```

Exemple de ligne de commande sur le clsuter de l'ifb (slurm) :

```
nohup snakemake --cores 400 -j 100 -T 3 \
--configfile /shared/projects/pmngs/projet_seq/benchmark/config_somatic_hg19.yml \
--use-singularity \
--singularity-prefix /shared/projects/pmngs/singularity_cache \
--cluster-config config/cluster2.yml \
--cluster "sbatch -N {cluster.host} -n {cluster.core} -t {cluster.time} --mem {cluster.mem} -o {cluster.stdout} -e {cluster.stderr}"&
```

Exemple de ligne de commande sur le clsuter de l'ifb (slurm) :

```
nohup snakemake --cores 400 -j 100 -T 3 \
--configfile /shared/projects/pmngs/projet_seq/benchmark/config_somatic_hg19.yml \
--use-singularity \
--singularity-prefix /shared/projects/pmngs/singularity_cache \
--cluster-config config/cluster2.yml \
--cluster "sbatch -N {cluster.host} -n {cluster.core} -t {cluster.time} --mem {cluster.mem} -o {cluster.stdout} -e {cluster.stderr}"& > nohupoutput &
```

## Développé Avec

* [Python](https://www.python.org/) - Language de programmation Python
* [bash](http://git.savannah.gnu.org/cgit/bash.git) - Language de programmation Bash
* [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) - Gestionnaire de workflow snakemake

## Auteurs

* **Arnaud Guille**


