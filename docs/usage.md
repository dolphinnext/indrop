# dolphinnext/indrop: Usage

## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run dolphinnext/indrop -profile docker --DOWNDIR /path/to/save/genome_data --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build mouse_mm10
```

If you're running for the first time, you need to download genome_data as follows:

```bash
# Please replace /path/to/save with your directory and execute following commands:
dnextdata="/path/to/save"
wget -r --no-parent -R "index.html*"  https://galaxyweb.umassmed.edu/pub/dnext_data/singleCell/ -P ${dnextdata}/singleCell -l inf -nc -nH --cut-dirs=3

# You can select the genome_build you want to use and download with the following commands:
dnextdata="/path/to/save"
wget https://galaxyweb.umassmed.edu/pub/genome_data/human/human_hg19/ -P ${dnextdata}/genome_data/human/human_hg19/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 
wget https://galaxyweb.umassmed.edu/pub/genome_data/human/hg38_gencode_v28/ -P ${dnextdata}/genome_data/human/hg38_gencode_v28/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 
wget https://galaxyweb.umassmed.edu/pub/genome_data/human/hg38_gencode_v34/ -P ${dnextdata}/genome_data/human/hg38_gencode_v34/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 
wget https://galaxyweb.umassmed.edu/pub/genome_data/mouse/mm10/ -P ${dnextdata}/genome_data/mouse/mm10/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 
wget https://galaxyweb.umassmed.edu/pub/genome_data/mouse/mm10_gencode_m25/ -P ${dnextdata}/genome_data/mouse/mm10_gencode_m25/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 
wget https://galaxyweb.umassmed.edu/pub/genome_data/rat/rn6/ -P ${dnextdata}/genome_data/rat/rn6/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 
wget https://galaxyweb.umassmed.edu/pub/genome_data/d_melanogaster/dm3/ -P ${dnextdata}/genome_data/d_melanogaster/dm3/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 

nextflow run dolphinnext/indrop -profile docker --DOWNDIR ${dnextdata}/genome_data --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build mouse_mm10 
```


This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results 
.nextflow_log   # Log file from Nextflow
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version. In order to download latest version of the pipeline you need to run following command:

```bash
nextflow pull dolpinnext/indrop
```

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker, test` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from Dockerhub: [`dolphinnext/indrop`](http://hub.docker.com/r/dolphinnext/indrop/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub


### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{R1,R2,R3}.fastq' --mate 'triple'
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with triple read data, the path must use `{R1,R2,R3}` notation to specify reads.


### `--mate`
Two options (single or triple) available for `--mate` parameter. If you have single-end data, you need to specify as 'single' and for triple read data, you need to specify as 'triple'. For example:

```bash
--reads 'path/to/data/sample_*_{R1,R2,R3}.fastq' --mate 'triple'
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```

It is not possible to run a mixture of single-end and triple read files in one run.


## Reference genomes

### `--genome_build` 
There are 4 different species supported in the UMMS-Biocore references. To run the pipeline, you must specify which to use with the `--genome_build` flag.

List of genomes that are supported are:

* Human
  * `--genome_build human_hg19`
  * `--genome_build human_hg38_gencode_v28`
  * `--genome_build human_hg38_gencode_v34`
* Mouse
  * `--genome_build mouse_mm10`
  * `--genome_build mouse_mm10_gencode_m25`
* Rat
  * `--genome_build rat_rn6_refseq`
* D. melanogaster
  * `--genome_build drosophila_melanogaster_dm3`

Note: For new genome requests, please send e-mail to UMMS-Biocore(biocore@umassmed.edu).


