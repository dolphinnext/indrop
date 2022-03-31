#!/bin/sh

# genome_build = 'human_hg38_gencode_v34' //* @dropdown @options:"human_hg19, human_hg38_gencode_v28, human_hg38_gencode_v34, mouse_mm10, mouse_mm10_gencode_m25, rat_rn6_refseq, drosophila_melanogaster_dm3, custom"
# indrop_version = 'NextSeq550' //*  @options:"NextSeq550, NextSeq2000"

# If you're running for the first time, you need to download genome_data as follows:

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

nextflow main.nf \
--DOWNDIR ${dnextdata} \
--reads '*_R{1,2}.fastq.gz' \ 
--mate 'pair' \
--genome_build 'human_hg38_gencode_v34' \
--run_Single_Cell_Module = 'yes' \
--run_Tophat = 'no' \
--run_STAR = 'yes' \
--run_HISAT2 = 'no' \
--cutoff_for_reads_per_cell = '10' \
--run_Split_Fastq = 'no' \
--indrop_version = 'NextSeq550' \
--cutoff_reads_for_valid_cell "100" \
--mate_split "single" \
--determined_fastq "./*_{R1,R2,R3}_001.fastq.gz" \
--cellBarcodeFile "./scripts/gel_barcode1_list.txt"


