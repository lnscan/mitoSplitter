<h2 align="center">mitoSplitter</h2>

<p align="center">A mitochondrial variants-based method for efficient demultiplexing of pooled single-cell RNA-seq.</p>

---

  - [Install](#Install)
  - [Usage](#Usage)
  - [Example](#Example)
  - [Citation](#Citation)


## Install
```bash
git clone https://github.com/lnscan/mitoSplitter.git
conda create -n mitoSplitter python==3.9.13
conda activate mitoSplitter
pip install -r requirement.txt
```
Make sure samtools can be used directly.

## Usage
```bash
sh mitoSplitter_pipeline.sh -h ## show help information
```
```txt
Usage:
	sh mitoSplitter_pipeline.sh <-i input.bam> <-b barcode.list> <-o out_dir> [-r bulk_matrix] [-l bulk_bam_list] [-s cor_value] [-m mito.fasta] [-t threads] [-f barcode_tag] [-q base_quality] [-a alignment_quality] [-g gold_file] [-p matrix_dir] [-d] [-h]
 
Description:
	A mitochondrial variants-based method for efficient demultiplexing of pooled single-cell RNA-seq.
 
Ordering options:
	-i	input bam file
	-b	barcode list file, one barcode per line
	-o	name of directory for mitoSplitter output files
 

Other options:
	-r	mitochondrial variant allele frequency matrix file for all samples, must be set unless -l is set
	-l	bulk RNA-seq mapping bam list file for all samples, one bam file per line, must be set unless -r is set
	-s	warning threshold for correlation between bulk samples, used to signal samples with similar genetic backgrounds, default = 0.65
	-m	mitochondrial reference genome fasta, default = GRCh38_MT.fasta
	-c	mitochondrial reference genome id, default = MT
	-t	max threads to use, default = 50
	-f	used to extract cell barcode from bam file, default = CB
	-q	minimum base quality to be considered in bam file, used for variant calling, default = 20
	-a	minimum alignment quality to be considered in bam file, used for variant calling, default = 20
	-g	benchmark file used for performance validation, one barcode and the corresponding cluster per line
	-p	filtered feature bc matrix directory for removing doublets using scrublet, python package scrublet needs to be installed, default = FALSE
	-d	remove doublets using Favg Gaussian fitted model, R package mixtools needs to be installed, default = FALSE
	-h	print this help and exit
```

## Example
```bash
cd example_data
sh test.sh ## result in ./example_result
```

## Citation
Lin X, Chen Y, Lin L, Yin K, Cheng R, Lin X, Wang X, Guo Y, Wu Z, Zhang Y, Li J, Yang C, Song J. mitoSplitter: A mitochondrial variants-based method for efficient demultiplexing of pooled single-cell RNA-seq. Proc Natl Acad Sci U S A. 2023 Sep 26;120(39):e2307722120. doi: 10.1073/pnas.2307722120. Epub 2023 Sep 19. PMID: 37725654; PMCID: PMC10523499.


