#!/bin/sh

## Command parameters
func() {
	echo "Usage:"
	echo "	sh $0 <-i input.bam> <-b barcode.list> <-o out_dir> [-r bulk_matrix] [-l bulk_bam_list] [-s cor_value] [-m mito.fasta] [-t threads] [-f barcode_tag] [-q base_quality] [-a alignment_quality] [-g gold_file] [-p matrix_dir] [-d] [-h]"
	echo "Description:"
	echo "	A mitochondrial variants-based method for sample-demultiplexing of single-cell sequencing data"
	echo "Ordering options:"
	echo "	-i	input bam file"
	echo "	-b	barcode list file, one barcode per line"
	echo "	-o	name of directory for mitoSplitter output files"
	echo ""
	echo "Other options:"
	echo "	-r	mitochondrial variant allele frequency matrix file for all samples, must be set unless -l is set"
	echo "	-l	bulk RNA-seq mapping bam list file for all samples, one bam file per line, must be set unless -r is set"
	echo "	-s	warning threshold for correlation between bulk samples, used to signal samples with similar genetic backgrounds, default = 0.65"
	echo "	-m	mitochondrial reference genome fasta, default = GRCh38_MT.fasta"
	echo "	-c	mitochondrial reference genome id, default = MT"
	echo "	-t	max threads to use, default = 50"
	echo "	-f	used to extract cell barcode from bam file, default = CB"
	echo "	-q	minimum base quality to be considered in bam file, used for variant calling, default = 20"
	echo "	-a	minimum alignment quality to be considered in bam file, used for variant calling, default = 20"
	echo "	-g	benchmark file used for performance validation, one barcode and the corresponding cluster per line"
	echo "	-p	filtered feature bc matrix directory for removing doublets using scrublet, python package scrublet needs to be installed, default = FALSE"
	echo "	-d	remove doublets using Favg Gaussian fitted model, R package mixtools needs to be installed, default = FALSE"
	echo "	-h	print this help and exit"
	echo ""
	exit -1
}

CURRENT_DIR=`dirname $0`
MITOFA="${CURRENT_DIR}/mito_fastas/GRCh38_MT.fasta"
BLANK="${CURRENT_DIR}/scripts/blank.list"
MITOCHR="MT"
PREFIX="."
THREADS=50
BARCODE_TAG="CB"
BASEQUAL=20
ALIGNQUAL=20
CORWAR=0.65
BULKAF=""
BULKLIST=""
SCRUBLET="FALSE"
FAVG="FALSE"

while getopts 'i:b:o:r:l:s:m:c:t:f:q:a:g:p:dh' OPT;do
	case $OPT in
		i) BAMFILE="$OPTARG";;
		b) BARCODE="$OPTARG";;
		o) PREFIX="$OPTARG";;
		r) BULKAF="$OPTARG";;
		l) BULKLIST="$OPTARG";;
		s) WARCOR="$OPTARG";;
		m) MITOFA="$OPTARG";;
		c) MITOCHR="$OPTARG";;
		t) THREADS="$OPTARG";;
		f) BARCODE_TAG="$OPTARG";;
		q) BASEQUAL="$OPTARG";;
		a) ALIGNQUAL="$OPTARG";;
	    g) GOLD="$OPTARG";;
	    p) SCRUBLET="$OPTARG";;
	    d) FAVG="TRUE";;
		h) func;;
		?) func;;
	esac
done

## Set path to samtools if needed
#samtools=/path/to/samtools

mkdir -p $PREFIX

## Set or generate bulk mtRNA matrix
if [[ -n "$BULKAF" ]];then
    echo "Set bulk mtRNA matrix : $BULKAF"
	python ${CURRENT_DIR}/scripts/check_bulk_cor.py $BULKAF $CORWAR
elif [[ -n "$BULKLIST" ]];then
	echo "Generate bulk mtRNA matrix..."
	sh ${CURRENT_DIR}/scripts/generate_bulk_mtrna_matrix.sh -l $BULKLIST -o $PREFIX &
	BULKAF="$PREFIX/bulk_all.af.txt"
else
	echo "ERROR:Bulk mtRNA matrix path and bulk RNA-seq bam list file cannot both be missing!"
	echo ""
	func
	exit -1
fi

echo "Start mitoSplitter pipeline at `date`"
## Extract mitochondrial bam file
if [[ ! -e "${BAMFILE}.bai" ]];then
	samtools sort -@ $THREADS $BAMFILE -o ${PREFIX}/sorted.bam
	samtools index -@ $THREADS ${PREFIX}/sorted.bam
	samtools view -@ $THREADS -h ${PREFIX}/sorted.bam $MITOCHR -b -o ${PREFIX}/mt.bam ## make sure $BAMFILE is sorted
else
	samtools view -@ $THREADS -h $BAMFILE $MITOCHR -b -o ${PREFIX}/mt.bam ## make sure $BAMFILE is sorted
fi

## Remapping
python ${CURRENT_DIR}/scripts/renamer.py --bam ${PREFIX}/mt.bam --barcodes $BARCODE --cell_tag ${BARCODE_TAG} --out ${PREFIX}/fq.fq
${CURRENT_DIR}/scripts/minimap2-2.24_x64-linux/minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no $MITOFA ${PREFIX}/fq.fq > ${PREFIX}/minimap.sam
python ${CURRENT_DIR}/scripts/retag.py --sam ${PREFIX}/minimap.sam --cell_tag ${BARCODE_TAG} --out ${PREFIX}/minitagged.bam

## Sort bam file by cell barcode
samtools sort -t ${BARCODE_TAG} ${PREFIX}/minitagged.bam -o ${PREFIX}/minitagged.bam.sorted_${BARCODE_TAG}

## Remove uninformative intermediate files
rm ${PREFIX}/mt.bam
rm ${PREFIX}/fq.fq
rm ${PREFIX}/minimap.sam
rm ${PREFIX}/minitagged.bam

## Extract bam file for each cell
split_bam="${PREFIX}/minitagged.bam.sorted_${BARCODE_TAG}"
block=10000
block_num=`wc -l $BARCODE | awk -v block=$block '{printf("%0.f\n" ,($1/block))}'`
for i in `seq 0 ${block_num}`
do
	mkdir -p ${PREFIX}/bam_split_by_bc${i}
done
python ${CURRENT_DIR}/scripts/bam_split_by_bc.py $split_bam ${PREFIX}/bam_split_by_bc $BARCODE_TAG

## Variant calling for each cell
job_count=0
for i in `seq 0 ${block_num}`
do
	if [[ "ls -A ${PREFIX}/bam_split_by_bc${i}" == "" ]];then
		rmdir ${PREFIX}/bam_split_by_bc${i}
		continue
	fi
	for f in `ls ${PREFIX}/bam_split_by_bc${i}/*.bam`
	do
		bc=`basename $f .bam`
		python ${CURRENT_DIR}/scripts/pileup_counts_sc.py $f ${PREFIX}/bam_split_by_bc${i}/${bc} $MITOFA $BASEQUAL $bc $ALIGNQUAL &
		let job_count=job_count+1
		if [[ $job_count%$THREADS -eq 0 ]];then
			wait
		fi
	done
	wait
	for base in A T C G coverage
	do
		cat ${PREFIX}/bam_split_by_bc${i}/*.${base}.txt > ${PREFIX}/block${i}_all.${base}.txt &
	done
done
wait

## Generate single-cell mtRNA matrix
for base in A T C G coverage
do
	cat ${PREFIX}/block*_all.${base}.txt > ${PREFIX}/whole_all.${base}.txt &
done
wait
sh ${CURRENT_DIR}/scripts/convert_txt2mtx.sh ${PREFIX}/whole_all 
gzip ${PREFIX}/whole_all*alt/matrix.mtx
python ${CURRENT_DIR}/scripts/cal_af_sc.py ${PREFIX}/whole_all 

## Remove useless files
rm -rf ${PREFIX}/bam_split_by_bc*
rm ${PREFIX}/block*_all.*.txt
rm ${PREFIX}/whole_all.[ATCG]*txt
rm ${PREFIX}/whole_all.coverage*txt

## Remove doublets, optional
wait

sc_af="${PREFIX}/whole_all.af.txt"
if [[ "$SCRUBLET" != "FALSE" ]];then
	python ${CURRENT_DIR}/scripts/remove_doublet_by_scrublet_mtx.py $PREFIX $BLANK $SCRUBLET
	sinlist="${PREFIX}/scrublet_singlet.list"
	if [[ "$FAVG" != "FALSE" ]];then
		Rscript ${CURRENT_DIR}/scripts/remove_doublet.R $sc_af $PREFIX ${PREFIX}/scrublet_doublet.list
		sinlist="${PREFIX}/Favg_Gaussian_fitted_singlet.list"
	fi
else
	sinlist=$BARCODE
	if [[ "$FAVG" != "FALSE" ]];then
		Rscript ${CURRENT_DIR}/scripts/remove_doublet.R $sc_af $PREFIX $BLANK
		sinlist="${PREFIX}/Favg_Gaussian_fitted_singlet.list"
	fi
fi
		
## mitoSplitter demultiplexing

resdir="${PREFIX}/mitoSplitter"
python ${CURRENT_DIR}/scripts/mitoSplitter_prob.py $BULKAF $sc_af $resdir $sinlist

## mitoSplitter perform validation (available if gold standard file provided)
if [[ -n "$GOLD" ]];then
	python ${CURRENT_DIR}/scripts/mitoSplitter_prob_info.py $GOLD $resdir
fi
