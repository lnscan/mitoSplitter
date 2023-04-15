#!/bin/sh

func() {
	echo "Usage:"
	echo "	sh $0 <-i input.bam> <-b barcode.list> <-o out_dir> <-r bulk_var_af.txt> [-m mito.fasta] [-t threads] [-f barcode_tag] [-q base_quality] [-a alignment_quality] [-h]"
	echo "Description:"
	echo "	mitoSpliter pipeline to multiplex barcoded single cell"
	echo "Ordering options:"
	echo "	-i	input bam file"
	echo "	-b	barcode list file, one barcode per line"
	echo "	-o	name of directory for mitoSpliter output files"
	echo "	-r	mitochondrial variant allele frequency matrix file for all samples"
	echo ""
	echo "Other options:"
	echo "	-m	mitochondrial reference genome, default = hg38_mito.fasta"
	echo "	-t	max threads to use, default = 1"
	echo "	-f	used to extract cell barcode from bam file, default = CB"
	echo "	-q	minimum base quality to be considered in bam file, used for variant calling, default = 20"
	echo "	-a	minimum alignment quality to be considered in bam file, used for variant calling, default = 20"
	echo "	-h	print this help and exit"
	exit -1
}

MITOFA="mito_fastas/hg38.fasta"
THREADS=1
BARCODE_TAG="CB"
BASEQUAL=20
ALIGNQUAL=20

while getopts 'i:b:o:r:mtfqah' OPT;do
	case $OPT in
		i) BAMFILE="$OPTARG";;
		b) BARCODE="$OPTARG";;
		o) PREFIX="$OPTARG";;
		r) BULKAF="$OPTARG";;
		m) MITOFA="$OPTARG";;
		t) THREADS="$OPTARG";;
		f) BARCODE_TAG="$OPTARG";;
		q) BASEQUAL="$OPTARG";;
		a) ALIGNQUAL="$OPTARG";;
		h) func;;
		?) func;;
	esac
done

#threads=50
#BARCODE_TAG="CB"

echo "start mitoSpliter pipeline at `date`"

samtools view -@ ${THREADS} -h $BAMFILE "MT" -b -o ${PREFIX}/mt.bam

python renamer.py --bam ${PREFIX}/mt.bam --barcodes $BARCODE --out ${PREFIX}/fq.fq
/home/songjia/bigdisk/linxr/software/minimap2-2.24_x64-linux/minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no $MITOFA ${PREFIX}/fq.fq > ${PREFIX}/minimap.sam
python retag.py --sam ${PREFIX}/minimap.sam --out ${PREFIX}/minitagged.bam
samtools sort -@ ${THREADS} ${PREFIX}/minitagged.bam -o ${PREFIX}/minitagged_sorted.bam
samtools index -@ ${THREADS} ${PREFIX}/minitagged_sorted.bam
samtools sort -t ${BARCODE_TAG} ${PREFIX}/minitagged.bam -o ${PREFIX}/minitagged.bam.sorted_${BARCODE_TAG}
split_bam="${PREFIX}/minitagged.bam.sorted_${BARCODE_TAG}"
block=10000
block_num=`wc -l $BARCODE | awk -v block=$block '{printf("%0.f\n" ,($1/block)-1)}'`
for i in `seq 0 ${block_num}`
do
	mkdir -p ${PREFIX}/bam_split_by_bc${i}
done
python bam_split_by_bc.py $split_bam ${PREFIX}/bam_split_by_bc $BARCODE_TAG

job_count=0
for i in `seq 0 ${block_num}`
do
	for f in `ls ${PREFIX}/bam_split_by_bc${i}/*.bam`
	do
		bc=`basename $f .bam`
		python 01_allbc_pileup_counts_sc.py $f ${PREFIX}/bam_split_by_bc${i}/${bc} $MITOFA $BASEQUAL $bc $ALIGNQUAL &
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

for base in A T C G coverage
do
	cat ${PREFIX}/block*_all.${base}.txt > ${PREFIX}/whole_all.${base}.txt &
done
wait
sh 02_convert_txt2mtx.sh ${PREFIX}/whole_all 
gzip ${PREFIX}/whole_q${qual}_all*alt/matrix.mtx
python 03_cal_af_var.py ${PREFIX}/whole_all 

mix_af="${PREFIX}/whole_all.var_af.txt"
resdir="${PREFIX}/mitoSpliter"
python mitoSpliter_prob.py $BULKAF $mix_af $resdir &
