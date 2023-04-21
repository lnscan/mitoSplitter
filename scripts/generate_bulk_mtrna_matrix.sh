#!/bin/sh

## Command parameters
func() {
    echo "Usage:"
    echo "  sh $0 <-l bulk_bam.list> <-o out_dir> [-m mito.fasta] [-t threads] [-q base_quality] [-a alignment_quality] [-h]"
    echo "Description:"
    echo "  mitoSplitter pipeline to generate bulk mtRNA matrix"
    echo "Ordering options:"
    echo "  -l  RNA-seq generated bam list file, one bam file per line"
	echo "  -o  name of directory for mitoSplitter output files"
    echo ""
    echo "Other options:"
    echo "  -m  mitochondrial reference genome fasta, default = GRCh38_MT.fasta"
    echo "  -t  max threads to use, default = 50"
    echo "  -q  minimum base quality to be considered in bam file, used for variant calling, default = 10"
    echo "  -a  minimum alignment quality to be considered in bam file, used for variant calling, default = 10"
    echo "  -h  print this help and exit"
	echo ""
    exit -1
}

CURRENT_DIR=`dirname $0`
MITOFA="${CURRENT_DIR}/../mito_fastas/GRCh38_MT.fasta"
PREFIX="."
THREADS=50
BASEQUAL=10
ALIGNQUAL=10
BULKLIST=""

while getopts 'l:o:m:t:q:a:h' OPT;do
    case $OPT in
        l) BULKLIST="$OPTARG";;
        o) PREFIX="$OPTARG";;
        m) MITOFA="$OPTARG";;
        t) THREADS="$OPTARG";;
        q) BASEQUAL="$OPTARG";;
        a) ALIGNQUAL="$OPTARG";;
        h) func;;
        ?) func;;
    esac
done

if [[ -z "$BULKLIST" ]];then
	echo "ERROR : No bam list file!"
	echo ""
	func;
fi

while read bamf
do
	echo "Start generate mtRNA matrix for $bamf at `date`"
	## Make output dir
	sam=`basename $bamf .bam`
	resdir="${PREFIX}/${sam}"
	mkdir -p $resdir

	## Remapping
	python ${CURRENT_DIR}/extract_bulk_fq.py --bam ${bamf} --out ${resdir}/${sam}.fq
	${CURRENT_DIR}/minimap2-2.24_x64-linux/minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no $MITOFA ${resdir}/${sam}.fq > ${resdir}/minimap_${sam}.sam
	samtools view -@ $THREADS -b ${resdir}/minimap_${sam}.sam > ${resdir}/minitagged_${sam}.bam

	## Remove uninformative intermediate files
	rm ${resdir}/${sam}.fq
	rm ${resdir}/minimap_${sam}.sam

	## Variant calling
	python ${CURRENT_DIR}/pileup_counts_bulk.py ${resdir}/minitagged_${sam}.bam ${resdir}/${sam} $MITOFA $BASEQUAL $sam $ALIGNQUAL 
	python ${CURRENT_DIR}/cal_af_bulk.py ${resdir}/${sam}
	echo "${resdir}/${sam}.af.txt" >> ${PREFIX}/bulk_af.list

done < $BULKLIST
python ${CURRENT_DIR}/merge_af_bulk.py ${PREFIX}/bulk_af.list $MITOFA $resdir


