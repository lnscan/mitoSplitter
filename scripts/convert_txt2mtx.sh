#!/bin/sh

prefix=$1
CURRENT_DIR=`dirname $0`

sed "s/ /_/" ${prefix}.coverage.txt | sort -k 1,1 > ${prefix}.coverage.sorted.txt &

for altbase in A T C G
do
	mkdir -p ${prefix}.${altbase}.alt
	grep -v " $altbase$" ${prefix}.${altbase}.txt | cut -f 1-3 -d ' ' > ${prefix}.${altbase}.alt.txt
	awk '{print $1"_"$2}' ${prefix}.${altbase}.alt.txt | sort > ${prefix}.${altbase}.alt.sorted.txt &
	python ${CURRENT_DIR}/convert_txt2mtx_ATCG.py $prefix $altbase &
done
wait

for altbase in A T C G
do
	mkdir -p ${prefix}.coverage.${altbase}.alt
	join ${prefix}.coverage.sorted.txt ${prefix}.${altbase}.alt.sorted.txt | sed "s/_/ /" > ${prefix}.coverage.${altbase}.alt.txt
	python ${CURRENT_DIR}/convert_txt2mtx_ATCG.py "${prefix}.coverage" $altbase &
done
wait

