#!/bin/bash

PICARD=/usr/bin/picard-tools

function bed_to_interval_list {
	echo '@HD VN:1.0 SO:coordinate' | tr " " "\t"
	samtools view -H $1 | grep '@SQ'
	awk -F "\t" '{ print $1 FS $2+1 FS $3 FS "+" FS $4 }' $2
}

function hsmetrics_queue {
	bed_to_interval_list ${RESULTS}/${3}.bam $1 > ${temp}/bait
	bed_to_interval_list ${RESULTS}/${3}.bam $2 > ${temp}/target

	${PICARD} CalculateHsMetrics \
		REFERENCE_SEQUENCE=${HG} \
		BAIT_INTERVALS=${temp}/bait \
		TARGET_INTERVALS=${temp}/target \
		INPUT=${RESULTS}/${3}.bam \
		OUTPUT=${temp}/${3}.picard_hsmetrics \
		PER_TARGET_COVERAGE=${temp}/${3}.picard_per_target_hsmetrics

	tail -n 4 ${temp}/${3}.picard_hsmetrics | head -n 2 | sed -e "s/^/${3}\t/" > ${RESULTS}/${3}.picard_hsmetrics
	cat ${temp}/${3}.picard_per_target_hsmetrics | sed -e "s/^/${3}\t/" > ${RESULTS}/${3}.picard_per_target_hsmetrics
}

hsmetrics_usage='<bait file> <target file>'
function hsmetrics {
	for i in ${SAMPLES}; do
		batch "hsmetrics-aligned-$i" hsmetrics_queue $1 $2 ${i}.aligned
		batch "hsmetrics-$i" hsmetrics_queue $1 $2 $i
	done
	wait_for_jobs

	for j in aligned.picard picard; do
		head -n 1 ${RESULTS}/$(set -- ${SAMPLES}; echo $1).${j}_hsmetrics | cut -f 2- | sed -e "s/^/SAMPLE_ID\t/" > ${RESULTS}/${j}_hsmetrics.csv
		head -n 1 ${RESULTS}/$(set -- ${SAMPLES}; echo $1).${j}_per_target_hsmetrics | cut -f 2- | sed -e "s/^/SAMPLE_ID\t/" > ${RESULTS}/${j}_per_target_hsmetrics.csv
		for i in ${SAMPLES}; do
			tail -n 1 ${RESULTS}/${i}.${j}_hsmetrics >> ${RESULTS}/${j}_hsmetrics.csv
			tail -n +2 ${RESULTS}/${i}.${j}_per_target_hsmetrics >> ${RESULTS}/${j}_per_target_hsmetrics.csv
		done
	done

	rm -rf ${RESULTS}/*.picard_*hsmetrics #hsmetrics-*.o*
}

