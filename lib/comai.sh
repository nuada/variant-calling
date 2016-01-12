#!/bin/bash

# Required parameters:
# * RUN - multiplexed FASTQs from samples will be demultiplexed
# Sample_id/barcode mapping should be in file: <RUN>.comai_barcodes.
# Barcode sequence length is infered from sample_id/barcode mapping file.
BARCODE_EXT=comai_barcodes
BARCODE_MISMATCHES=2

barcode_stats_usage='<run_id>'
function barcode_stats {
	log 'Barcode statistics'
	barcode_length=$(echo -n $(head -n 1 ${DATA}/${1}.${BARCODE_EXT} | cut -f 2) | wc --chars)
	for i in $(seq 2); do
		${FASTQ_CAT} ${DATA}/${1}.r${i}.${FASTQ_EXT} | \
			awk 'BEGIN { n=1 } { if (n == 2) print $0; if (n == 4) n=1; else n++; }' | \
			cut -c -${barcode_length} | sort | uniq -c | \
			sort -n -r > ${RESULTS}/${1}.r${i}.prefix_stats
		sort -k 2 ${RESULTS}/${1}.r${i}.prefix_stats | \
			join -j 2 <(sort -k 2 ${DATA}/${1}.${BARCODE_EXT}) - | \
			sort -k 3 -n -r > ${RESULTS}/${1}.r${i}.barcode_stats
	done
}

demultiplex_usage='<run_id>'
function demultiplex {
	log 'Demultiplex'
	fastq-multx -b -v ' ' -m ${BARCODE_MISMATCHES} -B ${DATA}/${1}.${BARCODE_EXT} \
		${DATA}/${1}.r{1,2}.${FASTQ_EXT} \
		-o ${DATA}/%.${1}.r1.${FASTQ_EXT} \
		-o ${DATA}/%.${1}.r2.${FASTQ_EXT}
}
