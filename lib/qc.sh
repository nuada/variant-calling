#!/bin/bash

function check_fastq_pair {
	if [[ ! -e ${DATA}/${1}.r2.${FASTQ_EXT} ]]; then
		log 'Reverse read not found!'
		return 0
	fi

	read1=$(${FASTQ_CAT} ${DATA}/${1}.r1.${FASTQ_EXT} | head -n 1 | cut -d ' ' -f 1)
	read2=$(${FASTQ_CAT} ${DATA}/${1}.r2.${FASTQ_EXT} | head -n 1 | cut -d ' ' -f 1)
	if [[ ${read1} != ${read2} ]]; then
		log 'PAIRED READS DO NOT MATCH!'
		exit 1
	fi
}

# TODO should run in parallel!
qc_usage='<sample_id>'
function qc {
	check_fastq_pair $*

	for i in $(seq 2); do
		if [[ -e ${DATA}/${1}.r${i}.${FASTQ_EXT} ]]; then
			fastqc -t ${CPUS} -o ${RESULTS} ${DATA}/${1}.r${i}.${FASTQ_EXT}
		fi
	done
}

function doc_hist {
	coverageBed -hist -abam ${RESULTS}/${1}.bam -b ${INTERVALS} | \
		grep ^all | sed -e "s/^/${1}\t/g" > ${RESULTS}/${1}.hist.doc
}

doc_fractions_usage=''
function doc_fractions {
	for i in ${SAMPLES}; do
		if [[ ! -e ${RESULTS}/${i}.hist.doc ]]; then
			batch "doc_hist-$i" doc_hist $i
		fi
	done
	wait_for_jobs
	
	for run in ${RUN}; do
		echo 'sample feature count bases length fraction' | tr ' ' '\t' > ${temp}/doc_fractions_${run}
		cat ${RESULTS}/*.${run}.hist.doc >> ${temp}/doc_fractions_${run}
		Rscript ../scripts/doc_fractions.R ${temp}/doc_fractions_${run}
	done
	cp -v ${temp}/doc_fractions* ${RESULTS}
}
