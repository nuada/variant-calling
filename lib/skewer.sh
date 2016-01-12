#!/bin/bash

function _trim_adapters {
	skewer --quiet \
		--threads ${CPUS} \
		-x ${ADAPTORS} \
		--compress \
		--min 50 \
		--end-quality 20 \
		--output ${temp}/${sample} \
		${DATA}/${sample}.r{1,2}.${FASTQ_EXT}
	mv ${temp}/${sample}-trimmed-pair1.fastq.gz ${temp}/${sample}.r1.fastq.gz
	mv ${temp}/${sample}-trimmed-pair2.fastq.gz ${temp}/${sample}.r2.fastq.gz
}

function _trim_adapters_single_end {
	skewer --quiet \
		--threads ${CPUS} \
		-x ${ADAPTORS} \
		--compress \
		--min 50 \
		--end-quality 20 \
		--output ${temp}/${sample} \
		${DATA}/${sample}.r1.${FASTQ_EXT}
	mv ${temp}/${sample}-trimmed.fastq.gz ${temp}/${sample}.r1.fastq.gz
}
