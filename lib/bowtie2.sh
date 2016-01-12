#!/bin/bash

BOWTIE2_INDEX=${RESOURCES}/${REFERENCE}/bowtie2/${REFERENCE_NAME}

function _align {
	bowtie2 -p ${CPUS} -x ${BOWTIE2_INDEX} -1 ${temp}/${sample}.r1.${FASTQ_EXT} -2 ${temp}/${sample}.r2.${FASTQ_EXT} \
		--rg-id "${sample}" \
		--rg "SM:${sample}" \
		--rg "LB:${sample}" \
		--rg "PL:ILLUMINA" \
		--rg "CN:OMICRON" | \
		samtools view -Sb - -o ${temp}/${sample}.aligned.bam
}

function _align_single_end {
	echo 'Single-end alignment not supported!'
	exit 1
}
