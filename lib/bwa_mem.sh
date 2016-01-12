#!/bin/bash

BWA_INDEX=${RESOURCES}/${REFERENCE}/bwa/0.7.5/${REFERENCE_NAME}

function _align {
	bwa mem -t ${CPUS} -M \
		-R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA\tCN:OMICRON" \
		${BWA_INDEX} ${temp}/${sample}.r{1,2}.${FASTQ_EXT} | \
		samtools view -Sb - -o ${temp}/${sample}.aligned.bam
}

function _align_single_end {
	bwa mem -t ${CPUS} -M \
		-R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA\tCN:OMICRON" \
		${BWA_INDEX} ${temp}/${sample}.r1.${FASTQ_EXT} | \
		samtools view -Sb - -o ${temp}/${sample}.aligned.bam
}
