#!/bin/bash

function _trim_adapters {
	TrimmomaticPE -threads ${CPUS} -phred33 \
		${DATA}/${sample}.r{1,2}.${FASTQ_EXT} \
		${temp}/${sample}.r1{,.unpaired}.${FASTQ_EXT} \
		${temp}/${sample}.r2{,.unpaired}.${FASTQ_EXT} \
		ILLUMINACLIP:${ADAPTORS}:2:25:10 \
		SLIDINGWINDOW:5:20 \
		MINLEN:100
}

function _trim_adapters_single_end {
	TrimmomaticSE -threads ${CPUS} -phred33 \
		${DATA}/${sample}.r1.${FASTQ_EXT} \
		${temp}/${sample}.r1.${FASTQ_EXT} \
		ILLUMINACLIP:${ADAPTORS}:2:25:10 \
		SLIDINGWINDOW:5:20 \
		MINLEN:100
}
