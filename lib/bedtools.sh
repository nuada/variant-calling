#!/bin/bash

# TODO function for bedtools protocols!

function offtarget {
	sample=$1
	bedtools slop -i ${INTERVALS} -g ${RESOURCES}/${REFERENCE}/$(echo ${REFERENCE_NAME} | sed -e 's/fasta/genome/') -b 150 | \
		bedtools sort -i - | \
		bedtools merge -i - -d 150 | \
		bedtools sort -i - > ${temp}/offtarget.intervals.bed
	gatk -T PrintReads --num_cpu_threads_per_data_thread ${CPUS} \
		--reference_sequence ${HG} \
		--input_file ${RESULTS}/${sample}.bam \
		-o ${sample}.offtarget.bam \
		--excludeIntervals ${temp}/offtarget.intervals.bed
}
