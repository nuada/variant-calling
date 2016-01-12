#!/bin/bash

function xhmm {
	${RESOURCES}/software/xhmm/xhmm "$@"
}

function pseq {
	${RESOURCES}/software/plinkseq-0.08/pseq "$@"
}

function xhmm_doc {
	sample=$1
	gatk -T DepthOfCoverage \
		--reference_sequence ${HG} \
		--input_file ${RESULTS}/${sample}.bam \
		-o ${temp}/${sample}.doc \
		--intervals ${INTERVALS} \
		-dt BY_SAMPLE -dcov 5000 -l INFO \
		--omitDepthOutputAtEachBase \
		--omitLocusTable \
		--minBaseQuality 0 \
		--minMappingQuality 20 \
		--start 1 --stop 5000 --nBins 200 \
		--includeRefNSites \
		--countType COUNT_FRAGMENTS

	cp -v ${temp}/${sample}.doc* ${RESULTS}/
}

cnv_usage=''
function cnv {
	for i in ${SAMPLES}; do
		batch "xhmm_doc-$i" xhmm_doc $i &
	done
	wait
#	rm -f xhmm_doc-*.o*

	# Merge per sample DOC files
	xhmm --mergeGATKdepths -o ${temp}/DATA.RD.txt $(for sample in ${SAMPLES}; do echo --GATKdepths=${RESULTS}/${sample}.doc.sample_interval_summary; done)

	# Targets with extreme GC content
	gatk -T GCContentByInterval \
		--reference_sequence ${HG} \
		-o ${temp}/DATA.locus_GC.txt \
		--intervals ${INTERVALS}
	awk '$2 < 0.1 || $2 > 0.9' ${temp}/DATA.locus_GC.txt | cut -f1 > ${temp}/extreme_gc_targets.txt

	# Calculate the fraction of repeat-masked bases in each target
	cat ${INTERVALS} | (echo -e "#CHR\tBP1\tBP2\tID"; awk -F '\t' '{ print $1 FS $2+1 FS $3 FS NR }') > ${temp}/EXOME.targets.reg
	pushd ${temp}
	pseq . loc-load  --locdb ${temp}/EXOME.targets.LOCDB --file ${temp}/EXOME.targets.reg --group targets
	pseq . loc-stats --locdb ${temp}/EXOME.targets.LOCDB --group targets --seqdb /resources/hg19/plinkseq/seqdb | \
		tail -n +2 | sort -k1 -g | awk -F '\t' '{print $4 FS $10}' | \
		sed -e 's/\.\./-/g' > ${temp}/DATA.locus_complexity.txt
	popd
	awk '$2 > 0.25' ${temp}/DATA.locus_complexity.txt | cut -f1 > ${temp}/low_complexity_targets.txt

	xhmm --matrix -r ${temp}/DATA.RD.txt --centerData --centerType target \
		-o ${temp}/DATA.filtered_centered.RD.txt \
		--outputExcludedTargets ${temp}/DATA.filtered_centered.RD.txt.filtered_targets.txt \
		--outputExcludedSamples ${temp}/DATA.filtered_centered.RD.txt.filtered_samples.txt \
		--excludeTargets ${temp}/extreme_gc_targets.txt \
		--excludeTargets ${temp}/low_complexity_targets.txt \
		--minTargetSize 10   --maxTargetSize 10000 \
		--minMeanTargetRD 10 --maxMeanTargetRD 500 \
		--minMeanSampleRD 25 --maxMeanSampleRD 200 \
		--maxSdSampleRD 150

	# PCA on mean-centered data
	xhmm --PCA -r ${temp}/DATA.filtered_centered.RD.txt --PCAfiles ${temp}/DATA.RD_PCA

	# Normalizes data using PCA information
	xhmm --normalize -r ${temp}/DATA.filtered_centered.RD.txt --PCAfiles ${temp}/DATA.RD_PCA \
		--normalizeOutput ${temp}/DATA.PCA_normalized.txt \
		--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

	# Filter and z-score center (by sample) the PCA-normalized data
	xhmm --matrix -r ${temp}/DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
		-o ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
		--outputExcludedTargets ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
		--outputExcludedSamples ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
		--maxSdTargetRD 30

	# Filter original read-depth data to be the same as filtered, normalized data
	xhmm --matrix -r ${temp}/DATA.RD.txt \
	--excludeTargets ${temp}/DATA.filtered_centered.RD.txt.filtered_targets.txt \
	--excludeTargets ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
	--excludeSamples ${temp}/DATA.filtered_centered.RD.txt.filtered_samples.txt \
	--excludeSamples ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
	-o ${temp}/DATA.same_filtered.RD.txt

	# Discover CNVs in normalized data
	echo "1e-08	6	70	-3	1	0	1	3	1" > ${temp}/params.txt
	xhmm --discover -p ${temp}/params.txt \
		-r ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ${temp}/DATA.same_filtered.RD.txt \
		-c ${temp}/DATA.xcnv -a ${temp}/DATA.aux_xcnv -s ${temp}/DATA

	# Genotype discovered CNVs in all samples
	xhmm --genotype -p ${temp}/params.txt \
		-r ${temp}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ${temp}/DATA.same_filtered.RD.txt \
		-g ${temp}/DATA.xcnv -F ${HG} \
		-v ${temp}/DATA.vcf
	cp -rv ${temp}/DATA.vcf ${RESULTS}/xhmm.vcf &

	# Annotate exome targets with their corresponding genes
	cat ${INTERVALS} | awk -F '\t' '{ split($4, name, "."); print $1 ":" $2+1 ".." $3 FS "1" FS name[1] }' > ${temp}/annotated_targets.refseq.loci

	mkdir ${temp}/xhmm_plots
	echo "library(gplots)
library(plotrix)
source('${RESOURCES}/software/xhmm/scripts/R_functions/sourceDir.R')
sourceDir('${RESOURCES}/software/xhmm/scripts/R_functions',trace = FALSE, recursive = TRUE)
setwd('${temp}')
XHMM_plots('xhmm_plots', 'DATA', 'annotated_targets.refseq.loci')" | Rscript -

	pushd ${temp}
	zip -r ${RESULTS}/xhmm_plots.zip xhmm_plots
	popd

	wait
}
