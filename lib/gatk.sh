#!/bin/bash

# Parameters
GATK_FILTERING_METHOD=_hard_filtering

GATK_HARD_FILTERING_SNPS='QDFilter=QD<2.0
	MQFilter=MQ<40.0
	FSFilter=FS>60.0
	HaplotypeScoreFilter=HaplotypeScore>13.0
	MQRankSumFilter=MQRankSum<-12.5
	ReadPosFilter=ReadPosRankSum<-8.0'

GATK_HARD_FILTERING_INDELS='QDFilter=QD<2.0
	ReadPosFilter=ReadPosRankSum<-20.0
	FSFilter=FS>200.0'

GATK_VARIANTS_PREFIX=variants

GATK_RESOURCES=${RESOURCES}/${REFERENCE}/GATK
DBSNP=${GATK_RESOURCES}/dbsnp_13?.${REFERENCE}.vcf
KG_OMNI=${GATK_RESOURCES}/1000G_omni2.5.${REFERENCE}.vcf
KG_INDELS=${GATK_RESOURCES}/1000G_phase1.indels.${REFERENCE}.vcf
KG_SNPS=${GATK_RESOURCES}/1000G_phase1.snps.high_confidence.${REFERENCE}.vcf
HAPMAP=${GATK_RESOURCES}/hapmap_3.3.${REFERENCE}.vcf
MILLS_N_KG_INDELS=${GATK_RESOURCES}/Mills_and_1000G_gold_standard.indels.${REFERENCE}.vcf

function gatk {
	java ${JAVA_OPTS} -jar ${RESOURCES}/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar "$@"
}

align_usage='<sample_id>'
function align {
	sample=$1

    # Test if file already exists
	if [[ -f ${RESULTS}/${sample}.bam ]]; then
        exit 0
    fi

	log 'Alignment'
	if [[ -e ${DATA}/${sample}.r2.${FASTQ_EXT} ]]; then
		_trim_adapters
		_align
	else
		_trim_adapters_single_end
		_align_single_end
	fi

	log 'Sort BAM'
	samtools sort ${temp}/${sample}.aligned.bam ${temp}/${sample}.sorted
	mv ${temp}/${sample}.sorted.bam ${temp}/${sample}.aligned.bam
	samtools index ${temp}/${sample}.aligned.bam
	samtools flagstat ${temp}/${sample}.aligned.bam > ${temp}/${sample}.aligned.flagstat
	cp -v ${temp}/${sample}.aligned.* ${RESULTS}/ &

	log 'Drop reads with low MQ'
	samtools view -b -q 10 ${temp}/${sample}.aligned.bam -o ${temp}/${sample}.filtered.bam
	samtools index ${temp}/${sample}.filtered.bam

	log 'Mapping stats'
	samtools flagstat ${temp}/${sample}.filtered.bam > ${temp}/${sample}.flagstat
	cp -v ${temp}/${sample}.flagstat ${RESULTS}/ &

	log 'Realignment'
	gatk -T RealignerTargetCreator --num_threads ${CPUS} \
		--reference_sequence ${HG} \
		--input_file ${temp}/${sample}.filtered.bam \
		-o ${temp}/${sample}.realign.intervals \
		--intervals ${INTERVALS} \
		--known ${KG_INDELS} \
		--known ${MILLS_N_KG_INDELS}

	gatk -T IndelRealigner \
		--reference_sequence ${HG} \
		--input_file ${temp}/${sample}.filtered.bam \
		--targetIntervals ${temp}/${sample}.realign.intervals \
		-o ${temp}/${sample}.realigned.bam \
		--intervals ${INTERVALS} \
		--knownAlleles ${KG_INDELS} \
		--knownAlleles ${MILLS_N_KG_INDELS}

	log 'Recalibration'
	gatk -T BaseRecalibrator --num_cpu_threads_per_data_thread ${CPUS} \
		--reference_sequence ${HG} \
		--input_file ${temp}/${sample}.realigned.bam \
		-o ${temp}/${sample}.recal.table \
		--intervals ${INTERVALS} \
		--knownSites ${KG_INDELS} \
		--knownSites ${MILLS_N_KG_INDELS}

	gatk -T PrintReads --num_cpu_threads_per_data_thread ${CPUS} \
		--reference_sequence ${HG} \
		--input_file ${temp}/${sample}.realigned.bam \
		--BQSR ${temp}/${sample}.recal.table \
		-o ${temp}/${sample}.recal.bam \
		--intervals ${INTERVALS}
	samtools index ${temp}/${sample}.recal.bam

	log 'Remove duplicates'
	samtools rmdup ${temp}/${sample}.recal.bam ${temp}/${sample}.bam
	samtools index ${temp}/${sample}.bam
	cp -v ${temp}/${sample}.bam* ${RESULTS}/

	wait
}

function _hard_filtering {
	log 'SNPs hard filtering'
	gatk -T VariantFiltration \
		--no_cmdline_in_header \
		--reference_sequence ${HG} \
		--variant ${temp}/${GATK_VARIANTS_PREFIX}.snps.vcf \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.snps.filtered.vcf \
		--clusterWindowSize 10 \
		$(for filter in ${GATK_HARD_FILTERING_SNPS}; do \
			name=$(echo $filter | cut -d = -f 1); \
			expression=$(echo $filter | cut -d = -f 2); \
			echo "--filterExpression \"${expression}\" --filterName \"${name}\""; \
		done) \
		--mask ${temp}/${GATK_VARIANTS_PREFIX}.indels.vcf \
		--maskName InDel

	log 'Indels hard filtering'
	gatk -T VariantFiltration \
		--no_cmdline_in_header \
		--reference_sequence ${HG} \
		--variant ${temp}/${GATK_VARIANTS_PREFIX}.indels.vcf \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.indels.filtered.vcf \
		--clusterWindowSize 10 \
		$(for filter in ${GATK_HARD_FILTERING_INDELS}; do \
			name=$(echo $filter | cut -d = -f 1); \
			expression=$(echo $filter | cut -d = -f 2); \
			echo "--filterExpression \"${expression}\" --filterName \"${name}\""; \
		done) \
	# TODO --maskName?!?
}

function _vqsr {
	log 'VQSR'
	gatk -T VariantRecalibrator --num_threads ${CPUS} \
		--reference_sequence ${HG} \
		--input ${temp}/${GATK_VARIANTS_PREFIX}.snps.vcf \
		--recal_file ${temp}/${GATK_VARIANTS_PREFIX}.snps.recal \
		--tranches_file ${temp}/${GATK_VARIANTS_PREFIX}.snps.tranches \
		--rscript_file ${temp}/${GATK_VARIANTS_PREFIX}.snps.plots.R \
		--maxGaussians 4 \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
		--resource:omni,known=false,training=true,truth=true,prior=12.0   ${KG_OMNI} \
		--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${KG_SNPS} \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0  ${DBSNP} \
		--use_annotation QD \
		--use_annotation MQRankSum \
		--use_annotation ReadPosRankSum \
		--use_annotation FS \
		--use_annotation DP \
		--mode SNP

	gatk -T VariantRecalibrator --num_threads ${CPUS} \
		--reference_sequence ${HG} \
		--input ${temp}/${GATK_VARIANTS_PREFIX}.indels.vcf \
		--recal_file ${temp}/${GATK_VARIANTS_PREFIX}.indels.recal \
		--tranches_file ${temp}/${GATK_VARIANTS_PREFIX}.indels.tranches \
		--rscript_file ${temp}/${GATK_VARIANTS_PREFIX}.indels.plots.R \
		--maxGaussians 4 \
		--resource:mills,known=false,training=true,truth=true,prior=12.0 ${MILLS_N_KG_INDELS} \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
		--use_annotation DP \
		--use_annotation FS \
		--use_annotation ReadPosRankSum \
		--use_annotation MQRankSum \
		--mode INDEL
	cp -v ${temp}/${GATK_VARIANTS_PREFIX}.*.plots.R ${RESULTS}/ &

	for variant in SNP INDEL; do
		var=$(echo ${variant} | tr '[:upper:]' '[:lower:]')
		gatk -T ApplyRecalibration --num_threads ${CPUS} \
			--reference_sequence ${HG} \
			--input ${temp}/${GATK_VARIANTS_PREFIX}.raw_annotated.vcf \
			--recal_file ${temp}/${GATK_VARIANTS_PREFIX}.${var}s.recal \
			--tranches_file ${temp}/${GATK_VARIANTS_PREFIX}.${var}s.tranches \
			--ts_filter_level 99.0 \
			--mode ${variant} \
			-o ${temp}/${GATK_VARIANTS_PREFIX}.${var}s.filtered.vcf
	done

	wait
}

function variant_calling_chunk {
	chunk=$1
	parts=$2
	split -n l/${chunk}/${parts} ${INTERVALS} > ${temp}/chunk${chunk}-intervals.bed
	INTERVALS=${temp}/chunk${chunk}-intervals.bed
	GATK_VARIANTS_PREFIX=${GATK_VARIANTS_PREFIX}-chunk${chunk}
	variant_calling
}

function variant_calling_merge_chunks {
	for i in filtered raw_annotated; do
		gatk -T CombineVariants --num_threads ${CPUS} \
			--reference_sequence ${HG} \
			$(for file in ${RESULTS}/${GATK_VARIANTS_PREFIX}-chunk*.${i}.vcf; do echo --variant ${file}; done) \
			-o ${temp}/${GATK_VARIANTS_PREFIX}.${i}.vcf
		cp -v ${temp}/${GATK_VARIANTS_PREFIX}.${i}.vcf* ${RESULTS}/ &
	done

	wait
}

variant_calling_usage=''
function variant_calling {
	log 'Raw variant calls'
	gatk -T UnifiedGenotyper --num_threads ${CPUS} \
		--no_cmdline_in_header \
		--reference_sequence ${HG} --dbsnp ${DBSNP} \
		$(for sample in ${SAMPLES}; do echo --input_file ${RESULTS}/${sample}.bam; done) \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.raw.vcf \
		--annotation Coverage \
		--group Standard \
		--standard_min_confidence_threshold_for_calling 30.0 \
		--standard_min_confidence_threshold_for_emitting 10.0 \
		--max_alternate_alleles 30 \
		--genotype_likelihoods_model BOTH \
		--intervals ${INTERVALS}

	log 'Variant annotation'
	pushd ${temp}
	snpeff eff -t -v -i vcf -o gatk ${REFERENCE_SNPEFF} ${temp}/${GATK_VARIANTS_PREFIX}.raw.vcf > ${temp}/${GATK_VARIANTS_PREFIX}.snpeff.vcf
	popd

	gatk -T VariantAnnotator --num_threads ${CPUS} \
		--no_cmdline_in_header \
		--reference_sequence ${HG} \
		--annotation SnpEff \
		--variant ${temp}/${GATK_VARIANTS_PREFIX}.raw.vcf \
		--snpEffFile ${temp}/${GATK_VARIANTS_PREFIX}.snpeff.vcf \
		--intervals ${temp}/${GATK_VARIANTS_PREFIX}.raw.vcf \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.raw_annotated.vcf
	cp -v ${temp}/${GATK_VARIANTS_PREFIX}.raw_annotated.vcf* ${RESULTS}/ &

	log 'Separate variants'
	for variant in SNP INDEL; do
		gatk -T SelectVariants --num_threads ${CPUS} \
			--reference_sequence ${HG} \
			--variant ${temp}/${GATK_VARIANTS_PREFIX}.raw_annotated.vcf \
			-o ${temp}/${GATK_VARIANTS_PREFIX}.$(echo ${variant} | tr '[:upper:]' '[:lower:]')s.vcf \
			--selectTypeToInclude ${variant}
	done

	${GATK_FILTERING_METHOD}

	log 'Combine filtered variants'
	gatk -T CombineVariants --num_threads ${CPUS} \
		--genotypemergeoption PRIORITIZE \
		--rod_priority_list snps,indels \
		--reference_sequence ${HG} \
		--variant:snps ${temp}/${GATK_VARIANTS_PREFIX}.snps.filtered.vcf \
		--variant:indels ${temp}/${GATK_VARIANTS_PREFIX}.indels.filtered.vcf \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.filtered.vcf
	cp -v ${temp}/${GATK_VARIANTS_PREFIX}.filtered.vcf* ${RESULTS}/

	wait
}

variant_stats_usage=''
function variant_stats {
	log 'Variant evaluation'
	gatk -T VariantEval --num_threads ${CPUS} \
		--reference_sequence ${HG} --dbsnp ${DBSNP} \
		--eval ${RESULTS}/${GATK_VARIANTS_PREFIX}.raw_annotated.vcf \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.eval.gatkreport \
		--stratificationModule FunctionalClass \
		--stratificationModule Filter

	log 'Filtered variants evaluation'
	gatk -T VariantEval --num_threads ${CPUS} \
		--reference_sequence ${HG} --dbsnp ${DBSNP} \
		--eval ${RESULTS}/${GATK_VARIANTS_PREFIX}.filtered.vcf \
		-o ${temp}/${GATK_VARIANTS_PREFIX}.filtered.eval.gatkreport \
		--stratificationModule FunctionalClass \
		--stratificationModule Filter

	cp -v ${temp}/${GATK_VARIANTS_PREFIX}.*.gatkreport ${RESULTS}/
}
