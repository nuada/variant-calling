#!/bin/bash

GEMINI_VARIANTS=${RESULTS}/variants.filtered.vcf
GEMINI_DB=${RESULTS}/variants.db
GEMINI_PED=

gemini_load_usage=''
function gemini_load {
	cat ${GEMINI_VARIANTS} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - -o ${temp}/decomposed.vcf
	vt normalize -r ${RESOURCES}/${REFERENCE}/${REFERENCE_NAME} ${temp}/decomposed.vcf -o ${temp}/normalized.vcf
	snpeff eff -hgvs -t -v ${REFERENCE_SNPEFF} ${temp}/normalized.vcf > ${temp}/annotated.snpeff.vcf

	if [[ ${GEMINI_PED} ]]; then
		ped="-p ${GEMINI_PED}"
	fi
	gemini load -v ${temp}/annotated.snpeff.vcf --cores ${CPUS} -t snpEff ${temp}/gemini.db ${ped}
	cp -v ${temp}/gemini.db ${GEMINI_DB}
}

gemini_query_usage=''
function gemini_query {
	gemini query --header --show-samples -q "$1" ${GEMINI_DB}
}
