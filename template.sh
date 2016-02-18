#!/bin/bash
set -e

LIBDIR=lib
source ${LIBDIR}/reference_hg19.sh

PROJECT=/data/template
INTERVALS=${PROJECT}/intervals.${REFERENCE}.bed
SAMPLES=$(cd ${PROJECT}; ls -1 s??*.fastq.gz | cut -d . -f 1-2 | sort -u)
RUN=

source ${LIBDIR}/common.sh
source ${LIBDIR}/qc.sh
source ${LIBDIR}/picard.sh
source ${LIBDIR}/skewer.sh
source ${LIBDIR}/bwa_mem.sh
source ${LIBDIR}/snpeff_3x.sh
source ${LIBDIR}/gatk.sh
source ${LIBDIR}/xhmm.sh
source ${LIBDIR}/gemini.sh

CPUS=12
CLUSTER_NODES=2

ADAPTORS=${RESOURCES}/qc/TruSeq3-PE-2.fa
SCRIPTS=scripts

queue_usage=''
function queue {
	for i in ${SAMPLES}; do
		batch "align-$i" align $i
		batch "qc-$i" qc $i
	done
	wait_for_jobs

    CPUS=2
	hsmetrics ${INTERVALS} ${INTERVALS}

    CPUS=12
	for i in $(seq 1 ${CLUSTER_NODES}); do
		batch "variant_calling-$i" variant_calling_chunk $i ${CLUSTER_NODES}
	done
	wait_for_jobs
	variant_calling_merge_chunks

	variant_stats | tee ${RESULTS}/variant_stats.log

	gemini_load
	variant_filtering
}

variant_filtering_usage=''
function variant_filtering {
	${SCRIPTS}/gemini_summarize.py 'select variant_id, chrom, start, end, ref, type,
			sub_type, aaf, num_hom_ref, num_het, num_hom_alt, num_unknown, hwe,
			gene, exon, transcript, codon_change, aa_change, aa_length, pfam_domain,
			biotype, impact_so, impact_severity, in_dbsnp,
			rs_ids, in_esp, aaf_esp_ea, aaf_esp_all, in_1kg, aaf_1kg_eur, aaf_1kg_all,
			in_exac, aaf_exac_all, aaf_adj_exac_all, aaf_adj_exac_nfe,
			in_omim, rmsk, in_cpg_island, in_segdup, cadd_raw, cadd_scaled, fitcons,
			clinvar_sig, clinvar_disease_name, clinvar_dbsource, clinvar_dbsource_id,
			clinvar_origin, clinvar_dsdb, clinvar_dsdbid, clinvar_disease_acc
		from variants
		where filter is NULL' \
		${RESULTS}/variants.db > ${RESULTS}/variants-all.csv
}

main $*
