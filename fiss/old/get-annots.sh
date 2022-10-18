#!/bin/bash

if (( $# < 3 )); then
	echo "usage: ${0##*/} <workspace> <pair_set> <file_name_tag>"
	exit 1
fi

workspace=$1
pair_set=$2
tag=$3

fiss annot_get $workspace samples \
	pset=$pair_set \
	clean_bam_file_capture \
	> brain-mets_${tag}_bams.tsv

sed 1d brain-mets_${tag}_bams.tsv \
	| cut -f 2 \
	| sed '/^\s*$/d' \
	> brain-mets_${tag}_bams.txt


fiss annot_get $workspace pair \
	pset=$pair_set \
	case_sample control_sample \
	> brain-mets_${tag}_pairs.tsv


# SNV calls
fiss annot_get $workspace pair \
	pset=$pair_set \
	somatic_mutation_coverage_capture \
	call_stats_capture \
	maf_file_oxoG3_capture \
	maf_file_ffpeBias_capture \
	call_stats_wex_pairs_novo_realign \
	maf_file_capture_realign_pairs_novo \
	maf_file_capture_master_filter_removed \
	> brain-mets_${tag}_snv.tsv

# INDEL calls
fiss annot_get $workspace pair \
	pset=$pair_set \
	strelka_passed_somatic_indel_maf_file_capture_pair \
	> brain-mets_${tag}_indel.tsv

# final MAF files
#fiss annot_get $workspace pair \
#	pset=$pair_set \
#	maf_file_capture_master_filter_removed \
#	maf_file_capture_indel_master_filter_removed \
#	union_maf_file_forcecalled \
#	> brain-mets_${tag}_maf.tsv

fiss annot_get $workspace samples \
	pset=$pair_set \
	contamination_percentage_consensus_capture \
	tumor_subtype \
	source_subtype_capture \
	processed_subtype_capture \
	recapseg_pre_qc \
	recapseg_post_qc \
	recapseg_num_seg \
	recapseg_segment_cn_sum \
	> brain-mets_${tag}_sample-stats.tsv


fiss annot_get $workspace pair \
	pset=$pair_set \
	somatic_mutation_covered_bases_capture \
	combined_deTiN_TiN_num \
	picard_oxoQ \
	oxoG3_reject_count_capture \
	oxoG3_pass_count_capture \
	ffpe_Q \
	ffpe_reject_count_capture \
	ffpe_pass_count_capture \
	> brain-mets_${tag}_pair-stats.tsv

fiss annot_get $workspace pair \
	pset=$pair_set \
	snowman_breakpoints_capture \
	snowman_germline_indel_vcf_capture \
	snowman_somatic_indel_vcf_capture \
	snowman_germline_vcf_capture \
	snowman_somatic_vcf_capture \
	snowman_assembled_contigs_capture \
	> brain-mets_${tag}_snowman.tsv

