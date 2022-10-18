#!/bin/bash

if (( $# < 3 )); then
	echo "usage: ${0##*/} <workspace> <pair_set> <file_name_tag>"
	exit 1
fi

workspace=$1
pair_set=$2
tag=$3

#fiss annot_get $workspace pair \
#	pset=$pair_set \
#	snowman_breakpoints_capture \
#	snowman_germline_indel_vcf_capture \
#	snowman_somatic_indel_vcf_capture \
#	snowman_germline_vcf_capture \
#	snowman_somatic_vcf_capture \
#	snowman_assembled_contigs_capture \
#	> ${tag}_snowman.tsv

#fiss annot_get $workspace pair \
#	pset=$pair_set \
#	union_maf_file_forcecalled \
#	maf_file_capture_forcecalls \
#	indel_maf_file_capture_forcecalled \
#	> ${tag}_union.tsv

#fiss annot_get $workspace pair \
#	pset=$pair_set \
#	delly_vcf \
#	lumpy_vcf \
#	platypus_wgs_all \
#	platypus_somatic_vcf \
#	platypus_somatic_indel_wgs \
#	> ${tag}_sv.tsv

fiss annot_get $workspace pair \
	pset=$pair_set \
	svaba_breakpoints_capture \
	svaba_germline_indel_vcf_capture \
	svaba_somatic_indel_vcf_capture \
	svaba_germline_vcf_capture \
	svaba_somatic_vcf_capture \
	svaba_contigs_capture \
	> ${tag}_svaba.tsv



