#!/bin/bash

inseg=brain-mets_pass_luad_absolute.lgr.seg
#inseg=brain-mets_pass_luad_absolute_bmet-only.lgr.seg
#inseg=brain-mets_pass_luad_absolute_prim-only.lgr.seg

refseg=brain-mets_gatk4cnv_no-pkb-pon_luad_normal-combined.lgr.cnv.seg

centromeres=~/data/ucsc/hg19/centromeres.seg

genomic=~/projects/genomic/build/genomic

# filter out segments with few probe/target counts
$genomic clean $inseg ${inseg%.*}.cntf.seg \
	--count 5 \
	--ref_state 0.0 --state_diff 0.15 --merge 1

# filter out segments overlapping highly with centromeres
$genomic filter ${inseg%.*}.cntf.seg $centromeres ${inseg%.*}.cntf.cenf.seg \
	--ref_state 0.0 --state_diff 0.15 --aberrant 1 --merge 0 \
	--threshold 0.5

# filter out segments overlapping highly with CNVs found in normals
$genomic filter ${inseg%.*}.cntf.cenf.seg $refseg ${inseg%.*}.cntf.cenf.cnvf.seg \
	--ref_state 0.0 --state_diff 0.15 --aberrant 0 --merge 1

