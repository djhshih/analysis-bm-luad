#!/bin/bash

indir=/seq/picard_aggregation/G94791/
outdir=qc

mkdir -p $outdir

samples=$(cut -f 1 bam_paths.txt | sed 1d)

for sample in ${samples[@]}; do
	cp $indir/$sample/current/*.*_metrics $outdir
done
