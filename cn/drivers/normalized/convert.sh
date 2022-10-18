#!/bin/bash

samples=cna-samples-any.vtr

indir=~/share/exigens/brain-mets/gatk4cnv/gathered/tumor_pcov
outdir=rds

mkdir -p $outdir

for s in $(cat $samples); do
	echo $s
	./convert.R $indir/${s}.tn.tsv --outdir $outdir
done

