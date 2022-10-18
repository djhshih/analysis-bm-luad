#!/bin/bash

indir=$1

printf "sample_id\tcount\n"

for f in $indir/*.maf; do
	count=$(grep -E "(KEEP)|(PASS)" $f | wc -l)
	fbase=${f##*/}
	fstem=${fbase%%.*}
	printf "$fstem\t$count\n"
done

