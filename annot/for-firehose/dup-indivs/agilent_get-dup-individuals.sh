#!/bin/bash

workspace=An_ALL_Agilent_BrainMets

for f in $(cat ../agilent_duplicated-individuals.txt); do

	fiss annot_get $workspace samples \
		indiv=$f \
		external_id_capture \
		> samples_${f}.tsv

done
