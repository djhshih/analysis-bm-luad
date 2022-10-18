#!/bin/bash

cut -f1 /home/unix/scarter/Projects/Brain_met_2/mutsig/mutsig_results/Bmet_unmatched_LUAD_and_carcinoma/patient_counts_and_rates.txt |
	sed 1d > sample-list_scarter.txt
