library(io);
library(dplyr);

scarter <- readLines("sample-list_scarter.txt");

pheno <- qread("sample-info_wes_stage2.tsv");
dshih <- qread("sample-info_wes_stage2_pass_luad.tsv");

filter(
	pheno,
	grepl("PB0441", pheno$sample_id)
);

grep("0052", scarter)

scarter.mod <- scarter;
scarter.mod[grep("441", scarter)] <- as.character(pheno$fh_pair_id[pheno$sample_id == "PB0441-M"]);

scarter.changed <- select(pheno[match(scarter.mod, pheno$fh_pair_id), ],
	sample_id,
	fh_sample_id,
	fh_pair_id,
	primary_cancer_type,
	primary_histotype,
	sample_type
) %>% mutate(
	included = sample_id %in% dshih$sample_id,
	scarter_id = scarter
);

table(as.character(scarter.changed$primary_histotype))

qwrite(scarter.changed, "sample-changes_scarter.tsv");

# Sample changes:
# /home/unix/dshih/exigens/brain-mets/annot/sample-changes_scarter.tsv
# Regarding PB0441-MT-PB0441-NB, this pair was duplicated (like a handful of
# others) on Firehose. I standardized to the pair with a FH unique identifier
# (e.g. PB-PB0441-TM-NT-SM-3WKYB-SM-3WKYC).

