library(io);
library(dplyr);

pheno <- qread("sample-info_wes_stage2_pass_luad.tsv");
pheno.all <- qread("sample-info_wes_stage2.tsv");

dshih.patients <- as.character(unique(pheno$clinical_id));

# scarter included one met per patient
scarter.met.samples <- qread("sample-list_scarter.txt", type="vtr");
scarter.met.samples[scarter.met.samples == "PB0441-MT-PB0441-N"] <- "PB-PB0441-TM-NT-SM-3WKYB-SM-3WKYC";

filter(pheno.all, grepl("PB0441", fh_pair_id))

scarter.patients <- as.character(pheno.all$clinical_id[match(scarter.met.samples, pheno.all$fh_pair_id)]);

any(duplicated(scarter.patients))
any(duplicated(dshih.patients))

# removed patients
length(setdiff(scarter.patients, dshih.patients))
setdiff(scarter.patients, dshih.patients)

# common patients
length(intersect(scarter.patients, dshih.patients))
intersect(scarter.patients, dshih.patients)

# added patients
length(setdiff(dshih.patients, scarter.patients))
setdiff(dshih.patients, scarter.patients)

