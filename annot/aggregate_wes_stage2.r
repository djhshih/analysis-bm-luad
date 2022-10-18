library(io);
library(dplyr);

source("df_helper.R");

options(stringsAsFactors=FALSE);

pheno <- qread("sample-info_wes.tsv");
sex <- qread("../sex/brain-mets_normals_sex.tsv");

picard.hs <- qread("../picard/selection-metrics.tsv") %>%
	normalize_df_field_names();

picard.dm <- qread("../picard/duplicate-metrics.tsv") %>%
	normalize_df_field_names();

#absolute <- qread("../exome-seq/SLC.ABSOLUTE.table.01.20.2016.tsv") %>%
#	normalize_df_field_names();
absolute <- qread("/xchip/scarter/dshih/exigens/brain-mets/absolute/ABSOLUTE_results/brain-mets_pass_luad_absolute-1-4/reviewed/absolute_brain-mets_pass_luad.tsv") %>%
	normalize_df_field_names();

abs.clin <- qread("Clinical_Data_Absolute_Outputs_ds_2016-12-05.tsv", comment="") %>%
	normalize_df_field_names();

samples.scarter <- qread("sample-list_scarter.txt", type="vtr");


# remove absolute statistics for low purity samples
idx <- absolute$call_status == "low purity" | absolute$call_status == "FAILED";
absolute[idx, -(1:3)] <- NA;

sample.stats.agilent <- qread("../fiss/brain-mets_agilent-clean_sample-stats.tsv") %>%
	normalize_df_field_names();
sample.stats.ice <- qread("../fiss/brain-mets_ice-20160822_sample-stats.tsv") %>%
	normalize_df_field_names();

pair.stats.agilent <- qread("../fiss/brain-mets_agilent-clean_pair-stats.tsv") %>%
	normalize_df_field_names();
pair.stats.ice <- qread("../fiss/brain-mets_ice-20160822_pair-stats.tsv") %>%
	normalize_df_field_names();

out.fname <- filename("brain-mets", tag="bams", ext="vtr");


sex <- mutate(sex,
	clinical_id = pheno$clinical_id[match(sample_id, pheno$sample_id)]
);

pheno <- mutate(pheno,
	sex = sex$sex[match(clinical_id, sex$clinical_id)]
);

stopifnot(with(pheno, !any(is.na(sex) & sample_type == "Normal")));


pheno <- left_join(pheno, 
	select(absolute, -array),
	by = c("sample_id" = "sample"));


pair.stats <- rbind(pair.stats.agilent, pair.stats.ice);
pheno <- left_join(pheno, pair.stats, by=c("fh_pair_id" = "pair_id"));

sample.stats <- rbind(sample.stats.agilent, sample.stats.ice) %>%
	select(sample_id, contamination_percentage_consensus_capture);

pheno <- left_join(pheno, sample.stats, by=c("fh_sample_id" = "sample_id"));

pheno <- left_join(pheno, select(picard.hs, sample_id, hs_target_coverage), by=c("sample_id" = "sample_id"));
pheno <- left_join(pheno, picard.dm, by=c("sample_id" = "sample_id"));


z <- pheno %>% group_by(clinical_id) %>% filter(sample_type %in% c("Primary", "Brain metastasis")) %>% summarize(n=n()) %>% 
	left_join(select(abs.clin, clinical_id, primary_met_relationship), by="clinical_id") %>% data.frame()

with(z, table(primary_met_relationship, useNA="always"))

with(z, table(n, primary_met_relationship, useNA="always"))

filter(z, n == 1, is.na(primary_met_relationship))
filter(pheno, clinical_id == "PB0064")
filter(pheno, clinical_id == "PB0325")
# mark these as Unmatched

filter(z, n == 2, is.na(primary_met_relationship))


filter(z, n == 6, primary_met_relationship == "Unrelated")$clinical_id
filter(pheno, clinical_id == "LS-025")
# NB  LS-025's primaries are actually local recurrences.
#     The original primary that gave rise to the met is not available...

unrelated <- filter(z, primary_met_relationship == "Unrelated")$clinical_id;
filter(pheno, clinical_id %in% unrelated) %>%
	select(sample_id, accession, primary_cancer_type, primary_histotype)

unmatched <- filter(z, primary_met_relationship == "Unmatched", n > 1)$clinical_id;
filter(pheno, clinical_id %in% unmatched) %>%
	select(sample_id, accession, primary_cancer_type, primary_histotype, call_status)


pheno <- mutate(
	pheno,
	pass_qc =
		contamination_percentage_consensus_capture < 10 &
		somatic_mutation_covered_bases_capture > 10e6  # removes normal samples
);

pheno.pass.tumours <- data.frame(filter(pheno, pass_qc), stringsAsFactors=TRUE);

# re-add normal samples from patients with at least one passing tumour
patients.pass <- unique(pheno.pass.tumours$clinical_id);
pheno.pass.normals <- filter(pheno, sample_type == "Normal", clinical_id %in% patients.pass);
pheno.pass.normals$pass_qc <- TRUE;
pheno.pass <- rbind(pheno.pass.tumours, pheno.pass.normals);
pheno.pass <- pheno.pass[with(pheno.pass, order(clinical_id, sample_type)), ];


length(unique(pheno$fh_pair_id))
length(unique(pheno.pass$fh_pair_id))

pheno.luad.pass <- filter(pheno.pass, primary_histotype == "Lung adenocarcinoma")
dim(pheno.luad.pass)
length(unique(pheno.luad.pass$clinical_id))
table(pheno.luad.pass$call_status, useNA="always")
summary(pheno.luad.pass)
length(unique(pheno.luad.pass$fh_pair_id))

pheno.lung.pass <- filter(pheno.pass, primary_cancer_type == "Lung cancer");
dim(pheno.lung.pass)
length(unique(pheno.lung.pass$clinical_id))
table(pheno.lung.pass$call_status, useNA="always")
summary(pheno.lung.pass)
length(unique(pheno.lung.pass$fh_pair_id))

pheno.nsclc.pass <- filter(pheno.pass,
	primary_cancer_type == "Lung cancer",
	primary_histotype != "Small cell lung carcinoma"
);
dim(pheno.nsclc.pass)
length(unique(pheno.nsclc.pass$clinical_id))
table(pheno.nsclc.pass$call_status, useNA="always")
summary(pheno.nsclc.pass)
length(unique(pheno.nsclc.pass$fh_pair_id))

patients.scarter <- filter(pheno, fh_pair_id %in% samples.scarter)$clinical_id;
samples.scarter[!(samples.scarter %in% pheno$fh_pair_id)]
filter(pheno, clinical_id == "PB0441")
patients.scarter  <- c(patients.scarter, "PB0441");

pheno.scarter <- filter(pheno, clinical_id %in% patients.scarter);

table(filter(pheno.scarter, sample_type == "Normal")$primary_histotype)

table(filter(pheno.luad.pass, sample_type == "Normal")$primary_histotype)
table(filter(pheno.nsclc.pass, sample_type == "Normal")$primary_histotype)


table(pheno.scarter$center)
table(pheno.luad.pass$center)

qwrite(pheno, filename("sample-info", tag=c("wes", "stage2"), ext="tsv"));
qwrite(pheno.pass, filename("sample-info", tag=c("wes", "stage2", "pass"), ext="tsv"));

qwrite(pheno.luad.pass, filename("sample-info", tag=c("wes", "stage2", "pass", "luad"), ext="tsv"));
qwrite(pheno.lung.pass, filename("sample-info", tag=c("wes", "stage2", "pass", "lung"), ext="tsv"));
qwrite(pheno.nsclc.pass, filename("sample-info", tag=c("wes", "stage2", "pass", "nsclc"), ext="tsv"));

qwrite(pheno.scarter, filename("sample-info", tag=c("wes", "stage2", "pass", "scarter"), ext="tsv"));

