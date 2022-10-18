library(io);
library(dplyr);
library(survival);
library(kmplot);

samples <- qread("sample-info_wes_stage2_pass_luad.tsv");

length(unique(samples$clinical_id))

length(unique(filter(samples, sample_type == "Primary")$clinical_id))

length(unique(filter(samples, sample_type == "Brain metastasis")$clinical_id))

length(unique(filter(samples, sample_type %in% c("Extracranial metastasis"))$clinical_id))
z <- unique(filter(samples, sample_type %in% c("Extracranial metastasis"))$clinical_id);
filter(samples, clinical_id %in% z);

length(unique(filter(samples, grepl("metastasis", sample_type))$clinical_id))

samples.luad <- filter(samples, primary_histotype == "Lung adenocarcinoma");

patients <- qread("patient-info.rds");

patients.luad <- filter(patients,
	primary_cancer_type == "Lung cancer",
	primary_histotype == "Lung adenocarcinoma",
	clinical_id %in% samples$clinical_id
) %>% as.data.frame();

days.to.bm <- as.numeric(with(patients.luad, date_of_bm_dx - date_of_primary_dx));
hist(days.to.bm, breaks=50)
abline(v=60)
abline(v=90)
table(days.to.bm < 60)
table(days.to.bm < 90)

survfit(with(patients.luad, Surv(bm_pps_years, death) ~ 1))

with(patients.luad, table(gender, useNA="always"))
with(patients.luad, table(synchronous_bm, useNA="always"))
with(patients.luad, table(brain_rt_prior_to_met_sx, useNA="always"))
with(patients.luad, table(smoker, useNA="always"))
with(patients.luad, table(tumor_stage, useNA="always"))
with(patients.luad, table(node_status, useNA="always"))
with(patients.luad, table(met_status, useNA="always"))
with(patients.luad, table(stage, useNA="always"))

filter(patients.luad, is.na(stage)) %>% select(clinical_id, tumor_stage, node_status, met_status)

with(patients.luad, summary(bm_pfs_years))
with(patients.luad, summary(os_years))

filter(patients.luad, is.na(bm_pfs_years)) %>% select(clinical_id)
filter(patients.luad, is.na(bm_pfs_years)) %>% select(clinical_id)

# probably not valid: since only brain met patients are observed
patients.luad$bm <- 1;
s <- with(patients.luad, Surv(bm_pfs_years, bm) ~ smoker);
#kmplot(s)

s <- with(patients.luad, Surv(os_years, death) ~ smoker);
kmplot(s, xtick.interval=1)

s <- with(patients.luad, Surv(os_years, death) ~ synchronous_bm);
kmplot(s, xtick.interval=1)
survfit(s)

#patients.luad$bm_pps_years <- pmax(patients.luad$bm_pps_years, 0.0001);
s <- with(patients.luad, Surv(bm_pps_years, death) ~ synchronous_bm);
kmplot(s, xtick.interval=1, ylab="Post brain met progression survival", xlab="Time (years)")

s <- with(patients.luad, Surv(os_years, bm) ~ gender);
kmplot(s, xtick.interval=1);

s <- with(patients.luad, Surv(bm_pfs_years, bm) ~ gender);
kmplot(s, xtick.interval=1);

s <- with(patients.luad, Surv(bm_pps_years, death) ~ gender);
kmplot(s, xtick.interval=1);

#s <- with(patients.luad, Surv(os_years, bm) ~ race);
#kmplot(s, xtick.interval=1);


