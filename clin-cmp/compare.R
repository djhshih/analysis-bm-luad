library(io);
library(dplyr);
library(ggplot2);
library(reshape2);
library(survival);
library(survminer);


pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3.tsv");
pkb.clin.all <- qread("~/exigens/brain-mets/annot/patient-info_stage2.tsv");

tcga.pheno.all <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");
tcga.clin.all <- qread("~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad_stage2.rds");

pdf.fname <- filename("pkb-tcga-compare", ext="pdf");
pkb.weights.fname <- filename("pkb-luad", tag="weights", ext="tsv");
tcga.weights.fname <- filename("tcga-luad", tag="weights", ext="tsv");

source("~/exigens/brain-mets/common/params.R");

####

pkb.pheno <- filter(pkb.pheno.all, primary_histotype == "Lung adenocarcinoma");
pkb.patients <- as.character(unique(pkb.pheno$clinical_id));

pkb.clin <- filter(pkb.clin.all, clinical_id %in% pkb.patients);


tcga.pheno <- filter(tcga.pheno.all, sample_type == "Tumor");
tcga.patients <- as.character(unique(tcga.pheno$clinical_id));

tcga.clin <- filter(tcga.clin.all, clinical_id %in% tcga.patients) %>%
	mutate(
		age_at_primary_dx = age_at_initial_pathologic_diagnosis,
		smoking_hx = number_pack_years_smoked
	)

setdiff(tcga.patients, tcga.clin.all$clinical_id)

pkb.n <- nrow(pkb.clin);
tcga.n <- nrow(tcga.clin);

pkb.clin$stage[pkb.clin$stage == "IIA/IIB"] <- "IIB";

stages.clean <- c("I", "II", "III", "IV");

pkb.clin <- mutate(pkb.clin,
	stage_clean = factor(gsub("A|B", "", stage), levels=stages.clean)
);

tcga.clin$stage_clean <- factor(toupper(gsub("a|b", "", gsub("stage ", "", tcga.clin$pathologic_stage))),
	levels=stages.clean);

tcga.clin <- mutate(tcga.clin,
	os_years = ifelse(is.na(days_to_last_followup), days_to_death, days_to_last_followup) / 365,
	death = as.integer(vital_status == "dead")
);

####

# days_to_last_known_alive appears to be useless
# (days_to_last_followup or days_to_death) are always longer
select(tcga.clin, days_to_last_followup, days_to_last_known_alive, days_to_death, vital_status) %>% 
	filter(!is.na(days_to_last_known_alive))

# days_to_initial_pathologic_diagnosis is always 0
# therefore, the origin of the time frame is at initial diagnosis

####

all.clin <- rbind(
	select(tcga.clin, clinical_id, os_years, death, stage=stage_clean) %>% mutate(cohort=cohorts[1]),
	select(pkb.clin, clinical_id, os_years, death, stage=stage_clean) %>% mutate(cohort=cohorts[2])
);

pkb.fit <- survfit(Surv(os_years, death) ~ stage_clean, pkb.clin);
ggsurvplot(pkb.fit, risk.table=TRUE);

tcga.fit <- survfit(Surv(os_years, death) ~ stage_clean, tcga.clin);
ggsurvplot(tcga.fit, risk.table=TRUE);

all.fit <- survfit(Surv(os_years, death) ~ stage + cohort, all.clin);
ggsurvplot(all.fit, color="stage", linetype="cohort", risk.table=TRUE);

