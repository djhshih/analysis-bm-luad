library(io);
library(dplyr);
library(ggplot2);
library(binom);
library(reshape2);

pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3.tsv");
pkb.clin.all <- qread("~/exigens/brain-mets/annot/patient-info_stage2.tsv");

tcga.pheno.all <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");
tcga.clin.all <- qread("~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad_stage2.rds");

pkb.weights <- qread("../pkb-luad_weights.tsv");
tcga.weights <- qread("../tcga-luad_weights.tsv");

pdf.fname <- filename("pkb-tcga-samples-compare", ext="pdf");


col.control <- "royalblue3";
col.case <- "darkorange";
cols.ch <- c(control=col.control, case=col.case);
cohorts <- c("TCGA", "Present");

cols.st <- c("Primary"="#56B4E9", "Brain metastasis"="#D55E00");

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

stopifnot(length(setdiff(tcga.patients, tcga.clin.all$clinical_id)) == 0)

pkb.n <- nrow(pkb.clin);
tcga.n <- nrow(tcga.clin);

pkb.clin$stage[pkb.clin$stage == "IIA/IIB"] <- "IIB";

pkb.clin <- mutate(pkb.clin,
	stage_clean = factor(gsub("A|B", "", stage), levels=c("I", "II", "III", "IV"))
);

tcga.clin$stage_clean <- toupper(gsub("a|b", "", gsub("stage ", "", tcga.clin$pathologic_stage)));

tcga.pheno <- mutate(tcga.pheno,
	covered_mb = somatic_mutation_covered_bases_capture / 1e6
);

pkb.pheno <- mutate(pkb.pheno,
	covered_mb = somatic_mutation_covered_bases_capture / 1e6
);

# merge in weights
pkb.clin <- left_join(pkb.clin, pkb.weights, by="clinical_id");
tcga.clin <- left_join(tcga.clin, tcga.weights, by="clinical_id");
pkb.pheno <- left_join(pkb.pheno, pkb.weights, by="clinical_id");
tcga.pheno <- left_join(tcga.pheno, tcga.weights, by="clinical_id");

tcga.pheno$sample_type <- as.character(tcga.pheno$sample_type);
tcga.pheno$sample_type[tcga.pheno$sample_type == "Tumor"] <- "Primary";
tcga.pheno$ffpe_q <- 100;
tcga.pheno$specimen_material <- "FF";
tcga.pheno$center <- "Broad";
