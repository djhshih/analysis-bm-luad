# Aggregate and organize all sample information tables

library(dplyr);
library(magrittr);
library(io);

source("df_helper.R");

options(stringsAsFactors=FALSE);

clin0 <- qread("Clinical_Data_Absolute_Outputs_ds_2016-04-28.tsv", comment="");
clin1 <- qread("bm_clinical_idj_2016-03-06.tsv", comment="");
clin2 <- qread("lung_bm_idj_2016-03-06.tsv", comment="");

# missing information from bm_clinical_idj:
# date of diagnosis
#   closest is lung biopsy date, but this field is mostly empty
# date of death or last follow-up
# date of brain met diagnosis: radiological or pathological?
# targeted treatment is based on primary bx?
# what are the dates of progression?
#   some start after brain met diagnosis, some before

normalize_smoking_hx <- function(d) {
	d$smoking_hx %<>%
		sub("(\\d+-)?<?>? ?(\\d+)\\+? (pack )?year.*", "\\2", .) %>%
		sub("never", "0", ., ignore.case=TRUE) %>%
		sub("unk", "", ., ignore.case=TRUE) %>%
		as.numeric
	d
}

clin0.n <- clin0 %>% normalize_df_field_names %>% 
	rename_df_field("primary_sid", "primary_pair_id") %>%
	rename_df_field("met_sid", "met_pair_id") %>%
	rename_df_field("smoking_hx", "smoker") %>%
	rename_df_field("smoking_pack_year_hx", "smoking_hx") %>%
	rename_df_field("date_of_dx_of_primary", "date_of_primary_dx") %>%
	mutate(
		death = as.numeric(!grepl("alive|visit", tolower(date_of_death))),
		date_of_primary_dx = normalize_date(date_of_primary_dx),
		date_of_death = normalize_date(date_of_death),
		os_years = as.numeric(date_of_death - date_of_primary_dx) / 365,
		date_of_bm_dx = date_of_primary_dx + bm_pfs_years * 356,
		primary_site = normalize_empty(primary_site),
		cancer_type = normalize_empty(cancer_type),
		histotype = normalize_empty(histotype),
		histosubtype = normalize_empty(histosubtype)
	) %>%
	mutate(
		gender = normalize_empty(gender)
	)


normalize_histotype <- function(histology) {
	x <- rep(NA, length(histology));
	x[grepl("NSCLC", histology)] <- "Non-small cell lung carcinoma";
	x[grepl("adenocarcinoma", histology)] <- "Lung adenocarcinoma";
	x[grepl("adenosquam", histology)] <- "Adenosquamous cell lung carcinoma";
	x[grepl("squamous", histology)] <- "Squamous cell lung carcinoma";
	x
}

clin1.n <- clin1 %>% normalize_df_field_names %>%
	rename_df_field("bm_dx", "date_of_bm_dx") %>%
	rename_df_field("tx_bm", "bm_tx") %>%
	normalize_smoking_hx %>%
	mutate(
		date_of_bm_dx = normalize_date(date_of_bm_dx),
		date_of_lung_bx = normalize_date(date_of_lung_bx),
		date_of_lung_rsxn = normalize_date(date_of_lung_rsxn),
		date_of_primary_dx = pmin(date_of_lung_bx, date_of_lung_rsxn, na.rm=TRUE),
		bm_pfs_years = as.numeric(date_of_bm_dx - date_of_primary_dx) / 365
	) %>%
	mutate(
		histotype = normalize_histotype(histology),
		race = ifelse(race == "unknown", NA, capitalize(race)),
		gender = normalize_empty(gender)
	) %>%
	select(-histology)


normalize_histotype2 <- function(histology) {
	x <- rep(NA, length(histology));
	x[grepl("7", histology)] <- "Small cell lung carcinoma";
	x[grepl("6", histology)] <- "Large cell lung carcinoma";
	x[grepl("5", histology)] <- "Non-small cell lung carcinoma";
	x[grepl("4", histology)] <- "Pleomorphic lung carcinoma";
	x[grepl("2", histology)] <- "Squamous cell lung carcinoma";
	x[grepl("1", histology)] <- "Lung adenocarcinoma";
	x[grepl("3", histology)] <- "Adenosquamous lung carcinoma";
	x
}

normalize_gender2 <- function(gender) {
	x <- rep(NA, length(gender));
	x[gender == 1] <- "Male";
	x[gender == 2] <- "Female";
}

clin2.n <- clin2 %>% normalize_df_field_names %>%
	rename_df_field("bm_dx", "date_of_bm_dx") %>%
	rename_df_field("tx_bm", "bm_tx") %>%
	normalize_smoking_hx %>%
	mutate(
		date_of_bm_dx = normalize_date(date_of_bm_dx),
		date_of_lung_bx = normalize_date(date_of_lung_bx),
		date_of_lung_rsxn = normalize_date(date_of_lung_rsxn),
		date_of_primary_dx = pmin(date_of_lung_bx, date_of_lung_rsxn, na.rm=TRUE),
		bm_pfs_years = as.numeric(date_of_bm_dx - date_of_primary_dx) / 365
	) %>%
	mutate(
		histotype = normalize_histotype2(histology),
		gender = as.character(factor(gender, levels=1:2, labels=c("M", "F"))),
		race = as.character(factor(race, levels=1:4,
			labels=c("White", "Black", "Asian", "Hispanic")))
	) %>%
	select(-histology)

# remove due to disagreement
clin2.n[clin2.n$clinical_id == "PB0147", "bm_pfs_years"] <- NA;


# check that there are no disagreements

clin.n <- clin0.n %>% 
	merge(clin1.n, by="clinical_id", all=TRUE) %>%
	merge(clin2.n, by="clinical_id", all=TRUE);

dim(clin.n)

intersect(names(clin0.n), names(clin1.n)) %>% intersect(names(clin2.n));

# check that tables do not disagree

select(clin.n, clinical_id, gender, gender.x, gender.y);
filter(clin.n, gender != gender.x);
filter(clin.n, gender != gender.y);

select(clin.n, clinical_id, histotype, histotype.x, histotype.y);
filter(clin.n, histotype != histotype.x);
filter(clin.n, histotype != histotype.y);

select(clin.n, clinical_id, smoking_hx, smoking_hx.x, smoking_hx.y);
filter(clin.n, smoking_hx != smoking_hx.x);
filter(clin.n, smoking_hx != smoking_hx.y);

select(clin.n, clinical_id, bm_pfs_years, bm_pfs_years.x, bm_pfs_years.y);
filter(clin.n, bm_pfs_years != bm_pfs_years.x);
filter(clin.n, bm_pfs_years != bm_pfs_years.y);

select(clin.n, clinical_id, date_of_primary_dx, date_of_primary_dx.x, date_of_primary_dx.y);
filter(clin.n, date_of_primary_dx != date_of_primary_dx.x);
filter(clin.n, date_of_primary_dx != date_of_primary_dx.y);

select(clin.n, clinical_id, date_of_bm_dx, date_of_bm_dx.x, date_of_bm_dx.y);
filter(clin.n, date_of_bm_dx != date_of_bm_dx.x);
filter(clin.n, date_of_bm_dx != date_of_bm_dx.y);


clin0.n$clinical_id %>%
	intersect(clin1.n$clinical_id) %>%
	intersect(clin2.n$clinical_id);

clin0.n$clinical_id %>%
	union(clin1.n$clinical_id) %>%
	union(clin2.n$clinical_id);

first_valid <- function(x) {
	x[which.max(!is.na(x))]
}

clin.n <- clin0.n %>% 
	full_join(clin1.n) %>%
	full_join(clin2.n) %>%
	group_by(clinical_id) %>%
	summarize_each(funs(first_valid));

dim(clin.n)



filter(clin.n, histotype == "Lung adenocarcinoma", is.na(bm_pfs_years));

filter(clin.n, histotype == "Lung adenocarcinoma", is.na(gender));

filter(clin.n, cancer_type == "Lung cancer", is.na(histotype));
# PB0225 does not even have gender

filter(clin.n, histotype == "Non-small cell lung carcinoma");
# PB0059 is listed in IDJ's table but has no information
# LS-014 is poorly differentiated
# PB0063 is poorly differentiated
# PB0390 is poorly differentiated
# PB0046 is ambiguous: Lung adenocarcinoma and Large cell lung carcinoma

nrow(filter(clin.n, rnaseq))

rnaseq <- qread("sample-info_rnaseq.tsv");

unique(rnaseq$patient_id) %in% clin.n$patient_id

filter(clin.n, patient_id %in% unique(rnaseq$patient_id))
filter(clin.n, patient_id %in% unique(rnaseq$patient_id), !rnaseq)


clin.n %<>% remove_df_field("other");

qwrite(clin.n, filename("patient-info", ext="tsv"), quote=TRUE);
qwrite(clin.n, filename("patient-info", ext="rds"));

