# Aggregate and organize all sample information tables

library(dplyr);
library(magrittr);
library(io);

source("df_helper.R");

clin0 <- qread("Clinical_Data_Absolute_Outputs_ds_2016-12-05.tsv", comment="", quote="\"");
clin1 <- qread("bm_clinical_idj_ds_ak_unc_2017-04-10.tsv", comment="", quote="\"");
clin2 <- qread("lung_bm_idj_ds_ak_unc_2017-04-03.tsv", comment="", quote="\"");
clin3 <- qread("MasterList-PB-deidentified_ds_2016-10-03.tsv", quote="\"");

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
		date_of_primary_dx = normalize_date(date_of_primary_dx),
		date_of_bm_dx = date_of_primary_dx + bm_pfs_years * 356,
		death = as.numeric(!grepl("alive|visit", tolower(date_of_death))),
		date_of_death = normalize_date(date_of_death),
		os_years = as.numeric(date_of_death - date_of_primary_dx) / 365,
		primary_site = normalize_empty(primary_site),
		bm_location = normalize_empty(bm_location),
		histotype = normalize_empty(histotype),
		histosubtype = normalize_empty(histosubtype)
	) %>%
	mutate(
		gender = normalize_empty(gender),
		age_at_primary_dx = as.integer(age_of_dx_of_primary),
		tumor_stage = normalize_empty(as.character(stage_of_dx_of_primary)),
		node_status = factor(nodal_involvement_at_dx_of_primary, levels=0:1, labels=c("N0", "N+")),
		met_status = ifelse(is.na(bm_pfs_years) | bm_pfs_years >= 0, NA, "M1b")
	) %>%
	select(-age_of_dx_of_primary)

normalize_histotype <- function(histology) {
	x <- rep(NA, length(histology));
	x[grepl("NSCLC", histology)] <- "Non-small cell lung carcinoma";
	x[grepl("adenocarcinoma", histology)] <- "Lung adenocarcinoma";
	x[grepl("adenosquam", histology)] <- "Adenosquamous cell lung carcinoma";
	x[grepl("squamous", histology)] <- "Squamous cell lung carcinoma";
	x
}

age_at_primary_dx <- function(date_of_birth, date_of_primary_dx) {
	floor(as.double(date_of_primary_dx - date_of_birth, units="days") / 365)
}

clin1.n <- clin1 %>% normalize_df_field_names %>%
	convert_df_na_fields_to_type("character") %>%
	rename_df_field("bm_dx", "date_of_bm_dx") %>%
	rename_df_field("tx_bm", "bm_tx") %>%
	normalize_smoking_hx %>%
	mutate(
		brain_rt_prior_to_met_sx = as.integer(brain_rt_prior_to_met_sx),
		date_of_birth = normalize_date(date_of_birth),
		date_of_bm_dx = normalize_date(date_of_bm_dx),
		date_of_lung_bx = normalize_date(date_of_lung_bx),
		date_of_lung_sx = normalize_date(date_of_lung_sx),
		date_of_primary_dx = pmin(normalize_date(date_of_primary_dx), date_of_lung_bx, date_of_lung_sx, na.rm=TRUE),
		death = as.numeric(!grepl("alive|visit", tolower(date_of_death))),
		date_of_death = normalize_date(date_of_death),
		os_years = as.numeric(date_of_death - date_of_primary_dx) / 365,
		bm_pfs_years = as.numeric(date_of_bm_dx - date_of_primary_dx) / 365
	) %>%
	mutate(
		cancer_type = "Lung cancer",
		bm_location = normalize_empty(bm_location),
		histotype = normalize_histotype(histology),
		histosubtype = normalize_empty(histosubtype),
		race = ifelse(race == "unknown", NA, trimws(capitalize(race))),
		gender = normalize_empty(gender),
		age_at_primary_dx = age_at_primary_dx(date_of_birth, date_of_primary_dx),
		bm_processed = 1
	) %>%
	select(-histology) %>%
	select(-contains("notes"))


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
	convert_df_na_fields_to_type("character") %>%
	rename_df_field("bm_dx", "date_of_bm_dx") %>%
	rename_df_field("tx_bm", "bm_tx") %>%
	normalize_smoking_hx %>%
	mutate(
		brain_rt_prior_to_met_sx = as.integer(brain_rt_prior_to_met_sx),
		date_of_birth = normalize_date(date_of_birth),
		date_of_bm_dx = normalize_date(date_of_bm_dx),
		date_of_lung_bx = normalize_date(date_of_lung_bx),
		date_of_lung_sx = normalize_date(date_of_lung_sx),
		date_of_primary_dx = pmin(normalize_date(date_of_primary_dx), date_of_lung_bx, date_of_lung_sx, na.rm=TRUE),
		death = as.numeric(!grepl("alive|visit", tolower(date_of_death))),
		date_of_death = normalize_date(date_of_death),
		os_years = as.numeric(date_of_death - date_of_primary_dx) / 365,
		bm_pfs_years = as.numeric(date_of_bm_dx - date_of_primary_dx) / 365
	) %>%
	mutate(
		cancer_type = "Lung cancer",
		bm_location = normalize_empty(bm_location),
		histotype = normalize_histotype2(histology),
		histosubtype = normalize_empty(histosubtype),
		gender = as.character(factor(gender, levels=1:2, labels=c("M", "F"))),
		age_at_primary_dx = age_at_primary_dx(date_of_birth, date_of_primary_dx),
		race = as.character(factor(race, levels=1:4,
			labels=c("White", "Black", "Asian", "Hispanic"))),
		bm_processed = 1
	) %>%
	select(-histology) %>%
	select(-contains("notes"))

# remove due to disagreement with clin2.n
clin0.n[clin0.n$clinical_id == "PB0147", "bm_pfs_years"] <- NA; 

pkbmid_to_pid <- function(x) {
	ifelse(is.na(x),
		x,
		paste0("PB", sub("-MT.*", "", x))
	)
}

grepli <- function(pattern, x, ...) {
	grepl(pattern, x, ignore.case=TRUE)
}

normalize_cancer_type <- function(x) {
	y <- rep(NA, length(x));
	y[grepli("lung", x)] <- "Lung cancer";
	y[grepli("small cell", x)] <- "Lung cancer";
	y[grepli("NSCLC", x)] <- "Lung cancer";
	y[grepli("ovarian", x)] <- "Ovarian cancer";
	y[grepli("pancreas", x)] <- "Pancreatic cancer";
	y[grepli("thyroid", x)] <- "Thyroid cancer";
	y[grepli("adrenocortical", x)] <- "Adrenal cancer";
	y[grepli("colon", x)] <- "Colorectal cancer";
	y[grepli("rectal", x)] <- "Colorectal cancer";
	y[grepli("endometrial", x)] <- "Endometrial cancer";
	y[grepli("breast", x)] <- "Breast cancer";
	y[grepli("esophageal", x)] <- "Esophageal cancer";
	y[grepli("laryngeal", x)] <- "Laryngeal cancer";
	y[grepli("melanoma", x)] <- "Skin cancer";
	y[grepli("rcc", x)] <- "Kidney cancer";
	y[grepli("renal", x)] <- "Kidney cancer";
	y[grepli("urothel", x)] <- "Urothelial cancer";
	y[grepli("testic", x)] <- "Testicular cancer";
	y[grepli("vaginal", x)] <- "Vaginal cancer";
	y[grepli("sarcoma", x)] <- "Sarcoma";
	y[grepli("submandibular", x)] <- "Salivary gland cancer";
	y[grepli("prostate", x)] <- "Prostate cancer";
	y[grepli("neck", x)] <- "Head and neck cancer";
	y[grepli(" or ", x)] <- NA;
	y
}

normalize_histotype3 <- function(x) {
	y <- rep(NA, length(x));
	y[grepli("small cell", x)] <- "Small cell lung carcinoma";
	y[grepli("non small cell", x)] <- "Non-small cell lung carcinoma";
	y[grepli("NSCLC", x)] <- "Non-small cell lung carcinoma";
	y[grepli("lung squamous", x)] <- "Squamous cell lung carcinoma";
	y[grepli("lung adeno", x)] <- "Lung adenocarcinoma";
	y[grepli("lung adeno/squam", x)] <- "Adenosquamous lung carcinoma";
	y[grepli("lung neuroendorcine", x)] <- "Neuroendocrine lung carcinoma";
	y
}

clin3.n <- clin3 %>% normalize_df_field_names %>%
	convert_df_na_fields_to_type("character") %>%
	mutate(
		clinical_id = pkbmid_to_pid(normalize_empty(brain_specimen_reference_id)),
		dx = normalize_empty(dx),
		cancer_type = normalize_cancer_type(dx),
		histotype = normalize_histotype3(dx)
	) %>%
	select(
		clinical_id, cancer_type, histotype,
		primary_processed = primary_sent_for_sequencing,
		bm_processed = bm_sent_for_sequencing
	) %>%
	filter(
		! (is.na(clinical_id) | is.na(cancer_type))
	)


# check that there are no disagreements

clin.n <- clin0.n %>% 
	merge(clin1.n, by="clinical_id", all=TRUE) %>%
	merge(clin2.n, by="clinical_id", all=TRUE) %>%
	merge(rename(clin3.n, histotype.z=histotype, cancer_type.z=cancer_type), 
		by="clinical_id", all=TRUE)

dim(clin.n)

intersect(names(clin0.n), names(clin1.n)) %>% intersect(names(clin2.n));

# check that tables do not disagree

select(clin.n, clinical_id, gender, gender.x, gender.y);
# should be empty
filter(clin.n, gender != gender.x);
filter(clin.n, gender != gender.y);

select(clin.n, clinical_id, cancer_type, cancer_type.x, cancer_type.y, cancer_type.z);
# should be empty
filter(clin.n, cancer_type != cancer_type.x);
filter(clin.n, cancer_type != cancer_type.y);
filter(clin.n, cancer_type != cancer_type.z);
filter(clin.n, cancer_type.x != cancer_type.z);
# TODO resolve discrepancies here...
# for now, x trumps z
filter(clin.n, cancer_type.y != cancer_type.z);

select(clin.n, clinical_id, histotype, histotype.x, histotype.y, histotype.z);
# should be empty
filter(clin.n, histotype != histotype.x);
filter(clin.n, histotype != histotype.y);
filter(clin.n, histotype != histotype.z);
# LUAD is the more specific diagnosis for PB0341;
# TTF-1 staining confirms diagnosis

select(clin.n, clinical_id, smoking_hx, smoking_hx.x, smoking_hx.y);
# should be empty
filter(clin.n, smoking_hx != smoking_hx.x);
filter(clin.n, smoking_hx != smoking_hx.y);
filter(clin.n, smoking_hx.x != smoking_hx.y);
# TODO resolve discrepancy here...
# For now, smoking hx 50 vs. 20 is not a big deal

select(clin.n, clinical_id, bm_pfs_years, bm_pfs_years.x, bm_pfs_years.y);
# should be empty
filter(clin.n, bm_pfs_years != bm_pfs_years.x) %>% select(clinical_id, bm_pfs_years, bm_pfs_years.x)
# discrepancy is due to rounding
filter(clin.n, bm_pfs_years != bm_pfs_years.y);
filter(clin.n, bm_pfs_years.x != bm_pfs_years.y);

select(clin.n, clinical_id, date_of_primary_dx, date_of_primary_dx.x, date_of_primary_dx.y);
# should be empty
filter(clin.n, date_of_primary_dx != date_of_primary_dx.x) %>% select(clinical_id, date_of_primary_dx, date_of_primary_dx.x)
filter(clin.n, date_of_primary_dx != date_of_primary_dx.y);
filter(clin.n, date_of_primary_dx.x != date_of_primary_dx.y);

select(clin.n, clinical_id, date_of_bm_dx, date_of_bm_dx.x, date_of_bm_dx.y);
# should be empty
filter(clin.n, date_of_bm_dx != date_of_bm_dx.x) %>% select(clinical_id, date_of_bm_dx, date_of_bm_dx.x)
filter(clin.n, date_of_bm_dx != date_of_bm_dx.y);
filter(clin.n, date_of_bm_dx.x != date_of_bm_dx.y);


clin0.n$clinical_id %>%
	intersect(clin1.n$clinical_id) %>%
	intersect(clin2.n$clinical_id);

clin0.n$clinical_id %>%
	union(clin1.n$clinical_id) %>%
	union(clin2.n$clinical_id);


# summarize_each(funs(.)) cannot take NA as output for some reason...
first_valid <- function(x) {
	# find first non-null value
	x[which.max(!is.na(x))]
}

clin.n <- clin0.n %>% 
	full_join(clin1.n) %>%
	full_join(clin2.n) %>%
	full_join(clin3.n) %>%
	group_by(clinical_id) %>%
	summarize_each(funs(first_valid));

clin.n %<>%
	mutate(
		histotype = factor(normalize_empty(histotype)),
		cancer_type = factor(normalize_empty(cancer_type)),
		gender = factor(gender, levels = c("F", "M"), labels = c("Female", "Male"))
	) %>%
	rename(
		primary_histotype = histotype,
		primary_histosubtype = histosubtype,
		primary_cancer_type = cancer_type
	);

dim(clin.n)

filter(clin.n, primary_histotype == "Lung adenocarcinoma", is.na(bm_pfs_years));

filter(clin.n, primary_histotype == "Lung adenocarcinoma")$bm_processed %>% table();

filter(clin.n, primary_histotype == "Lung adenocarcinoma", is.na(gender));
# *PB0225 does not even have gender

filter(clin.n, primary_cancer_type == "Lung cancer", is.na(primary_histotype));
# PB0052 is likely not a lung cancer??? It is an adrenocortical carcinoma.

filter(clin.n, primary_histotype == "Non-small cell lung carcinoma") %>% as.data.frame()
# LS-014 is poorly differentiated
# PB0036 is poorly differentiated
# PB0046 is ambiguous: Lung adenocarcinoma and Large cell lung carcinoma --
#        reclassified as Lung adenocarcinoma
# PB0063 is poorly differentiated
# PB0254 is a external NSCLC with no additional information
# PB0390 is poorly differentiated
# PB0441 is poorly differentiated
# PB0425 has no information or sample


# Calculate derivative values

smoker.orig <- clin.n$smoker;

clin.n <- mutate(clin.n,
	# brain met is synchronous if brain met occurs within 3 months
	synchronous_bm = factor(as.integer((date_of_bm_dx - date_of_primary_dx) < 90), levels=0:1, labels=c("Metachronous", "Synchronous")),
	# never smoker: pack-year == 0
	smoker = smoking_hx > 0,
	# post brain met progression survival
	bm_pps_years = as.numeric(date_of_death - date_of_bm_dx) / 365,
	# all patients in cohort have brain met
	bm = 1,
	# set NA to its own level
	tumor_stage = ifelse(is.na(tumor_stage), "TX", as.character(tumor_stage)),
	node_status = ifelse(is.na(node_status), "NX", as.character(node_status)),
	met_status = ifelse(is.na(met_status), "MX", as.character(met_status)),
	# set disease stage to IV if presentation with brain met
	stage = ifelse(is.na(bm_pfs_years) | bm_pfs_years >= 0, as.character(stage), "IV")
) %>% mutate(
	# set disease stage to IA if T1/T1a/Tb and N0 (assuming M0)
	stage = ifelse(tumor_stage %in% c("T1", "T1a", "T1b") & node_status == "N0", "IA", stage)
);
clin.n$smoker[smoker.orig == 1] <- TRUE;
clin.n$smoker[smoker.orig == 0] <- FALSE;
clin.n$smoker <- factor(clin.n$smoker, levels=c(FALSE, TRUE), labels=c("Never", "Ever"));


nrow(filter(clin.n, rnaseq))

# NB rnaseq uses slightly different names... sigh
#rnaseq.df <- qread("sample-info_rnaseq.tsv");

#unique(rnaseq.df$patient_id) %in% clin.n$patient_id

#filter(clin.n, patient_id %in% unique(rnaseq.df$patient_id))
#filter(clin.n, patient_id %in% unique(rnaseq.df$patient_id), !rnaseq)

#stopifnot(!any(is.na(clin.n$primary_cancer_type)));

filter(clin.n, primary_cancer_type == "Lung cancer") %>%
	group_by(primary_histotype) %>% summarize(count=n())

filter(clin.n, primary_histotype == "Adenosquamous lung carcinoma") %>% select(clinical_id)

filter(clin.n, primary_histotype == "Non-small cell lung carcinoma") %>% select(clinical_id)


clin.n.lung <- filter(clin.n, primary_cancer_type == "Lung cancer");

table(clin.n.lung$primary_histotype, useNA="always");

clin.n.luad <- filter(clin.n.lung, primary_histotype == "Lung adenocarcinoma");

qwrite(clin.n, filename("patient-info", ext="tsv"), quote=TRUE);
qwrite(clin.n, filename("patient-info", ext="rds"));

qwrite(clin.n.lung, filename("patient-info", tag="lung", ext="tsv"), quote=TRUE);
qwrite(clin.n.lung, filename("patient-info", tag="lung", ext="rds"));

qwrite(clin.n.luad, filename("patient-info", tag="luad", ext="tsv"), quote=TRUE);
qwrite(clin.n.luad, filename("patient-info", tag="luad", ext="rds"));

