# Aggregate and organize all WES sample tables

library(magrittr);
library(io);
library(dplyr);

source("df_helper.R");

options(stringsAsFactors=FALSE);

# PB0403-N is superceded by PB0403-N-RE
#blacklist <- c("PB0403-N", "BS-532-P");

ice.in.fname <- as.filename("../fiss/brain-mets_ice-20160822");

ice.bams <- qread(insert(ice.in.fname, tag="bams", ext="tsv"));
ice.pairs <- qread(insert(ice.in.fname, tag="pairs", ext="tsv"));

ice <- left_join(ice.bams, ice.pairs,
		by = c("sample_id" = "case_sample")
	) %>%
	rename(fh_sample_id = sample_id, fh_pair_id = pair_id);
stopifnot(!any(duplicated(ice$fh_sample_id)))


agilent.in.fname <- as.filename("../fiss/brain-mets_agilent-clean");

agilent.bams <- qread(insert(agilent.in.fname, tag="bams", ext="tsv"));
agilent.pairs <- qread(insert(agilent.in.fname, tag="pairs", ext="tsv"));

agilent <- left_join(agilent.bams, agilent.pairs,
		by = c("sample_id" = "case_sample")
	) %>%
	rename(fh_sample_id = sample_id, fh_pair_id = pair_id);
stopifnot(!any(duplicated(agilent$fh_sample_id)))

patients <- qread("patient-info.rds");

samples <- qread("sample_types_2016-07-26.tsv");

#absolute <- qread("SLC.ABSOLUTE.table.07.24.2016.txt", type="tsv") %>% normalize_df_field_names;

samples.pkb <- qread("MasterList-PB-deidentified_ds_2016-10-03.tsv", quote="\"") %>% normalize_df_field_names;
samples.nn <- qread("brain-mets_ice_samples_nn_2016-08-04.tsv") %>% normalize_df_field_names;

patients.pkb <- qread("Brastianos_Carter_Supplementary_Table_1_7_22_15.tsv") %>% normalize_df_field_names %>%
	mutate(
		primary_frozen = ifelse(is.na(primary_frozen), 0, primary_frozen),
		bm_frozen = ifelse(is.na(bm_frozen), 0, bm_frozen)
	);

materials.nn <- qread("sample-info_missing-material_NN.tsv");

out.fname <- filename("sample-info", tag="wes", ext="tsv");
bams.fname <- filename("sample-bams", tag="wes", ext="vtr");
normal.fname <- filename("normal-bams", ext="vtr");

wes <- rbind(
	mutate(ice, wes_platform="ICE"),
	mutate(agilent, wes_platform="Agilent")
);

# derive sample_ids based on BAM file name

extract_name_from_bam_path <- function(x) {
	x %>%
		sub(".*/", "", .) %>%
		sub("\\.bam$", "", .)
}

# derive processing center based on BAM file name

extract_center_from_bam_path <- function(x) {
	ifelse(grepl("ccgd", x), "CCGD", "Broad")
}

wes$sample_id <- extract_name_from_bam_path(wes$clean_bam_file_capture);
wes$center <- extract_center_from_bam_path(wes$clean_bam_file_capture);

# put sample_id first
wes <- select(wes, sample_id, fh_sample_id:center);

# fix the sample_ids to make them more uniform and readily understandable

fix_bam_path <- function(x) {
	# remove tier3b subdirectory
	x %>%
	sub("tier3b/", "", .) %>%
	sub("picard_aggregation2", "picard_aggregation", .)
}

wes <- mutate(wes,
	clean_bam_file_capture = normalize_empty(
		fix_bam_path(as.character(clean_bam_file_capture)))
);

# fix sample names that match the dedup.cleaned pattern
dedup.idx <- grep("\\.dedup\\.cleaned", wes$sample_id);
# extract sample_id from fh_sample_id
wes[dedup.idx, "sample_id"] <- sub(".+SM-(PB\\d+-[NPM])", "\\1", wes[dedup.idx, "fh_sample_id"]);

# fix sample names that match the DAS pattern
das.idx <- grep("^DAS\\d+", wes$sample_id);
wes[das.idx, "sample_id"] <- sub(".+SM-(PB\\d+-[NPM])", "\\1", wes[das.idx, "fh_sample_id"]);

fix_generic_sample_id <- function(x) {
	x %>%
		# replace "RE" suffix
		sub("\\._?RE", "-RE", . ) %>%
		# fix "RE" suffix
		sub("(-RE)\\.(\\d+)", "\\1\\2", .) %>%
		# replace ".[ABCDE]" with "-[0-9]" pattern
		sub("(P?B?\\d+-[NPMT]+)\\.([A-Z0-9])", "\\1-\\2", .) %>%
		# replace "MT" designation with "M"
		sub("-MT", "-M", .)
}

insert_prefix <- function(x, prefix, sep="") {
	ifelse(is.na(x), NA, paste(prefix, x, sep=sep))
}

remove_na <- function(x) {
	x[!is.na(x)]
}


wes <- mutate(wes, sample_id = fix_generic_sample_id(sample_id));

# fix a specific sample_id
wes$sample_id[wes$sample_id == "PB00176-M-B"] <- "PB0176-M-B";

# set empty sample_id to NA
wes$sample_id[wes$sample_id == ""] <- NA;

# derive clinical_id
wes <- mutate(wes,
	clinical_id = 
		sample_id %>%
			sub("(((LS)|(BS)|(MM)|(MS))-\\d+)-.*", "\\1", .) %>%
			sub("(PB\\d+)-.*", "\\1", .)
);

wes.a <- wes %>%
	left_join(
		select(patients, clinical_id, gender, primary_cancer_type, primary_histotype),
		by = "clinical_id"
	) %>%
	left_join(
		select(samples, fh_pair_id = pair_id, sample_type, sample_description),
		by = "fh_pair_id"
	);

wes.a$sample_type <- as.character(wes.a$sample_type);
wes.a$sample_type[grepl("-N(-RE)?$", wes.a$sample_id)] <- "Normal";

# assume for now that all primary designated samples are primary tumours
# (NB  LS-005-P is a local recurrence of a primary tumour)

pkbid_to_sid <- function(x) {
	normalize_empty(as.character(x)) %>%
		strsplit("/", fixed=TRUE) %>% unlist %>%
		insert_prefix("PB") %>%
		fix_generic_sample_id
}

expand_field <- function(x) {
	normalize_empty(as.character(x)) %>%
		strsplit("/", fixed=TRUE) %>% unlist
}

omit_missing_key <- function(x, key.col=1) {
	key <- x[, key.col]
	x[!is.na(key), ]
}

options(stringsAsFactors=FALSE);

pkb.n.sids <- pkbid_to_sid(samples.pkb$normal_reference_id);
pkb.n.formats <- expand_field(samples.pkb$normal_source);
c(length(pkb.n.sids), length(pkb.n.formats))
stopifnot(length(pkb.n.sids) == length(pkb.n.formats));
pkb.n <- data.frame(
	sample_id = pkb.n.sids,
	specimen_material = pkb.n.formats,
	sample_type = "Normal"
) %>% omit_missing_key;

pkb.bm.sids <- pkbid_to_sid(samples.pkb$brain_specimen_reference_id);
pkb.bm.formats <- expand_field(samples.pkb$brain_specimen_tissue_format);
c(length(pkb.bm.sids), length(pkb.bm.formats))
stopifnot(length(pkb.bm.sids) == length(pkb.bm.formats));
pkb.bm <- data.frame(
	sample_id = pkb.bm.sids,
	specimen_material = pkb.bm.formats,
	sample_type = "Brain metastasis"
) %>% omit_missing_key;


pkb.pri.sids <- pkbid_to_sid(samples.pkb$primary_reference_id);
pkb.pri.formats <- expand_field(samples.pkb$primary_tissue_format);
c(length(pkb.pri.sids), length(pkb.pri.formats))
stopifnot(length(pkb.pri.sids) == length(pkb.pri.formats));
pkb.pri <- data.frame(
	sample_id = pkb.pri.sids,
	specimen_material = pkb.pri.formats,
	sample_type = "Primary"
) %>% omit_missing_key;

pkb.om.sids <- pkbid_to_sid(samples.pkb$other_specimen_reference_id);
pkb.om.formats <- expand_field(samples.pkb$other_specimen_tissue_format);
c(length(pkb.om.sids), length(pkb.om.formats))
stopifnot(length(pkb.om.sids) == length(pkb.om.formats));
pkb.om <- data.frame(
	sample_id = pkb.om.sids,
	specimen_material = pkb.om.formats,
	sample_type = NA
) %>% omit_missing_key;


pkb.n.sids <- pkbid_to_sid(samples.pkb$normal_reference_id) %>% remove_na;
pkb.bm.sids <- pkbid_to_sid(samples.pkb$brain_specimen_reference_id) %>% remove_na;
pkb.pri.sids <- pkbid_to_sid(samples.pkb$primary_reference_id) %>% remove_na;
pkb.m.sids <- pkbid_to_sid(samples.pkb$other_specimen_reference_id) %>% remove_na;

pkb.df <- rbind(pkb.n, pkb.bm, pkb.pri, pkb.om);

# join sample information from PKB's table
wes.a2 <- wes.a %>%
	left_join(pkb.df, by = "sample_id") %>%
	mutate(sample_type = ifelse(is.na(sample_type.y), sample_type.x, sample_type.y)) %>%
	select(-sample_type.x, -sample_type.y);

# join sample information from NN's table
samples.nn <- mutate(samples.nn, sample_id = fix_generic_sample_id(sample_id));
wes.a3 <- wes.a2 %>%
	left_join(select(samples.nn, -wes_platform, -notes), by = "sample_id") %>%
	mutate(specimen_material = ifelse(is.na(specimen_material.y), specimen_material.x, specimen_material.y)) %>%
	select(-specimen_material.x, -specimen_material.y) %>%
	mutate(sample_type = ifelse(is.na(sample_type.y), sample_type.x, sample_type.y)) %>%
	select(-sample_type.x, -sample_type.y) %>%
	mutate(sample_description = ifelse(is.na(sample_description.y), sample_description.x, sample_description.y)) %>%
	select(-sample_description.x, -sample_description.y);


# join material information from NN's table
wes.a4 <- wes.a3 %>%
	left_join(select(materials.nn, sample_id, specimen_material, accession),
		by = "sample_id") %>%
	mutate(
		specimen_material = ifelse(is.na(specimen_material.y), specimen_material.x, specimen_material.y)) %>%
	select(-specimen_material.x, -specimen_material.y) %>%
	mutate(
		accession = ifelse(is.na(accession.y), accession.x, accession.y)) %>%
	select(-accession.x, -accession.y);

wes.a.orig <- wes.a;
wes.a <- wes.a4;

#wes.a <- wes.a %>%
#	filter(! sample_id %in% blacklist);

# set x to y if y is not NA
set_if_not_na <- function(x, y) {
	ifelse(is.na(y), x, y)
}

idx_to_logical <- function(x, n) {
	y <- rep(FALSE, n);
	y[x] <- TRUE;
	y
}


# FIXME temporarily infer cancer type from ID
# Most BS-.* samples are breast cancers; however, there are some mislabelings
# All MM-.* samples should be multiple myeloma (unless there are mislabelings)
# All MS-.* samples should be melanoma (unless there are mislabelings)
idx <- grepl("^BS-.*", wes.a$sample_id) & is.na(wes.a$primary_cancer_type);
wes.a$primary_cancer_type[idx] <- "Breast cancer";
idx <- grepl("^MM-.*", wes.a$sample_id) & is.na(wes.a$primary_cancer_type);
wes.a$primary_cancer_type <- as.character(wes.a$primary_cancer_type);
wes.a$primary_histotype <- as.character(wes.a$primary_histotype);
wes.a$primary_cancer_type[idx] <- "Blood cancer";
wes.a$primary_histotype[idx] <- "Multiple myeloma";
idx <- grepl("^MS-.*", wes.a$sample_id) & is.na(wes.a$primary_cancer_type);
wes.a$primary_cancer_type[idx] <- "Skin cancer";
wes.a$primary_histotype[idx] <- "Melanoma";

wes.a$sample_type[wes.a$sample_type == "Lymph node"] <- "Extracranial metastasis";


filter(wes.a, is.na(primary_cancer_type))
# PB0027 has an unknown primary
# PB0238 is ambiguous: breast or endometrial cancer

message("Missing BAM files (", sum(is.na(wes.a$clean_bam_file_capture)), ")");
message(paste(wes.a$fh_sample_id[is.na(wes.a$clean_bam_file_capture)], collapse="\n"))
table(wes.a$primary_cancer_type, useNA="always");
filter(wes.a, is.na(primary_cancer_type))

#wes.a$sample_type[idx.pb] <- set_if_not_na(wes.a$sample_type[idx.pb], types.pb);

#wes.a[is.na(wes.a$sample_type) & idx_to_logical(idx.pb, nrow(wes.a)), ]

wes.lung <- filter(wes.a, primary_cancer_type == "Lung cancer");
table(wes.lung$sample_type, useNA="always");
unique(filter(wes.lung, sample_type == "Extracranial metastasis")$clinical_id) %>% length()
table(wes.lung$primary_histotype, useNA="always");
wes.lung %>%
	group_by(clinical_id) %>% summarize(primary_histotype=first(primary_histotype)) %>%
	group_by(primary_histotype) %>% summarize(count=n());
wes.lung %>%
	group_by(clinical_id) %>% summarize(primary_histotype=first(primary_histotype)) %>%
	group_by(remove_rare_levels(primary_histotype, 9)) %>% summarize(count=n());
wes.lung[is.na(wes.lung$sample_type), ];
length(unique(wes.lung$clinical_id));

table(wes.lung$primary_histotype)

wes.luad <- filter(wes.a, primary_histotype == "Lung adenocarcinoma");
nrow(wes.luad)
table(wes.luad$sample_type, useNA="always");
unique(filter(wes.luad, sample_type == "Extracranial metastasis")$clinical_id) %>% length()
wes.luad[is.na(wes.luad$sample_type), ];
wes.luad[is.na(wes.luad$gender), ];
length(unique(wes.luad$clinical_id))

wes.nsclc <- filter(wes.lung, primary_cancer_type != "Small cell lung carcinoma")
length(unique(wes.nsclc$clinical_id))

#pid.missing <- wes.luad[is.na(wes.luad$sample_type), "fh_pair_id"];
#idx <- pid.missing %in% absolute$sample;
#table(idx);
#absolute[absolute$sample %in% pid.missing, ]
#table(absolute[absolute$sample %in% pid.missing, "call_status"])

# samples with missing information
wes.a[which(is.na(wes.a$sample_type) & !is.na(wes.a$sample_id)), "sample_id"]

filter(wes.a, is.na(sample_type) & primary_cancer_type == "Lung cancer")
#filter(wes.a, is.na(sample_type))

stopifnot(!any(duplicated(wes.a$fh_sample_id)));
wes.a[ duplicated(wes.a$fh_sample_id) | duplicated(wes.a$fh_sample_id, fromLast=TRUE), ]

# it is unclear which PB0126-P-RE to use; for now, rename the later sample
wes.a$sample_id[wes.a$fh_sample_id == "PB-PB0126-Tumor-SM-7LYMY"] <- "PB0126-P-RE2";

stopifnot(!any(duplicated(remove_na(wes.a$sample_id))));
wes.a[ duplicated(wes.a$sample_id) | duplicated(wes.a$sample_id, fromLast=TRUE), ]

qwrite(wes.a, out.fname);


all.bams <- wes.a$clean_bam_file_capture;
normal.bams <- filter(wes.a, sample_type == "Normal")$clean_bam_file_capture;

qwrite(all.bams, bams.fname);
qwrite(normal.bams, normal.fname);


missing.material <- filter(wes.a, is.na(specimen_material)) %>% 
	select(sample_id, fh_sample_id, fh_pair_id, clinical_id, primary_cancer_type,
		specimen_material, accession);

qwrite(missing.material, "sample-info_missing-material.tsv");

