# Aggregate and organize all WES sample tables

library(magrittr);
library(io);
library(dplyr);

source("df_helper.R");

# PB0403-N is superceded by PB0403-N-RE
# BS-532-P was blacklisted
blacklist <- c("PB0403-N", "BS-532-P");

# missing BAM file:

ice <- qread("brain-mets_ice_samples.tsv");
stopifnot(!any(duplicated(ice$fh_sample_id)))

agilent <- qread("brain-mets_agilent-master_samples.tsv");
# remove duplicate samples (same bam file or deleted bam file)
# always remove the first sid, unless the second sid is unusual or the BAM file is not available
#fhsids.agilent.remove <- c("PB-PB0184-Tumor-SM-79ZXD", "PB-BS-532-Normal-SM-3WL9D", "PB-BS-532-Tumor-SM-3WL9E", "PB-BS-536-Tumor-SM-79ZX8", "PB-PB0190-Normal-SM-3WLA8", "PB-PB0238-Tumor-SM-3WLAJ", "PB0441-MT", "PB0441-N", "PB0441-P.A", "PB0441-P.B");
agilent <- filter(agilent, ! fh_sample_id %in% fhsids.agilent.remove);
stopifnot(!any(duplicated(agilent$fh_sample_id)))

patients <- qread("patient-info.rds");

samples <- qread("sample_types_2016-07-26.tsv");

#absolute <- qread("SLC.ABSOLUTE.table.07.24.2016.txt", type="tsv") %>% normalize_df_field_names;

samples.pkb <- qread("MasterList-PB-deidentified_ds_2016-08-04.tsv", quote="\"") %>% normalize_df_field_names;
samples.nn <- qread("brain-mets_ice_samples_nn_2016-08-04.tsv") %>% normalize_df_field_names;

out.fname <- filename("sample-info", tag="wes", ext="tsv");
normal.fname <- filename("normal-bams", ext="vtr");

wes <- rbind(
	mutate(ice, wes_platform="ICE"),
	mutate(agilent, wes_platform="Agilent")
);

# sample_ids were preliminarily assigned name based on BAM file name

# fix the sample_ids to make them more uniform and readily understandable

fix_bam_path <- function(x) {
	# remove tier3b subdirectory
	x %>%
	sub("tier3b/", "", .) %>%
	sub("picard_aggregation2", "picard_aggregation", .)
}

wes <- mutate(wes,
	sample_id = as.character(wes$sample_id),
	clean_bam_file_capture = normalize_empty(fix_bam_path(as.character(clean_bam_file_capture)))
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
			sub("([LB]S-\\d+)-.*", "\\1", .) %>%
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

# joint sample information from NN's table
samples.nn <- mutate(samples.nn, sample_id = fix_generic_sample_id(sample_id));
wes.a3 <- wes.a2 %>%
	left_join(select(samples.nn, -wes_platform, -notes), by = "sample_id") %>%
	mutate(specimen_material = ifelse(is.na(specimen_material.y), specimen_material.x, specimen_material.y)) %>%
	select(-specimen_material.x, -specimen_material.y) %>%
	mutate(sample_type = ifelse(is.na(sample_type.y), sample_type.x, sample_type.y)) %>%
	select(-sample_type.x, -sample_type.y) %>%
	mutate(sample_description = ifelse(is.na(sample_description.y), sample_description.x, sample_description.y)) %>%
	select(-sample_description.x, -sample_description.y);

wes.a.orig <- wes.a;
wes.a <- wes.a3;

wes.a <- wes.a %>%
	filter(! sample_id %in% blacklist);

# set x to y if y is not NA
set_if_not_na <- function(x, y) {
	ifelse(is.na(y), x, y)
}

idx_to_logical <- function(x, n) {
	y <- rep(FALSE, n);
	y[x] <- TRUE;
	y
}

message("Missing BAM files (", sum(is.na(wes.a$clean_bam_file_capture)), ")");
message(paste(wes.a$fh_sample_id[is.na(wes.a$clean_bam_file_capture)], collapse="\n"))

#wes.a$sample_type[idx.pb] <- set_if_not_na(wes.a$sample_type[idx.pb], types.pb);

#wes.a[is.na(wes.a$sample_type) & idx_to_logical(idx.pb, nrow(wes.a)), ]

wes.lung <- filter(wes.a, primary_cancer_type == "Lung cancer");
table(wes.lung$sample_type, useNA="always");
wes.lung[is.na(wes.lung$sample_type), ];

table(wes.lung$primary_histotype)

wes.luad <- filter(wes.a, primary_histotype == "Lung adenocarcinoma");
table(wes.luad$sample_type, useNA="always");
wes.luad[is.na(wes.luad$sample_type), ];

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
stopifnot(!any(duplicated(remove_na(wes.a$sample_id))));


qwrite(wes.a, out.fname);

normal.bams <- filter(wes.a, sample_type == "Normal")$clean_bam_file_capture;

qwrite(normal.bams, normal.fname);


