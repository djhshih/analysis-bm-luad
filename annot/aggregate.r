# Aggregate and organize all sample information tables

library(magrittr);
library(io);
library(dplyr);

source("df_helper.R");

clin.n <- qread("patient-info.tsv");
phylogic <- qread("all_cases_phylogic_SIF_09_01_2015.tsv")
rna <- qread("RNA_iset.tsv");
unpaired <- qread("unpaired_met_SIF_07.09.2015.tsv");

phylogic.n <- phylogic %>%
	normalize_df_field_names %>%
	rename_df_field("sample", "sample_description") %>%
	rename_df_field("sample_type", "met_type");

unpaired.n <- unpaired[, -1] %>%
	normalize_df_field_names %>%
	rename_df_field("met_sid", "pair_id");

rna.n <- rna[, 2:3] %>%
	normalize_df_field_names %>%
	rename_df_field("individual_id", "rnaseq_id");

merged <- phylogic.n %>%
	#merge(unpaired.n[! unpaired.n$pair_id %in% phylogic.n$pair_id, ], all=TRUE) %>%
	merge(rna.n, all=TRUE);

# add missing samples from clin.n into merged
merged$met_type <- as.character(merged$met_type);
idx <- ! clin.n$patient_id %in% merged$patient_id;
for (i in which(idx)) {
	xi <- merged[1, ];
	xi[1, ] <- NA;
	if (!is.na(clin.n$primary_pair_id[i]) && ! clin.n$primary_pair_id[i] %in% merged$pair_id) {
		xi$pair_id = clin.n$primary_pair_id[i];
		xi$patient_id = clin.n$patient_id[i];
		xi$met_type = "Primary";
		merged <- rbind(merged, xi);
	}
	if (!is.na(clin.n$met_pair_id[i]) && ! clin.n$met_pair_id[i] %in% merged$pair_id) {
		xi$pair_id = clin.n$met_pair_id[i];
		xi$patient_id = clin.n$patient_id[i];
		# Problem: do not know the metastatic type of these missing samples
		xi$met_type = "M";
		merged <- rbind(merged, xi);
	}
}


pair_ids.missing <- as.character(merged$pair_id[!is.na(merged$rnaseq_id) & is.na(merged$met_type)]);

merged$patient_id <- as.character(merged$patient_id);

i <- merged$pair_id == "PB-LS-004-TP-NT-SM-7M185-SM-7M186";
merged$patient_id[i] <- "PB-LS-004";
merged$met_type[i] <- "Primary";

i <- merged$pair_id == "PB0067-TP-NT-SM-3WL9V-SM-PB0067-N";
merged$patient_id[i] <- "0067-MT";
merged$met_type[i] <- "Primary";

i <- merged$pair_id == "PB0067-TP-NT-SM-3WL9W-SM-PB0067-N";
merged$patient_id[i] <- "0067-MT";
merged$met_type[i] <- "Primary";

i <- merged$pair_id == "PB0300-TM-NT-SM-2XUEX-SM-2XUEY";
merged$patient_id[i] <- "0300-MT";
merged$met_type[i] <- "BM"; # ??

i <- merged$pair_id == "PB0387-TP-NT-SM-2XUIF-SM-2XUIE";
merged$patient_id[i] <- "PB0387";
merged$met_type[i] <- "BM";

#histology.df <- clin.n[match(d$patient_id, clin.n$patient_id), c("patient_id", "histology_group")]

# hideous hack to derive tumour ids (used in exome-seq output)
merged$tumor_id <-
	gsub("-SM-[^-]+-N$", "",
	gsub("-SM-[^-]+$", "",
	gsub("-T.-NT-", "-Tumor-", merged$pair_id)));

#exome.seq <- qread("../exome-seq/luad-bm_sig-snvs_exome-seq_mutsigcv.tsv", quote="");
#stopifnot(all(!is.na(match(exome.seq$Tumor_Sample_Barcode, merged$tumor_id))))

#cna <- qread("../cn/combined_gene_corrected_CN.txt", type="mtx");
#idx <- is.na(match(colnames(cna), merged$pair_id));
#colnames(cna)[idx]

#idx <- merged.rnaseq$pair_id%in% colnames(cna)
#merged.rnaseq$pair_id[!idx]

merged.rnaseq <- merged[!is.na(merged$rnaseq_id), ];

# mets inherit histotype from primary tumour...
merged$histotype <- clin.n[match(merged$patient_id, clin.n$patient_id), "histotype"];

stopifnot(all(!duplicated(merged$pair_id)));

qwrite(merged, filename("sample-info", ext="tsv"));
#qwrite(clin.n, filename("patient-info", ext="tsv"));
qwrite(merged[!is.na(merged$rnaseq_id), ], filename("sample-info_rnaseq", ext="tsv"));

