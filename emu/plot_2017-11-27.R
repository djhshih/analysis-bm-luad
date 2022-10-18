library(io);
library(ggplot2);
library(reshape2);
library(dplyr);

out.fname <- filename("tcga-pkb-luad");

pheno.fname <- "~/exigens/brain-mets/annot/sample-info_wes_stage2.tsv";
clin.fname <- "~/exigens/brain-mets/annot/patient-info.rds";

tcga.clin.fname <- "~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad.tsv";
tcga.pheno.fname <- "~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage2.tsv";


# TCGA smoking history classes:
# Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
# Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
# Current reformed smoker for greater than 15 years = 3
# Current reformed smoker for less than or equal to 15 years = 4
# Current reformed smoker, duration not specified = 5
# Smoking History not documented = 7
# (Source: https://groups.google.com/forum/#!topic/cbioportal/irEXZRj9Wh)

# for k = 6, signatures A and B are both correlated with smoking history
#spectra.fname <- "out/force/pkb_6_ml_spectra.txt";
#contribs.fname <- "out/force/pkb_6_map_activities.txt";

spectra.fname <- "out/force/pkb_5_ml_spectra.txt";
contribs.fname <- "out/force/pkb_5_map_activities.txt";
assigned.fname <- "out/force/pkb_5_assigned.txt";
samples.fname <- "in/pkb_all_samples.vtr";

#spectra.fname <- "out/force/pkb_4_ml_spectra.txt";
#contribs.fname <- "out/force/pkb_4_map_activities.txt";
#assigned.fname <- "out/force/pkb_4_assigned.txt";
#samples.fname <- "in/pkb_all_samples.vtr";

#spectra.fname <- "out/force/pkb_3_ml_spectra.txt";
#contribs.fname <- "out/force/pkb_3_map_activities.txt";
#assigned.fname <- "out/force/pkb_3_assigned.txt";
#samples.fname <- "in/pkb_all_samples.vtr";

pheno <- qread(pheno.fname);
clin <- qread(clin.fname);
tcga.clin <- qread(tcga.clin.fname);
tcga.pheno <- qread(tcga.pheno.fname);
samples <- qread(samples.fname);

substitutions <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
nucleotides <- c("A", "C", "G", "T");
context <- paste0(
	rep(nucleotides, each=length(nucleotides)),
	rep("X", times=length(nucleotides) * length(nucleotides)),
	rep(nucleotides, times=length(nucleotides))
);

read_matrix <- function(fname) {
	as.matrix(read.table(fname, header=FALSE));
}

spectra <- read_matrix(spectra.fname);
colnames(spectra) <- NULL;
rownames(spectra) <- LETTERS[1:nrow(spectra)];

spectra.df <- melt(spectra);

spectra_matrix_to_df <- function(x) {
	xt <- t(x);
	nsig <- ncol(xt);
	data.frame(
		signature = rep(colnames(xt), each=96),
		context = rep(context, times=6*nsig),
		substitution = rep(substitutions, each=16),
		activity = c(xt)
	)
}

spectra.df <- spectra_matrix_to_df(spectra);

qdraw(
	ggplot(spectra.df, aes(x=context, y=activity, fill=substitution)) +
		facet_grid(signature ~ substitution) +
		geom_bar(stat="identity") + theme_bw()
	,
	width = 8,
	height = 6,
	file = insert(out.fname, tag="spectra", ext="pdf")
);


contribs <- read_matrix(contribs.fname);
colnames(contribs) <- LETTERS[1:ncol(contribs)];
rownames(contribs) <- samples;
contribs <- contribs / rowSums(contribs);

contribs.pkb <- contribs[rownames(contribs) %in% pheno$sample_id, ];

contribs.pkb.mdf <- melt(contribs.pkb);
names(contribs.pkb.mdf) <- c("sample_id", "signature", "contribution");

contribs.pkb.df <- mutate(contribs.pkb.mdf,
	clinical_id = pheno$clinical_id[match(sample_id, pheno$sample_id)]
) %>% left_join(
	select(clin, clinical_id, smoker, smoking_hx, age_at_primary_dx),
	by = "clinical_id"
);

ggplot(filter(contribs.pkb.df), aes(x=smoker, y=contribution, colour=signature)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw() + coord_flip() + facet_grid(signature ~ .)

ggplot(filter(contribs.pkb.df), aes(x=contribution, fill=smoker)) +
	geom_density(alpha=0.3) + theme_bw() + facet_grid(signature ~ ., scale="free_y")

ggplot(filter(contribs.pkb.df), aes(x=smoking_hx, y=contribution)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y")

ggplot(filter(contribs.pkb.df), aes(x=age_at_primary_dx, y=contribution)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y")


pkb.samples.exposed <- filter(contribs.pkb.df, signature == "E", contribution > 0.2)$sample_id;
pkb.patients <- as.character(unique(pheno$clinical_id[match(rownames(contribs.pkb), pheno$sample_id)]));
pkb.patients.exposed <- as.character(unique(pheno$clinical_id[pheno$sample_id %in% pkb.samples.exposed]));
length(pkb.patients.exposed) / length(pkb.patients)

d <- data.frame(
	clinical_id = pkb.patients,
	smoking_signature_high = pkb.patients %in% pkb.patients.exposed,
	smoker = clin$smoker[match(pkb.patients, clin$clinical_id)]
);
ct <- with(d, table(smoker, smoking_signature_high));
fisher.test(ct)


contribs.tcga <- contribs[grep("TCGA", rownames(contribs)), ];

contribs.mdf <- melt(contribs.tcga);
names(contribs.mdf) <- c("sample_id", "signature", "contribution");

tcga_sample_id_to_clinical_id <- function(x) {
	sub(".*(TCGA-[^-]+-[^-]+).*", "\\1", x)
}

tcga.clin.smoking <- mutate(tcga.clin,
		age_at_primary_dx = age_at_initial_pathologic_diagnosis,
		smoking_hx_class = as.factor(tobacco_smoking_history),
		smoking_hx = ifelse(!is.na(smoking_hx_class) & smoking_hx_class == "1", 0, number_pack_years_smoked),
		smoker = smoking_hx > 0
	) %>%
	select(
		clinical_id, smoking_hx_class, smoking_hx, smoker, age_at_primary_dx
	);

contribs.df <- mutate(contribs.mdf,
		clinical_id = tcga_sample_id_to_clinical_id(sample_id)
	) %>%
	left_join(
		tcga.clin.smoking,
		by = "clinical_id"
	);

with(tcga.clin.smoking, table(smoker, smoking_hx_class, useNA="always"))
with(tcga.clin.smoking, table(is.na(smoking_hx_class), is.na(smoking_hx)))

ggplot(tcga.clin.smoking, aes(x=smoking_hx_class, y=smoking_hx)) +
	geom_jitter(width=0.2) + theme_bw()


ggplot(contribs.df, aes(x=smoker, y=contribution, colour=signature)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw() + coord_flip() + facet_grid(signature ~ .)

ggplot(contribs.df, aes(x=smoking_hx_class, y=contribution, colour=signature)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw() + coord_flip() + facet_grid(signature ~ .)

# drop class 5, since there are too few data points
d <- select(contribs.df, contribution, smoking_hx_class, signature) %>%
	mutate(smoking_hx_class = factor(smoking_hx_class, levels=1:4));
d <- na.omit(d);

ggplot(contribs.df, aes(x=contribution, fill=smoking_hx_class)) +
	geom_density(alpha=0.3) + theme_bw() + facet_grid(signature ~ .)

quantile(filter(contribs.df, signature == "E", smoking_hx_class == 1)$contribution, seq(0, 1, 0.05))

ggplot(contribs.df, aes(x=smoking_hx, y=contribution)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y")

ggplot(contribs.df, aes(x=age_at_primary_dx, y=contribution)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y")

#ggplot(filter(contribs.df, signature=="D"), aes(x=age_at_primary_dx, y=contribution)) +
#	geom_point() + theme_bw() + ylim(0, 1e-04)

tcga.samples.exposed <- filter(contribs.df, signature == "E", contribution > 0.2)$sample_id;
tcga.patients <- as.character(unique(tcga_sample_id_to_clinical_id(rownames(contribs.tcga))));
tcga.patients.exposed <- as.character(unique(tcga_sample_id_to_clinical_id(tcga.samples.exposed)));
length(tcga.patients.exposed) / length(tcga.patients)

tcga.d <- data.frame(
	clinical_id = tcga.patients,
	smoking_signature_contrib = factor(tcga.patients %in% tcga.patients.exposed, levels=c(FALSE, TRUE), labels=c("Low", "High")),
	smoker = factor(tcga.clin.smoking$smoking_hx_class[match(tcga.patients, tcga.clin.smoking$clinical_id)] > 1, levels=c(FALSE, TRUE), labels=c("Never", "Ever"))
);
tcga.ct <- with(tcga.d, table(smoker, smoking_signature_contrib));
tcga.ct
fisher.test(tcga.ct)


####

assigned <- read_matrix(assigned.fname);
colnames(assigned) <- LETTERS[1:ncol(assigned)];
rownames(assigned) <- samples;
assigned.n <- assigned / rowSums(assigned);

plot(c(contribs), c(assigned.n))
abline(a=0, b=1)

# signature D normalized assigned value differs systematically from
# contributions... possible bug in EMu?
plot(contribs[, "D"], assigned.n[, "D"])
abline(a=0, b=1)

plot(contribs[, "E"], assigned.n[, "E"])
abline(a=0, b=1)


assigned.mdf <- melt(assigned[grep("TCGA", rownames(assigned)), ]);
names(assigned.mdf) <- c("sample_id", "signature", "assigned");

assigned.df <- mutate(assigned.mdf,
		clinical_id = tcga_sample_id_to_clinical_id(sample_id)
	) %>%
	left_join(
		tcga.clin.smoking,
		by = "clinical_id"
	);

ggplot(filter(assigned.df), aes(x=smoker, y=assigned, colour=signature)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw() + coord_flip() + facet_grid(signature ~ .) +
	scale_y_log10()

ggplot(filter(assigned.df), aes(x=smoking_hx_class, y=assigned, colour=signature)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw() + coord_flip() + facet_grid(signature ~ .) +
	scale_y_log10()

g <- ggplot(filter(assigned.df), aes(x=smoking_hx, y=assigned)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y");
g
g + scale_y_log10()


g <- ggplot(filter(assigned.df), aes(x=age_at_primary_dx, y=assigned)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y");
g
g + scale_y_log10()
g + ylim(0, 100)


####


# FIXME for mut, want individual mmaf
#mut <- qread("../comut/pkb-luad_brain-mets_black-f_20171031T112959.pset.mmaf", type="tsv", quote="");
tcga.mut <- qread("../comut/tcga-luad/tcga-luad_pass_black-f_20171113T105434.pset.mmaf", type="tsv", quote="");

#mut$sample_id <- pheno$sample_id[match(mut$Tumor_Sample_Barcode, pheno$fh_sample_id)];
tcga.mut$sample_id <- tcga.pheno$sample_id[match(tcga.mut$Tumor_Sample_Barcode, tcga.pheno$fh_sample_id)];

tcga.df <- data.frame(
	sample_id = rownames(contribs.tcga),
	clinical_id = tcga_sample_id_to_clinical_id(rownames(contribs.tcga)),
	contribs.tcga
);
colnames(tcga.df)[-(1:2)] <- paste0("mutsig_k5_", tolower(colnames(tcga.df[-(1:2)])));

lof.muts <- c("Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "In_Frame_Ins");

tcga.df$brca_mut <- FALSE;
brca.genes <- c("BRCA1", "BRCA2");
brca.maf <- filter(tcga.mut, Hugo_Symbol %in% brca.genes, Variant_Classification %in% lof.muts);
brca.mut.samples <- as.character(brca.maf$sample_id);
tcga.df$brca_mut[tcga.df$sample_id %in% brca.mut.samples] <- TRUE;

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_b)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_c)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_a)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_b)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_c)) +
	geom_boxplot() + theme_bw()


mmr.genes <- c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2");
tcga.df$mmr_mut <- FALSE;
mmr.maf <- filter(tcga.mut, Hugo_Symbol %in% mmr.genes, Variant_Classification %in% lof.muts);
mmr.mut.samples <- as.character(mmr.maf$sample_id);
tcga.df$mmr_mut[tcga.df$sample_id %in% mmr.mut.samples] <- TRUE;

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_b)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_c)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_a)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_b)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_c)) +
	geom_boxplot() + theme_bw()

wilcox.test(mutsig_k5_c ~ mmr_mut, data = tcga.df)


apobec.genes <- c("APOBEC1", "APOBEC3A", "APOBEC3B");
tcga.df$apobec_mut <- FALSE;
apobec.maf <- filter(tcga.mut, Hugo_Symbol %in% apobec.genes, Variant_Classification %in% lof.muts);
apobec.mut.samples <- as.character(apobec.maf$sample_id);
tcga.df$apobec_mut[tcga.df$sample_id %in% apobec.mut.samples] <- TRUE;

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_b)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_c)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_a)) +
	geom_boxplot()+ theme_bw()
wilcox.test(mutsig_k5_a ~ apobec_mut, data = tcga.df)

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_b)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_c)) +
	geom_boxplot() + theme_bw()


tcga.df <- left_join(tcga.df, tcga.clin.smoking, by="clinical_id");

ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_b)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_c)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

#ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_e)) +
#	geom_jitter() + theme_bw()

ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_a)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_b)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_c)) +
	geom_boxplot() + theme_bw()

#ggplot(filter(tcga.df), aes(x=as.factor(smoking_hx_class), y=mutsig_k5_e)) +
#	geom_boxplot() + theme_bw()

