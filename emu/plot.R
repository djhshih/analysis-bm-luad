library(io);
library(ggplot2);
library(reshape2);
library(dplyr);

source("~/exigens/brain-mets/common/colours.R");

out.fname <- filename("tcga-pkb-luad");
pkb.fname <- filename("pkb-luad");
tcga.fname <- filename("tcga-luad");

pheno.fname <- "~/exigens/brain-mets/annot/sample-info_wes_stage2.tsv";
clin.fname <- "~/exigens/brain-mets/annot/patient-info.rds";

tcga.clin.fname <- "~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad.tsv";
tcga.pheno.fname <- "~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage2.tsv";

# COSMIC
# LUAD
# Signatures observed:
# 1 (5'methylcytosine deamination)
# 2 (APOBEC)
# 4 (Smoking)
# 5 (Unknown)
# 6 (MMR)
# 13 (APOBEC)
# 17 (Unknown)

# TCGA smoking history classes:
# Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
# Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
# Current reformed smoker for greater than 15 years = 3
# Current reformed smoker for less than or equal to 15 years = 4
# Current reformed smoker, duration not specified = 5
# Smoking History not documented = 7
# (Source: https://groups.google.com/forum/#!topic/cbioportal/irEXZRj9Wh)

samples.fname <- "in/pkb_all_samples.vtr";
cmut.fname <- "in/pkb_all.mut.matrix";
opp.fname <- "in/pkb_all.opp.matrix";

# for k = 6, signatures C and E are both correlated with smoking history
#spectra.fname <- "out/force/pkb_6_ml_spectra.txt";
#contribs.fname <- "out/force/pkb_6_map_activities.txt";

#spectra.fname <- "out/force/pkb_5_ml_spectra.txt";
#contribs.fname <- "out/force/pkb_5_map_activities.txt";
#assigned.fname <- "out/force/pkb_5_assigned.txt";

#spectra.fname <- "out/force/pkb_4_ml_spectra.txt";
#contribs.fname <- "out/force/pkb_4_map_activities.txt";
#assigned.fname <- "out/force/pkb_4_assigned.txt";
#k <- 4;

spectra.fname <- "out/force/pkb_3_ml_spectra.txt";
contribs.fname <- "out/force/pkb_3_map_activities.txt";
#assigned.fname <- "out/force/pkb_3_assigned.txt";
k <- 3;


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

entropy <- function(p) {
	x <- p * log(p);
	x[p == 0] <- 0;
	-sum(x, na.rm=TRUE)
}

order_signatures_old <- function(spectra) {
	starts <- seq(1, 96, by=16);
	ends <- seq(16, 96, by=16);
	sums <- mapply(
		function(start, end) {
			rowSums(spectra[, start:end])
		},
		starts,
		ends
	);
	maxs <- apply(sums, 1, max);
	order(maxs, decreasing=TRUE)
}

order_signatures <- function(spectra) {
	apply(spectra, 1, entropy)
	order(apply(spectra, 1, entropy))
}

spectra <- read_matrix(spectra.fname);
sig.idx <- order_signatures(spectra);
spectra <- spectra[sig.idx, ];
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
	ggplot(spectra.df, aes(x=context, y=activity, fill=signature)) +
		facet_grid(signature ~ substitution) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(type="qual", palette="Set1")
	,
	width = 8,
	height = 6,
	file = insert(out.fname, tag="spectra", ext="pdf")
);


contribs <- read_matrix(contribs.fname);
contribs <- contribs[, sig.idx];
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

qdraw(
	ggplot(filter(contribs.pkb.df, !is.na(smoker)), aes(x=smoker, y=contribution, colour=signature)) +
		geom_jitter(width=0.2, alpha=0.3) + theme_bw() + facet_grid(. ~ signature) + ylim(0, 1) +
		scale_colour_brewer(type="qual", palette="Set1")
	,
	file = insert(pkb.fname, c("mutsig", "cor", "smoker"), ext="pdf")
);

qdraw(
	ggplot(filter(contribs.pkb.df, !is.na(smoker)), aes(x=contribution, fill=smoker)) +
		geom_density(alpha=0.3, bw=0.05) + theme_bw() + facet_grid(signature ~ ., scale="free_y") + xlim(0, 1)
	,
	file = insert(pkb.fname, c("mutsig", "cor", "smoker", "density"), ext="pdf")
);


logit <- function(x) {
	log(x) - log(1 - x)
}

qdraw(
	ggplot(filter(contribs.pkb.df), aes(x=smoking_hx, y=logit(contribution))) +
		geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y") + 
		geom_smooth(method = "lm")
	,
	width = 9, height = 3,
	file = insert(pkb.fname, c("mutsig", "lm", "smoking-hx"), ext="pdf")
);

qdraw(
	ggplot(filter(contribs.pkb.df), aes(x=age_at_primary_dx, y=logit(contribution))) +
		geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y") +
		geom_smooth(method = "lm")
	,
	width = 9, height = 3,
	file = insert(pkb.fname, c("mutsig", "lm", "age"), ext="pdf")
);

ggplot(filter(contribs.pkb.df), aes(x=age_at_primary_dx, y=contribution)) +
	geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y") +
	geom_smooth(method = "glm", method.args=list(family="binomial"))


pkb.samples.exposed <- filter(contribs.pkb.df, signature == "C", contribution > 0.2)$sample_id;
pkb.patients <- as.character(unique(pheno$clinical_id[match(rownames(contribs.pkb), pheno$sample_id)]));
pkb.patients.exposed <- as.character(unique(pheno$clinical_id[pheno$sample_id %in% pkb.samples.exposed]));
length(pkb.patients.exposed) / length(pkb.patients)

d <- data.frame(
	clinical_id = pkb.patients,
	smoking_signature_high = factor(pkb.patients %in% pkb.patients.exposed, levels=c(FALSE, TRUE), labels=c("Low", "High")),
	smoker = clin$smoker[match(pkb.patients, clin$clinical_id)]
);
ct <- with(d, table(smoker, smoking_signature_high));
ct
fisher.test(ct)
(ct[1,1] + ct[2,2]) / sum(ct)

qwrite(d, insert(pkb.fname, "smoking", ext="tsv"));


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
		smoker = factor(smoking_hx > 0, levels=c(FALSE, TRUE), labels=c("Never", "Ever"))
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


qdraw(
	ggplot(filter(contribs.df, !is.na(smoker)), aes(x=smoker, y=contribution, colour=signature)) +
		geom_jitter(width=0.2, alpha=0.3) + theme_bw() + facet_grid(. ~ signature) + ylim(0, 1) +
		scale_colour_brewer(type="qual", palette="Set1")
	,
	file = insert(tcga.fname, c("mutsig", "cor", "smoker"), ext="pdf")
);

qdraw(
	ggplot(filter(contribs.df, !is.na(smoker)), aes(x=contribution, fill=smoker)) +
		geom_density(alpha=0.3, bw=0.05) + theme_bw() + facet_grid(signature ~ ., scale="free_y")
	,
	file = insert(tcga.fname, c("mutsig", "cor", "smoker", "density"), ext="pdf")
);

qdraw(
	ggplot(filter(contribs.df), aes(x=smoking_hx, y=logit(contribution))) +
		geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y") + 
		geom_smooth(method = "lm")
	,
	height = 3, width = 9,
	file = insert(tcga.fname, c("mutsig", "lm", "smoking-hx"), ext="pdf")
);

qdraw(
	ggplot(filter(contribs.df), aes(x=age_at_primary_dx, y=logit(contribution))) +
		geom_point(alpha=0.3) + theme_bw() + facet_wrap(~signature, scale="free_y") +
		geom_smooth(method = "lm")
	,
	height = 3, width = 9,
	file = insert(tcga.fname, c("mutsig", "lm", "age"), ext="pdf")
);


#ggplot(filter(contribs.df, signature=="D"), aes(x=age_at_primary_dx, y=contribution)) +
#	geom_point() + theme_bw() + ylim(0, 1e-04)

tcga.samples.exposed <- filter(contribs.df, signature == "C", contribution > 0.2)$sample_id;
tcga.patients <- as.character(unique(tcga_sample_id_to_clinical_id(rownames(contribs.tcga))));
tcga.patients.exposed <- as.character(unique(tcga_sample_id_to_clinical_id(tcga.samples.exposed)));
length(tcga.patients.exposed) / length(tcga.patients)

tcga.d <- data.frame(
	clinical_id = tcga.patients,
	smoking_signature_high = factor(tcga.patients %in% tcga.patients.exposed, levels=c(FALSE, TRUE), labels=c("Low", "High")),
	smoker = factor(tcga.clin.smoking$smoking_hx_class[match(tcga.patients, tcga.clin.smoking$clinical_id)] != "1", levels=c(FALSE, TRUE), labels=c("Never", "Ever"))
);
tcga.ct <- with(tcga.d, table(smoker, smoking_signature_high));
tcga.ct
fisher.test(tcga.ct)

qwrite(tcga.d, insert(tcga.fname, "smoking", ext="tsv"));


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

tcga.cn <- qread("../comut/tcga-luad/tcga-luad_pass_absolute-1-4_gene_corrected_CN.txt", type="tsv");

lof.muts <- c("Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "In_Frame_Ins");
gof.muts <- c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins");

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

ggplot(filter(tcga.df), aes(x=brca_mut, y=mutsig_k5_d)) +
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

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_d)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_a)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_b)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_c)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=mmr_mut, y=mutsig_k5_d)) +
	geom_boxplot() + theme_bw()

wilcox.test(mutsig_k5_c ~ mmr_mut, data = tcga.df)


apobec.genes <- c("APOBEC1", "APOBEC3A", "APOBEC3B");
tcga.df$apobec_mut <- FALSE;
apobec.maf <- filter(tcga.mut, Hugo_Symbol %in% apobec.genes, Variant_Classification %in% gof.muts);
apobec.mut.samples <- as.character(apobec.maf$sample_id);
tcga.df$apobec_mut[tcga.df$sample_id %in% apobec.mut.samples] <- TRUE;

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_b)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_c)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_d)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_a)) +
	geom_boxplot()+ theme_bw()
wilcox.test(mutsig_k5_a ~ apobec_mut, data = tcga.df)

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_b)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_c)) +
	geom_boxplot() + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_mut, y=mutsig_k5_d)) +
	geom_boxplot() + theme_bw()

amp.cut <- 8;
apobec.amp.samples <- names(which(colSums(tcga.cn[apobec.genes, ] > amp.cut) > 1));
tcga.df$apobec_amp <- FALSE;
tcga.df$apobec_amp[tcga.df$sample_id %in% apobec.amp.samples] <- TRUE;

filter(tcga.df, apobec_amp)
tcga.cn[apobec.genes, apobec.amp.samples]

ggplot(filter(tcga.df), aes(x=apobec_amp, y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=apobec_amp, y=mutsig_k5_a)) +
	geom_boxplot(width=0.2, alpha=0.3) + theme_bw()

# KRAS is close to APOBEC...
kras.amp.samples <- colnames(tcga.cn)[which(tcga.cn["KRAS", ] > amp.cut)];
tcga.df$kras_amp <- FALSE;
tcga.df$kras_amp[tcga.df$sample_id %in% kras.amp.samples] <- TRUE;

with(tcga.df, table(kras_amp, apobec_amp))
# but amplifications do not overlap

filter(tcga.df, kras_amp)

ggplot(filter(tcga.df), aes(x=kras_amp, y=mutsig_k5_a)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()
ggplot(filter(tcga.df), aes(x=kras_amp, y=mutsig_k5_a)) +
	geom_boxplot(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=kras_amp, y=mutsig_k5_b)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()
ggplot(filter(tcga.df), aes(x=kras_amp, y=mutsig_k5_b)) +
	geom_boxplot(width=0.2, alpha=0.3) + theme_bw()

ggplot(filter(tcga.df), aes(x=kras_amp, y=mutsig_k5_c)) +
	geom_jitter(width=0.2, alpha=0.3) + theme_bw()
ggplot(filter(tcga.df), aes(x=kras_amp, y=mutsig_k5_c)) +
	geom_boxplot(width=0.2, alpha=0.3) + theme_bw()


####

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


####

#' @param A  an N by P matrix of barycentric coordinatess of N points in P
#'           dimensional space
#' @param V  a P by 2 matrix of 2-dimensional cartesian coordinates of the P
#'           vertices
#' @return a N by 2 matrix of points in 2-dimensional cartesian space
bary_to_cartesian2d <- function(A, V) {
	stopifnot(ncol(V) == 2);
	R <- A %*% V;
	colnames(R) <- c("x", "y");
	R
}

# vertices are ordered counterclockwise:
# {(1, -1), (1, 1), (-1, 1), (-1, -1)}
#V <- matrix(c(1, -1, 1, 1, -1, 1, -1, -1), byrow=TRUE, ncol=2);
#rownames(V) <- colnames(contribs);

#V <- matrix(c(0, 0, 1, 0, 0, 1), byrow=TRUE, ncol=2);
#rownames(V) <- colnames(contribs);

V <- matrix(c(0, 0, 1, 0, 0.5, sin(pi/3)), byrow=TRUE, ncol=2);
rownames(V) <- colnames(contribs);
colnames(V) <- c("x", "y");

clin.smoking.df <- rbind(
	select(clin, clinical_id, smoker),
	select(tcga.clin.smoking, clinical_id, smoker)
);

contribs.coord <- data.frame(
	sample_id = rownames(contribs),
	bary_to_cartesian2d(contribs, V)
) %>% mutate(
	cohort = factor(grepl("TCGA", sample_id), levels=c(TRUE, FALSE), labels=cohorts),
	clinical_id = ifelse(cohort == "TCGA", tcga_sample_id_to_clinical_id(sample_id), as.character(pheno$clinical_id[match(sample_id, pheno$sample_id)]))
) %>% left_join(
	clin.smoking.df, by="clinical_id"
) %>% left_join(
	select(pheno, sample_id, sample_type, specimen_material, wes_platform), by="sample_id"
);

contribs.coord$wes_platform[is.na(contribs.coord$wes_platform)] <- "Agilent";
contribs.coord$sample_type[is.na(contribs.coord$sample_type)] <- "Primary";
contribs.coord$specimen_material[is.na(contribs.coord$specimen_material)] <- "FF";

theme_blank <- function() {
      theme(
	    axis.line=element_blank(),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
	    axis.ticks=element_blank(),
	    axis.title.x=element_blank(),
	    axis.title.y=element_blank(),
	    #legend.position="none",
	    panel.background=element_blank(),
	    panel.border=element_blank(),
	    panel.grid.major=element_blank(),
	    panel.grid.minor=element_blank(),
	    plot.background=element_blank()
      )
}	


library(RColorBrewer)

baryplot <- function(d, V, mapping, alpha=0.6) {
	pal <- brewer.pal(3, "Set1");
	v.coord <- data.frame(V, signature=rownames(V));
	v.path <- as.data.frame(rbind(V, V[1,]));
	ggplot(d, aes(x=x, y=y)) +
		geom_point(mapping, alpha=alpha) +
		theme_blank() + xlim(-0.05, 1.05) + ylim(-0.05, 1.05) +
		geom_path(aes(x=x, y=y), data=v.path, colour="black") +
		geom_point(aes(x=x, y=y),colour=pal, size=2, data=v.coord) +
		annotate("text", x=c(0, 1, 0.5), y=c(0, 0, sin(pi/3)), label=c("APOBEC", "5-mC deamination", "Smoking"), colour=pal, vjust=c(2, 2, -1), hjust=c(0, 1, 0.5))
}

qdraw(
	baryplot(contribs.coord, V, aes(colour=cohort, shape=cohort)) +
	    scale_colour_manual(values=col.cohorts)
	,
	width = 5.3,
	file = insert(out.fname, c("mutsig", "bary", "cohort"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord, V, aes(colour=smoker, shape=smoker)) +
	    scale_colour_manual(values=col.smoker)
	,
	width = 5.3,
	file = insert(out.fname, c("mutsig", "bary", "smoker"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord, V, aes(colour=specimen_material, shape=specimen_material)) + 
		labs(colour="specimen", shape="specimen") 
	,
	width = 5.5,
	file = insert(out.fname, c("mutsig", "bary", "specimen-material"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord, V, aes(colour=sample_type, shape=sample_type))
	,
	width = 6.3,
	file = insert(out.fname, c("mutsig", "bary", "sample-type"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord, V, aes(colour=wes_platform, shape=wes_platform))
	,
	width = 5.5,
	file = insert(out.fname, c("mutsig", "bary", "wes-platform"), ext="pdf")
);


contribs.coord.pkb <- filter(contribs.coord, cohort == "Present");

qdraw(
	baryplot(contribs.coord.pkb, V, aes(colour=smoker, shape=smoker))
	,
	width = 5.3,
	file = insert(pkb.fname, c("mutsig", "bary", "smoker"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord.pkb, V, aes(colour=specimen_material, shape=specimen_material)) + 
		labs(colour="specimen", shape="specimen") 
	,
	width = 5.5,
	file = insert(pkb.fname, c("mutsig", "bary", "specimen-material"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord.pkb, V, aes(colour=sample_type, shape=sample_type))
	,
	width = 6.3,
	file = insert(pkb.fname, c("mutsig", "bary", "sample-type"), ext="pdf")
);

# add arrows to connect primary sample to other samples from the same patient
sample_arrows <- function(d) {
	ds <- split(d, d$clinical_id);
	e <- ds[[1]];

	arrows.df <- do.call(rbind, lapply(ds,
		function(e) {
			idx <- which(e$sample_type == "Primary");
			if (length(idx) > 0) {
				idx.other <- setdiff(1:nrow(e), idx);
				n.other <- length(idx.other);
				if (n.other > 0) {
					x <- rep(e$x[idx], n.other);
					y <- rep(e$y[idx], n.other);
					xend <- e$x[idx.other];
					yend <- e$y[idx.other];
					data.frame(x, y, xend, yend)
				} else {
					NULL
				}
			} else {
				NULL
			}
		}
	));

	geom_segment(
		aes(x=x, xend=xend, y=y, yend=yend),
		data = arrows.df,
		arrow = arrow(length=unit(6, "points")),
		alpha = 0.3
	)
}

qdraw(
	baryplot(contribs.coord.pkb, V, aes(colour=sample_type, shape=sample_type), alpha=0.8) +
		sample_arrows(contribs.coord.pkb)
	,
	height = 7.5,
	width = 9.45,
	file = insert(pkb.fname, c("mutsig", "bary", "sample-type", "arrows"), ext="pdf")
);

qdraw(
	baryplot(contribs.coord.pkb, V, aes(colour=wes_platform, shape=wes_platform))
	,
	width = 5.5,
	file = insert(pkb.fname, c("mutsig", "bary", "wes-platform"), ext="pdf")
);


#### Below code only applies for EMu k = 4

if (k == 4) {

	outliers <- as.character(filter(contribs.coord, x < 0.3, y < 0.3)$sample_id);

	filter(pheno, sample_id %in% outliers);
	filter(tcga.pheno, sample_id %in% outliers);

	contribs[outliers, ]

	outlier.mut.genes <- as.character(filter(tcga.mut, sample_id %in% outliers, Variant_Classification %in% lof.muts)$Hugo_Symbol);

	match(brca.genes, outlier.mut.genes)
	match(mmr.genes, outlier.mut.genes)
	match("TP53", outlier.mut.genes)
	match(apobec.genes, outlier.mut.genes)
	match("POLE", outlier.mut.genes)
	match(c("TP53", "KRAS", "EGFR", "KEAP1"), outlier.mut.genes)

	filter(tcga.mut, sample_id == "LUAD-TCGA-05-4396-TP-NB-SM-23JGI-SM-23JIH", Hugo_Symbol == "POLE")
	table(filter(tcga.mut, sample_id == "LUAD-TCGA-05-4396-TP-NB-SM-23JGI-SM-23JIH")$Variant_Classification)


	cmut <- as.matrix(read.table(cmut.fname, header=FALSE));
	rownames(cmut) <- samples;
	colnames(cmut) <- paste(rep(substitutions, each=16), context);

	opp <- as.matrix(read.table(opp.fname, header=FALSE));
	rownames(opp) <- samples;
	colnames(opp) <- paste(rep(substitutions, each=16), context);

	rates <- cmut / opp;

	barplot(cmut[outliers, ])

	which(contribs[, "C"] > 0.7)
	barplot(cmut["PB0262-P", ], las=2, cex.names=0.5)
	barplot(rates["PB0262-P", ], las=2, cex.names=0.5)

	which(contribs[, "B"] > 0.7)
	barplot(cmut["PB0324-P", ], las=2, cex.names=0.5)
	barplot(rates["PB0324-P", ], las=2, cex.names=0.5)
}
