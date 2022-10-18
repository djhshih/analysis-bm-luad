library(io);
library(dplyr);
library(ggplot2);
library(reshape2);
library(tidyr);

pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
pkb.clin.all <- qread("~/exigens/brain-mets/annot/patient-info_stage2.tsv");
pkb.weights <- qread("~/exigens/brain-mets/matchit/pkb-luad_weights.tsv");

tcga.pheno.all <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");
tcga.clin.all <- qread("~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad_stage2.rds");
tcga.weights <- qread("~/exigens/brain-mets/matchit/tcga-luad_weights.tsv");

pkb.seg <- qread("~/exigens/brain-mets/gistic/brain-mets_pass_luad_absolute.lgr.cntf.cenf.cnvf.seg");
tcga.seg <- qread("~/exigens/brain-mets/gistic/tcga-luad_pass_absolute.lgr.cntf.cenf.cnvf.seg");

pdf.fname <- filename("pkb-luad", ext="pdf");

source("~/exigens/brain-mets/common/params.R");
source("~/exigens/brain-mets/common/theme.R");
cols <- col.cohorts;

gain.cut <- log2(3/2);
loss.cut <- log2(1/2);

chrom.gain.cut <- log2(5/4);
chrom.loss.cut <- log2(3/4);

#amp.cut <- log2(8/2);
#del.cut <- log2(0.5/2);
amp.cut <- 8;
del.cut <- 0.5;

####

# filter segments

filter_segments <- function(seg) {
	lens <- with(seg, end - start + 1);
	seg[lens > 1e5, ]	
}

pkb.seg.f <- filter_segments(pkb.seg);
nrow(pkb.seg)
nrow(pkb.seg.f)

tcga.seg.f <- filter_segments(tcga.seg);
nrow(tcga.seg)
nrow(tcga.seg.f)

pkb.seg <- pkb.seg.f;
tcga.seg <- tcga.seg.f;

chromosome_segments <- function(seg) {
	group_by(seg, sample, chromosome) %>%
		mutate(length = end - start + 1) %>%
		summarize(state = sum(length * state) / sum(length));
}

pkb.chrom.seg <- chromosome_segments(pkb.seg);
tcga.chrom.seg <- chromosome_segments(tcga.seg);

####

pkb.pheno <- filter(pkb.pheno.all, primary_histotype == "Lung adenocarcinoma") %>%
	left_join(pkb.weights);
pkb.patients <- as.character(unique(pkb.pheno$clinical_id));

#pkb.pheno.pick <- filter(pkb.pheno, pick, sample_type == "Brain metastasis");

pkb.clin <- filter(pkb.clin.all, clinical_id %in% pkb.patients);

# add purity
pkb.clin <- left_join(pkb.clin, select(pkb.pheno, sample_id, clinical_id, purity, ploidy));

tcga.pheno <- filter(tcga.pheno.all, sample_type == "Tumor") %>%
	left_join(tcga.weights);
tcga.patients <- as.character(unique(tcga.pheno$clinical_id));

tcga.clin <- filter(tcga.clin.all, clinical_id %in% tcga.patients) %>%
	mutate(
		age_at_primary_dx = age_at_initial_pathologic_diagnosis,
		smoking_hx = number_pack_years_smoked
	)

tcga.clin <- left_join(tcga.clin, select(tcga.pheno, sample_id, clinical_id, purity, ploidy));

pkb.seg <- left_join(pkb.seg, select(pkb.pheno, sample=sample_id, ploidy));
pkb.seg <- mutate(pkb.seg, tcn = (2^state) * ploidy);

####

setdiff(tcga.patients, tcga.clin.all$clinical_id)

pkb.n <- nrow(pkb.clin);
tcga.n <- nrow(tcga.clin);

pkb.clin$stage[pkb.clin$stage == "IIA/IIB"] <- "IIB";

pkb.clin <- mutate(pkb.clin,
	stage_clean = factor(gsub("A|B", "", stage), levels=c("I", "II", "III", "IV"))
);

tcga.clin$stage_clean <- toupper(gsub("a|b", "", gsub("stage ", "", tcga.clin$pathologic_stage)));

####

test_to_string <- function(h, show.statistic=FALSE) {
	if (show.statistic) {
		s0 <- sprintf("%s = %.2f", names(h$statistic), h$statistic);
	} else {
		s0 <- NULL;
	}
	s1 <- ifelse(h$p.value < 0.01,
		"p < 0.01",
		sprintf("p = %.2f", h$p.value)
	);
	if (is.null(s0)) {
		s <- s1;
	} else {
		s <- paste(s0, s1, sep=",  ");
	}
	s
}

p_to_string <- function(p) {
	ifelse(p < 0.001,
		"p < 0.001",
		ifelse(p < 0.01,
			sprintf("p = %.3f", p),
			sprintf("p = %.2f", p)
		)
	)
}


####

pkb.cna <- mutate(pkb.seg,
		amp = tcn > amp.cut,
		del = tcn < del.cut,
		hcna = amp | del,
		gain = state > gain.cut,
		loss = state < loss.cut,
		cna = gain | loss
	) %>%
	group_by(sample) %>%
	summarize(
		count_hcna = sum(hcna),
		count_hamp = sum(amp),
		count_hdel = sum(del),
		count_cna = sum(cna),
		count_gain = sum(gain),
		count_loss = sum(loss)
	);

variables <- c("count_hcna", "count_hamp", "count_hdel", "count_cna", "count_gain", "count_loss");
#variables.broad <- c("count_cna", "count_gain", "count_loss");

cna.labeller <- function(value) {
	x <- list(
		count_cna = "Gains or losses",
		count_hcna = "High-level amplifications or deep deletions",
		count_hamp = "High level amplifications",
		count_hdel = "Deep deletions",
		count_gain = "Gains",
		count_loss = "Losses"
	)
	x[value]
}

pkb <- left_join(pkb.pheno, pkb.cna, by = c("sample_id"="sample")) %>%
	filter(pick, sample_type != "Normal");

pkb.sel <- select(pkb, sample_id, purity, specimen_material, count_hcna, count_hamp, count_hdel, count_cna, count_gain, count_loss);

pkb.tall <- pkb.sel %>% gather(key="variable", value="value", -sample_id, -purity, -specimen_material)

pkb.tall$variable <- factor(pkb.tall$variable, levels=variables);

pkb.tests <- lapply(variables, function(v) {
	wilcox.test(pkb[[v]] ~ pkb$specimen_material)
});

pkb.text <- data.frame(
	variable = variables,
	test = unlist(lapply(pkb.tests, test_to_string, show.statistic=FALSE))
);

qdraw(
	ggplot(pkb.tall, aes(x=specimen_material, y=value)) +
		geom_boxplot(outlier.shape=NA) + theme_box() +
		geom_jitter(aes(alpha=purity)) +
		geom_text(data = pkb.text, aes(x=1, y=Inf, label=test), hjust=0.3, vjust=6) +
		facet_wrap(~ variable, nrow=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = cna.labeller)
		) +
		ylab("Count of copy-number events") + xlab("")
	,
	width = 8, height = 4,
	file = insert(pdf.fname, c("cn-ffpe", "box"))
);

qdraw(
	ggplot(pkb.tall, aes(x=specimen_material, y=value)) +
		geom_violin(outlier.shape=NA) + theme_box() +
		geom_jitter(aes(alpha=purity)) +
		geom_text(data = pkb.text, aes(x=1, y=Inf, label=test), hjust=0.3, vjust=6) +
		facet_wrap(~ variable, nrow=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = cna.labeller)
		) +
		ylab("Count of copy-number events") + xlab("")
	,
	width = 8, height = 4,
	file = insert(pdf.fname, c("cn-ffpe", "violin"))
);


####

variables <- c("mutation_burden", "mutation_burden_indel");

mut.labeller <- function(value) {
	x <- list(
		mutation_burden = "SNVs",
		mutation_burden_indel = "indels"
	)
	x[value]
}

pkb.sel <- select(pkb, sample_id, purity, specimen_material, mutation_burden, mutation_burden_indel);

pkb.tall <- pkb.sel %>% gather(key="variable", value="value", -sample_id, -purity, -specimen_material)

pkb.tall$variable <- factor(pkb.tall$variable, levels=variables);

pkb.tests <- lapply(variables, function(v) {
	wilcox.test(pkb[[v]] ~ pkb$specimen_material)
});

pkb.text <- data.frame(
	variable = variables,
	test = unlist(lapply(pkb.tests, test_to_string, show.statistic=FALSE))
);

qdraw(
	ggplot(pkb.tall, aes(x=specimen_material, y=value)) +
		geom_boxplot(outlier.shape=NA) + theme_box() +
		geom_jitter(aes(alpha=purity)) +
		# having no vertical jiter allow too many points to overlap
		#geom_jitter(aes(alpha=purity), height0) +
		geom_text(data = pkb.text, aes(x=1, y=Inf, label=test), hjust=0.3, vjust=4) +
		facet_wrap(~ variable, nrow=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = mut.labeller)
		) +
		scale_y_log10() +
		ylab("Mutations per Mbp") + xlab("")
	,
	width = 3.75, height = 4,
	file = insert(pdf.fname, c("mut-ffpe", "box"))
);

