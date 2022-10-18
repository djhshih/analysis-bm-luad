library(io);
library(dplyr);
library(ggplot2);
library(mmalign);
library(binom);
library(reshape2);
library(MatchIt);

options(plot.device=NULL)

#pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3.tsv");
pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
pkb.clin.all <- qread("~/exigens/brain-mets/annot/patient-info_stage2.tsv");

tcga.pheno.all <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");
tcga.clin.all <- qread("~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad_stage2.rds");

pdf.fname <- filename("pkb-tcga-compare", ext="pdf");
pkb.weights.fname <- filename("pkb-luad", tag="weights", ext="tsv");
tcga.weights.fname <- filename("tcga-luad", tag="weights", ext="tsv");

source("~/exigens/brain-mets/common/params.R");
source("~/exigens/brain-mets/common/theme.R");
cols <- col.cohorts;

####

pkb.pheno <- filter(pkb.pheno.all, primary_histotype == "Lung adenocarcinoma");
pkb.pheno.pick <- filter(pkb.pheno, pick, sample_type == "Brain metastasis");
pkb.patients <- as.character(unique(pkb.pheno$clinical_id));

pkb.clin <- filter(pkb.clin.all, clinical_id %in% pkb.patients);

# add purity of picked representative sample
pkb.clin <- left_join(pkb.clin, select(pkb.pheno.pick, sample_id, clinical_id, purity, ploidy));
filter(pkb.pheno, clinical_id == "PB0164")



tcga.pheno <- filter(tcga.pheno.all, sample_type == "Tumor");
tcga.patients <- as.character(unique(tcga.pheno$clinical_id));

tcga.clin <- filter(tcga.clin.all, clinical_id %in% tcga.patients) %>%
	mutate(
		age_at_primary_dx = age_at_initial_pathologic_diagnosis,
		smoking_hx = number_pack_years_smoked
	)

tcga.clin <- left_join(tcga.clin, select(tcga.pheno, sample_id, clinical_id, purity, ploidy));

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

summary(pkb.clin$age_at_primary_dx)
summary(tcga.clin$age_at_primary_dx)

summary(pkb.clin$sex)
summary(filter(pkb.clin, exac_pop != "EAS")$sex)
summary(tcga.clin$sex)

summary(pkb.clin$smoking_signature_high)
summary(pkb.clin$smoking_signature_high) / pkb.n
summary(tcga.clin$smoking_signature_high)
summary(tcga.clin$smoking_signature_high) / tcga.n

summary(pkb.clin$smoker)
summary(pkb.clin$smoker) / pkb.n
summary(tcga.clin$smoker)
summary(tcga.clin$smoker) / tcga.n

summary(pkb.clin$race) / sum(!is.na(pkb.clin$race))
summary(tcga.clin$race) / sum(!is.na(tcga.clin$race))

summary(pkb.clin$exac_pop)
summary(pkb.clin$exac_pop) / pkb.n
summary(tcga.clin$exac_pop)
summary(tcga.clin$exac_pop) / tcga.n

summary(filter(pkb.clin, exac_pop != "EAS")$sex) 
summary(filter(tcga.clin, exac_pop != "EAS")$sex) 


# Hung 2016
# Brain metastasis from LUAD, Chinese
# 39 male, 42 female
#
# Han 2016
# Brain metastasis from LUAD, Chinese
# 13 male, 26 female
#
# Zhao 2014
# Brain metastasis from NSCLC, Chinese
# 8 male, 9 female
#
# Li 2015
# Brain metastasis from non-squamous NSCLC, Chinese
# 57 male, 43 female
#
# Rau 2016
# Brain metastasis from LUAD, Chinese
# 27 male, 22 female
#
# Lee 2014
# Brain metastasis from LUAD, Korean
# 127 male, 165 female
#
# Riihimaki 2014
# Brain metastasis from LUAD, White (Swedish)
# 1817 male, 2001 female
#
# Renaud 2016
# Brain metastasis from LUAD, White (French)
# 27 male, 24 female
#
# Lohinai 2016
# Brain metastasis from LUAD, White (Hungarian)
# 36 male, 48 female
#
# Oh 2009
# Brain metastasis from NSCLC, White (MD Anderson)
# 378 male, 291 female
#
# Welsh 2013
# Brain metastasis from LUAD, White (MD Anderson)
# 17 male, 23 female
#
# Lin 2015
# Metastatic EGFR+ LUAD, White (DFCI)
# 31 male, 106 female
#
# Eichler 2010
# Brain metastasis from LUAD, White (MGH)
# 31 male, 62 female
#
# PKB
# Brain metastasis from LUAD, White
# 24 male, 49 female
#
# TCGA
# LUAD, White
# 239 male, 264 female

fisher.test(matrix(c(24, 239, 49, 264), nrow=2))

binom.confint(24, 24 + 49, conf.level=0.95, method="exact")
binom.confint(239, 239 + 264, conf.level=0.95, method="exact")



meta.gender <- qread("meta_gender.tsv");

confints <- do.call(rbind, apply(meta.gender[, c("female", "male")], 1,
	function(x) {
		binom.confint(x[1], sum(x), conf.level=0.90, method="exact")
	}
));

meta.gender.df <- cbind(meta.gender, confints);

meta.gender.df$label <- as.character(meta.gender.df$study);
meta.gender.df$label <- sprintf("%s (%d)", meta.gender.df$study, meta.gender.df$n);
meta.gender.df$label <- factor(meta.gender.df$label, levels=meta.gender.df$label);

qdraw(
	ggplot(meta.gender.df, aes(x=label, y=mean, ymin=lower, ymax=upper, colour=race, alpha=n)) + coord_flip() +
		facet_grid(patient ~ ., scales="free", space="free") +
		scale_alpha_continuous(trans="log10", limits=c(10, 1000)) +
		theme_bw() +
		theme(panel.grid.major.y = element_blank()) +
		theme(strip.text.y = element_text(angle = 0)) +
		geom_point() + geom_errorbar(width=0.2) + 
		geom_hline(yintercept=0.5) +
		xlab("") + ylab("Proportion of females")
	,
	width = 7,
	file = insert(pdf.fname, c("age", "meta"))
);

####

pkb.clin.sel <- select(pkb.clin, clinical_id, sample_id, age_at_primary_dx, sex, smoking_signature_high, exac_pop, stage = stage_clean) %>% mutate(cohort="BM-LUAD");
tcga.clin.sel <- select(tcga.clin, clinical_id, sample_id, age_at_primary_dx, sex, smoking_signature_high, exac_pop, stage = stage_clean) %>% mutate(cohort="TCGA-LUAD");

clin.sel <- rbind(pkb.clin.sel, tcga.clin.sel) %>%
	mutate( age_group = factor(ifelse(age_at_primary_dx < 45, "Young", ifelse(age_at_primary_dx < 65, "Old", "Older")),
			levels = c("Young", "Old", "Older")),
		age_group2 = factor(ifelse(age_at_primary_dx < 60, "Not older", "Older"))
	);

clin.sel$cohort <- factor(clin.sel$cohort, levels=cohorts);
clin.sel$cohorti  <- as.integer(clin.sel$cohort) - 1;

clin.sel.nona <- na.omit(select(clin.sel, -clinical_id, -sample_id));

exac.pop.design <- model.matrix(~ exac_pop, clin.sel.nona);

clin.mat <- select(clin.sel.nona, age=age_at_primary_dx, sex, smoking_signature_high) %>% 
	mutate(sex = as.integer(sex) - 1, smoking_signature_high = as.integer(smoking_signature_high) - 1, age=scale(age)) %>%
	cbind(exac.pop.design[, -1]) %>%
	as.matrix();

clin.mat.pca <- pca(t(clin.mat));

qdraw(
	mmalign:::pca_plot_base(jitter(clin.mat.pca$Z, factor=300), pheno=clin.sel.nona, mapping=aes(colour=cohort), vars=NULL, dims=1:2) +
		geom_point(alpha=0.6) +
		scale_colour_manual(values=cols)
	,
	height = 6, width = 7,
	file = insert(pdf.fname, c("pca", "jitter", "dims-1-2"))
);

qdraw(
	mmalign:::pca_plot_base(jitter(clin.mat.pca$Z, factor=300), pheno=clin.sel.nona, mapping=aes(colour=cohort), vars=NULL, dims=2:3) +
		geom_point(alpha=0.6) +
		scale_colour_manual(values=cols)
	,
	height = 6, width = 7,
	file = insert(pdf.fname, c("pca", "jitter", "dims-2-3"))
);

qdraw(
	mmalign:::pca_plot_base(jitter(clin.mat.pca$Z, factor=300), pheno=clin.sel.nona, mapping=aes(colour=cohort), vars=NULL, dims=3:4) +
		geom_point(alpha=0.6) +
		scale_colour_manual(values=cols)
	,
	height = 6, width = 7,
	file = insert(pdf.fname, c("pca", "jitter", "dims-3-4"))
);


pca.cor <- matrix(NA, nrow=nrow(clin.mat.pca$Z), ncol=ncol(clin.mat),
	dimnames=list(pc=rownames(clin.mat.pca$Z), covariate=colnames(clin.mat)));
for (i in 1:nrow(clin.mat.pca$Z)) {
	for (j in 1:ncol(clin.mat)) {
		pca.cor[i,j] <- cor(clin.mat.pca$Z[i,], clin.mat[, j]);
	}
}

pca.cor.m <- melt(pca.cor);

qdraw(
	ggplot(pca.cor.m, aes(x=pc, y=value, fill=covariate)) +
		geom_col() + theme_bw() +
		xlab("Principle component") + ylab("Pearson correlation")
	,
	file = insert(pdf.fname, c("pca", "cor"))
);


qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=sex, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(smoking_signature_high ~ .) + 
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "sex", "smoking"))
);

qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=sex, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(smoking_signature_high ~ age_group) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 12, height = 5,
	file = insert(pdf.fname, c("exac-pop", "sex", "smoking", "age-group"))
);

qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=sex, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(smoking_signature_high ~ age_group2) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 12, height = 5,
	file = insert(pdf.fname, c("exac-pop", "sex", "smoking", "age-group2"))
);

qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=smoking_signature_high, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(sex ~ .) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "smoking", "sex"))
);

qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=smoking_signature_high, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(sex ~ age_group) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 12, height = 5,
	file = insert(pdf.fname, c("exac-pop", "smoking", "sex", "age_group"))
);

qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=age_at_primary_dx, shape=sex, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(smoking_signature_high ~ .) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=c(45, 65), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "age", "smoking"))
);

mi0 <- matchit(cohorti ~ age_at_primary_dx + sex + smoking_signature_high + exac_pop, data=clin.sel.nona, method="cem");
mi0

mi1 <- matchit(cohorti ~ sex + exac_pop + smoking_signature_high + age_group, data=clin.sel.nona, method="cem");
mi1

unmatched.imbal <- cem::imbalance(clin.sel.nona$cohorti, clin.sel.nona, drop=c("cohort", "cohorti"));

cem::imbalance(clin.sel.nona$cohorti, clin.sel.nona, drop=c("cohort", "cohorti"), weights=mi1$weights)

filter(clin.sel.nona, mi1$weight == 0, cohort == "BM-LUAD")

mi2 <- matchit(cohorti ~ sex + exac_pop + smoking_signature_high + age_group2, data=clin.sel.nona, method="cem");
mi2

cem::imbalance(clin.sel.nona$cohorti, clin.sel.nona, drop=c("cohort", "cohorti"), weights=mi2$weights)

filter(clin.sel.nona, mi2$weight == 0, cohort == "BM-LUAD")

summary(mi2$weight[clin.sel.nona$cohort == "BM-LUAD"])

mi <- matchit(cohorti ~ sex + exac_pop + smoking_signature_high, data=select(clin.sel, sex, smoking_signature_high, exac_pop, cohorti), method="cem");

mi.imbal <- cem::imbalance(clin.sel$cohorti, select(clin.sel, -clinical_id, -sample_id), drop=c("cohort", "cohorti"), weights=mi$weights);


l1.df <- rbind(
	with(unmatched.imbal$tab, data.frame(
		matching = "Before",
		covariate = rownames(unmatched.imbal$tab),
		type = type,
		l1 = L1
	)),
	with(mi.imbal$tab, data.frame(
		matching = "After",
		covariate = rownames(mi.imbal$tab),
		type = type,
		l1 = L1
	))
);

qdraw(
	ggplot(filter(l1.df, type == "(Chi2)"), aes(x=covariate, y=l1)) +
		geom_col() + facet_grid(. ~ matching) + theme_bw() +
		coord_flip() +
		ylab("L1 imbalance between frequencies in case vs. control") +
		xlab("Covariate")
	,
	width = 6, height = 2,
	file = insert(pdf.fname, tag=c("l1-imbalance", "cem"))
);

#plot(mi);
#plot(mi, type="hist");
#plot(mi, type="jitter");

####

clin.sel.m <- mutate(clin.sel,
	weight = mi$weights
);

qdraw(
	ggplot(clin.sel.m, aes(x=exac_pop, y=smoking_signature_high, colour=cohort, alpha=log(weight+1))) +
		geom_jitter() + theme_bw() +
		facet_grid(sex ~ .) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols) +
		xlab("ExAC genetic population") + ylab("Smoking exposure") +
		theme(strip.text.y = element_text(angle = 0))
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "smoking", "sex", "cem", "alpha"))
);

qdraw(
	ggplot(clin.sel.m, aes(x=exac_pop, y=smoking_signature_high, colour=cohort, size=weight)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(sex ~ .) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols) +
		xlab("ExAC genetic population") + ylab("Smoking exposure") +
		theme(strip.text.y = element_text(angle = 0))
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "smoking", "sex", "cem", "size"))
);

qdraw(
	ggplot(clin.sel.m, aes(x=exac_pop, y=smoking_signature_high, colour=cohort, size=weight, alpha=log(weight+1))) +
		geom_jitter() + theme_bw() +
		facet_grid(sex ~ .) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		geom_hline(yintercept=seq(1.5, 2, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols) +
		xlab("ExAC genetic population") + ylab("Smoking exposure") +
		theme(strip.text.y = element_text(angle = 0))
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "smoking", "sex", "cem", "size", "alpha"))
);

####

summary(lm(age_at_primary_dx ~ cohort, clin.sel))

lapply(split(clin.sel.m, clin.sel.m$cohort), summary)

####

clin.sel.unmatched <- clin.sel;

clin.sel.m.pkb <- filter(clin.sel.m, cohort == "BM-LUAD") %>%
	mutate(weight_norm = weight / sum(weight));
clin.sel.m.tcga <- filter(clin.sel.m, cohort == "TCGA-LUAD") %>%
	mutate(weight_norm = weight / sum(weight));

set.seed(1234);
boot.idx <- sample.int(nrow(clin.sel.m.tcga), 2000, replace=TRUE, prob=clin.sel.m.tcga$weight_norm);
clin.sel.m.tcga.boot <- clin.sel.m.tcga[boot.idx, ];
summary(clin.sel.m.tcga);
summary(clin.sel.m.pkb)
summary(clin.sel.m.tcga.boot)

clin.sel.matched <- rbind(clin.sel.m.pkb, clin.sel.m.tcga.boot);

clin.sel.matching <- list(
	before = clin.sel.unmatched,
	after = clin.sel.matched
);

###

alphav <- 0.8;


for (matching in names(clin.sel.matching)) {

	clin.sel <- clin.sel.matching[[matching]];

	qdraw(
		ggplot(clin.sel, aes(x=cohort, y=age_at_primary_dx, fill=cohort)) +
			geom_violin(alpha=alphav) + theme_bw() +
			scale_fill_manual(values=cols) +
			theme(legend.position="none") +
			xlab("") + ylab("Age at primary diagnosis")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "age", "violin"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, y=age_at_primary_dx, fill=cohort)) +
			geom_boxplot(alpha=alphav, notch=TRUE) + theme_bw() +
			scale_fill_manual(values=cols) +
			theme(legend.position="none") +
			xlab("") + ylab("Age at primary diagnosis")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "age", "box"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, y=age_at_primary_dx, colour=cohort)) +
			geom_jitter(width=0.2, alpha=alphav) + theme_bw() +
			scale_colour_manual(values=cols) +
			theme(legend.position="none") +
			xlab("") + ylab("Age at primary diagnosis")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "age", "jitter"))
	);


	wilcox.test(age_at_primary_dx ~ cohort, data=clin.sel)

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=sex)) +
			geom_bar() + theme_bw() +
			xlab("") + ylab("Count")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "sex", "bar"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=sex)) +
			geom_bar(position="fill") + theme_bw() +
			xlab("") + ylab("Proportion")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "sex", "bar", "fill"))
	);

	sex.df <- do.call(rbind, lapply(
		split(clin.sel, clin.sel$cohort),
		function(d) {
			z <- with(d, table(sex));
			binom.confint(z[1], sum(z), conf.level=0.9, method="ac")
		}
	));
	sex.df$cohort <- factor(rownames(sex.df), levels=rownames(sex.df));

	qdraw(
		ggplot(sex.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) +
			geom_col(alpha=alphav) + geom_errorbar(width=0.2) + theme_bar() +
			xlab("") + ylab("Proportion of females") +
			scale_fill_manual(values=cols) + ylim(0, 1)
		,
		width = 1.5, height = 3,
		file = insert(pdf.fname, c(matching, "sex", "female-prop"))
	);

	qdraw(
		ggplot(sex.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, colour=cohort)) +
			geom_point() + 
			geom_errorbar(width=0.1) + theme_vdot() +
			xlab("") + ylab("Proportion of females") +
			scale_colour_manual(values=cols) + ylim(0, 1)
		,
		width = 1.2, height = 3,
		file = insert(pdf.fname, c(matching, "sex", "female-prop", "vdot"))
	);

	sex.ct <- with(clin.sel, table(cohort, sex));
	fisher.test(sex.ct)

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=exac_pop)) +
			geom_bar() + theme_bw() +
			xlab("") + ylab("Count")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "exac-pop", "bar"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=exac_pop)) +
			geom_bar(position="fill") + theme_bw() +
			xlab("") + ylab("Proportion")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "exac-pop", "bar", "fill"))
	);

	smoking.df <- do.call(rbind, mapply(
		function(cohort, d) {
			z <- with(d, table(smoking_signature_high));
			data.frame(
				cohort = cohort,
				smoking_signature_high = names(z)[1],
				binom.confint(z[1], sum(z), conf.level=0.9, method="ac")
			)
		},
		levels(clin.sel$cohort),
		split(clin.sel, clin.sel$cohort),
		SIMPLIFY = FALSE
	));
	smoking.df$lower <- pmax(smoking.df$lower, 0);
	rownames(smoking.df) <- NULL;

	qdraw(
		ggplot(smoking.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) +
			geom_col(alpha=alphav) + geom_errorbar(width=0.2) + theme_bar() +
			xlab("") + ylab("Proportion of individuals with high smoking exposure") +
			scale_fill_manual(values=cols) + ylim(0, 1)
		,
		width = 1.5, height = 3,
		file = insert(pdf.fname, c(matching, "smoking", "prop"))
	);

	qdraw(
		ggplot(smoking.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, colour=cohort)) +
			geom_point() + geom_errorbar(width=0.1) + theme_vdot() +
			xlab("") + ylab("Proportion of individuals with high smoking exposure") +
			scale_colour_manual(values=cols) + ylim(0, 1)
		,
		width = 1.2, height = 3,
		file = insert(pdf.fname, c(matching, "smoking", "prop", "vdot"))
	);

	exac.df <- do.call(rbind, mapply(
		function(cohort, d) {
			z <- with(d, table(exac_pop));
			data.frame(
				cohort = cohort,
				do.call(rbind, lapply(1:length(z), function(i) {
					data.frame(
						exac_pop = names(z)[i],
						binom.confint(z[i], sum(z), conf.level=0.9, method="ac")
					)
				}))
			)
		},
		levels(clin.sel$cohort),
		split(clin.sel, clin.sel$cohort),
		SIMPLIFY = FALSE
	));
	exac.df$lower <- pmax(exac.df$lower, 0);
	rownames(exac.df) <- NULL;

	qdraw(
		ggplot(exac.df, aes(x=exac_pop, y=mean, ymin=lower, ymax=upper, fill=cohort)) +
			geom_bar(stat="identity", position=position_dodge(), alpha=alphav) + 
			geom_errorbar(width=0.2, position=position_dodge(0.9)) + theme_bw() +
			xlab("") + ylab("Proportion") +
			theme(panel.grid = element_blank(), legend.position = "bottom") +
			scale_fill_manual(values=cols)
		,
		width = 6, height = 3,
		file = insert(pdf.fname, c(matching, "exac-pop", "prop"))
	);

	qdraw(
		ggplot(exac.df, aes(x=exac_pop, y=mean, ymin=lower, ymax=upper, colour=cohort)) +
			geom_point(position=position_dodge(0.6)) + 
			geom_errorbar(width=0.1, position=position_dodge(0.6)) + theme_bw() +
			xlab("") + ylab("Proportion") +
			theme(panel.grid = element_blank(), legend.position = "bottom") +
			scale_colour_manual(values=cols) + ylim(0, NA)
		,
		width = 4, height = 3,
		file = insert(pdf.fname, c(matching, "exac-pop", "prop", "vdot"))
	);

	exac.ct <- with(clin.sel, table(cohort, exac_pop));
	exac.ct
	fisher.test(exac.ct)

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=stage)) +
			geom_bar() + theme_bw() +
			xlab("") + ylab("Count")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "stage", "bar"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=stage)) +
			geom_bar(position="fill") + theme_bw() +
			xlab("") + ylab("Proportion")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "stage", "bar", "fill"))
	);

	stage.wide.df <- do.call(rbind, mapply(
		function(cohort, d) {
			z <- with(d, table(stage));
			data.frame(
				cohort = cohort,
				do.call(rbind, lapply(1:length(z), function(i) {
					data.frame(
						stage = names(z)[i],
						binom.confint(z[i], sum(z), conf.level=0.9, method="ac")
					)
				}))
			)
		},
		levels(clin.sel$cohort),
		split(clin.sel, clin.sel$cohort),
		SIMPLIFY = FALSE
	));
	stage.wide.df$lower <- pmax(stage.wide.df$lower, 0);
	rownames(stage.wide.df) <- NULL;

	qdraw(
		ggplot(stage.wide.df, aes(x=stage, y=mean, ymin=lower, ymax=upper, fill=cohort)) +
			geom_bar(stat="identity", position=position_dodge(), alpha=alphav) + 
			geom_errorbar(width=0.2, position=position_dodge(0.9)) + theme_bw() +
			xlab("Stage") + ylab("Proportion") +
			theme(panel.grid = element_blank(), legend.position = "bottom") +
			scale_fill_manual(values=cols)
		,
		width = 4, height = 3,
		file = insert(pdf.fname, c(matching, "stage", "prop"))
	);

	qdraw(
		ggplot(stage.wide.df, aes(x=stage, y=mean, ymin=lower, ymax=upper, colour=cohort)) +
			geom_point(position=position_dodge(0.6)) + 
			geom_errorbar(width=0.1, position=position_dodge(0.6)) + theme_bw() +
			xlab("Stage") + ylab("Proportion") +
			theme(panel.grid = element_blank(), legend.position = "bottom") +
			scale_colour_manual(values=cols) + ylim(0, NA)
		,
		width = 3, height = 3,
		file = insert(pdf.fname, c(matching, "stage", "prop", "vdot"))
	);

	#stage.ct <- with(stage.df, table(cohort, stage));
	#stage.ct
	#fisher.test(stage.ct)
}

###

# FIXME This code is fragile

#stopifnot(all(clin.sel.m.tcga$age_at_primary_dx == tcga.clin$age_at_primary_dx, na.rm=TRUE))
#stopifnot(all(clin.sel.m.tcga$sex == tcga.clin$sex, na.rm=TRUE))
#stopifnot(all(as.character(clin.sel.m.tcga$exac_pop) == as.character(tcga.clin$exac_pop), na.rm=TRUE))
#
#tcga.weights <- data.frame(
#	clinical_id = tcga.clin$clinical_id,
#	weight = clin.sel.m.tcga$weight,
#	weight_norm = clin.sel.m.tcga$weight_norm
#);
#
#stopifnot(all(clin.sel.m.pkb$age_at_primary_dx == pkb.clin$age_at_primary_dx, na.rm=TRUE))
#stopifnot(all(clin.sel.m.pkb$sex == pkb.clin$sex, na.rm=TRUE))
#stopifnot(all(as.character(clin.sel.m.pkb$exac_pop) == as.character(pkb.clin$exac_pop), na.rm=TRUE))
#
#pkb.weights <- data.frame(
#	clinical_id = pkb.clin$clinical_id,
#	weight = clin.sel.m.pkb$weight,
#	weight_norm = clin.sel.m.pkb$weight_norm
#);


tcga.weights <- left_join(select(tcga.clin, clinical_id), select(clin.sel.m.tcga, clinical_id, weight, weight_norm)) %>%
	mutate(weight_norm = ifelse(is.na(weight_norm), 0, weight_norm));

pkb.weights <- left_join(select(pkb.clin, clinical_id), select(clin.sel.m.pkb, clinical_id, weight, weight_norm)) %>%
	mutate(weight_norm = ifelse(is.na(weight_norm), 0, weight_norm));

####

qwrite(tcga.weights, tcga.weights.fname);
qwrite(pkb.weights, pkb.weights.fname);

####

qdraw(
	{
		plot(ecdf(tcga.clin$smoking_hx[tcga.clin$smoking_hx > 0]), col=col.control,
			las = 1,
			main = "Smoking history among ever smokers",
			xlab = "Smoking pack-year",
			ylab = "Cumulative probability"
		);
		lines(ecdf(pkb.clin$smoking_hx[pkb.clin$smoking_hx > 0]), col=col.case);
		legend("bottomright", legend=cohorts, pch=19, col=cols, bty="n", inset=0.05);
	},
	file = insert(pdf.fname, tag=c("ecdf", "smoking"))
);

qdraw(
	{
		plot(ecdf(tcga.clin$age_at_primary_dx), col=col.control,
			las = 1,
			main = "Age at primary diagnosis",
			xlab = "Age (year)",
			ylab = "Cumulative probability"
		);
		lines(ecdf(pkb.clin$age_at_primary_dx), col=col.case);
		legend("bottomright", legend=cohorts, pch=19, col=cols, bty="n", inset=0.05);
	},
	file = insert(pdf.fname, tag=c("ecdf", "age"))
);

####

with(tcga.clin, table(smoker, smoking_signature_high))
with(pkb.clin, table(smoker, smoking_signature_high))

smoking.df <- rbind(
	transmute(tcga.clin, cohort="TCGA-LUAD", clinical_id, smoker, smoking_signature_high, smoking_hx),
	transmute(pkb.clin, cohort="BM-LUAD", clinical_id, smoker, smoking_signature_high, smoking_hx)
);

qdraw(
	ggplot(smoking.df, aes(x=smoker, fill=smoking_signature_high)) +
		geom_bar(position="fill") + theme_bw() +
		xlab("Smoker") + ylab("Proportion") + labs(fill="Smoking exposure")
	,
	width = 3.5,
	file = insert(pdf.fname, c("smoker", "smoking", "combined"))
);

qdraw(
	ggplot(smoking.df, aes(x=smoker, fill=smoking_signature_high)) +
		geom_bar(position="fill") + theme_bw() +
		facet_grid(. ~ cohort) +
		xlab("Smoker") + ylab("Proportion") + labs(fill="Smoking exposure")
	,
	width = 6,
	file = insert(pdf.fname, c("smoker", "smoking"))
);

