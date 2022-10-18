library(io);
library(dplyr);
library(ggplot2);
library(mmalign);
library(binom);
library(reshape2);
library(MatchIt);

#pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3.tsv");
pkb.pheno.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
pkb.clin.all <- qread("~/exigens/brain-mets/annot/patient-info_stage2.tsv");

tcga.pheno.all <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");
tcga.clin.all <- qread("~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad_stage2.rds");

pdf.fname <- filename("pkb-tcga-compare", ext="pdf");
pkb.weights.fname <- filename("pkb-luad", tag="weights", ext="tsv");
tcga.weights.fname <- filename("tcga-luad", tag="weights", ext="tsv");

source("~/exigens/brain-mets/common/params.R");
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

tmp <- rbind(
	select(pkb.pheno.pick, clinical_id, sample_id, purity, ploidy) %>% mutate(cohort = "BM-LUAD"),
	select(tcga.pheno, clinical_id, sample_id, purity, ploidy) %>% mutate(cohort = "TCGA-LUAD")
);

ggplot(tmp, aes(x = cohort, y = purity)) +
	geom_jitter()

ggplot(tmp, aes(x = cohort, y = purity)) +
	geom_boxplot()

ggplot(tmp, aes(x = cohort, y = purity)) +
	geom_violin()

wilcox.test(purity ~ cohort, tmp)
t.test(purity ~ cohort, tmp)

# purity is on the backdoor path and it is not a collider
# we should control for purity


ggplot(tmp, aes(x = cohort, y = ploidy)) +
	geom_jitter()

ggplot(tmp, aes(x = cohort, y = ploidy)) +
	geom_boxplot()

ggplot(tmp, aes(x = cohort, y = ploidy)) +
	geom_violin()

wilcox.test(ploidy ~ cohort, tmp)
t.test(ploidy ~ cohort, tmp)

# therefore, no need to corect to ploidy


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

pkb.clin.sel <- select(pkb.clin, clinical_id, sample_id, age_at_primary_dx, sex, smoking_signature_high, exac_pop, purity, stage = stage_clean) %>% mutate(cohort="BM-LUAD");
tcga.clin.sel <- select(tcga.clin, clinical_id, sample_id, age_at_primary_dx, sex, smoking_signature_high, exac_pop, purity, stage = stage_clean) %>% mutate(cohort="TCGA-LUAD") %>%
	filter(stage == "III");

clin.sel <- rbind(pkb.clin.sel, tcga.clin.sel) %>%
	mutate( age_group = factor(ifelse(age_at_primary_dx < 45, "Young", ifelse(age_at_primary_dx < 65, "Old", "Older")),
			levels = c("Young", "Old", "Older")),
		age_group2 = factor(ifelse(age_at_primary_dx < 60, "Not older", "Older"))
	);

clin.sel$cohort <- factor(clin.sel$cohort, levels=cohorts);
clin.sel$cohorti  <- as.integer(clin.sel$cohort) - 1;

qdraw(
	ggplot(clin.sel, aes(x=exac_pop, y=purity, shape=sex, colour=cohort)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(smoking_signature_high ~ .) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols)
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "purity", "smoking"))
);

nona.idx <- complete.cases(select(clin.sel, sex, smoking_signature_high, exac_pop, cohorti, purity));
clin.sel.nona <- clin.sel[nona.idx, ];

mi <- matchit(cohorti ~ sex + exac_pop + smoking_signature_high + purity,
	data = select(clin.sel.nona, sex, smoking_signature_high, exac_pop, cohorti, purity),
	method="cem");
print(mi)
# only 44 cases are matched
# one case sample (PB0164-M) had no purity estimate due to complex copy-number that precluded
# ABSOLUTE analysis

# Therefore, we cannot adjust for purity due to lack of power.

nona.idx <- complete.cases(select(clin.sel, sex, smoking_signature_high, exac_pop, cohorti));
clin.sel.nona <- clin.sel[nona.idx, ];

mi <- matchit(cohorti ~ sex + exac_pop + smoking_signature_high,
	data = select(clin.sel.nona, sex, smoking_signature_high, exac_pop, cohorti),
	method="cem");
print(mi)


mi.imbal <- cem::imbalance(clin.sel$cohorti, select(clin.sel, -clinical_id, -sample_id), drop=c("cohort", "cohorti"), weights=mi$weights);
print(mi.imbal)

unmatched.imbal <- cem::imbalance(clin.sel$cohorti, select(clin.sel, -clinical_id, -sample_id), drop=c("cohort", "cohorti"));

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

clin.sel.m <- mutate(clin.sel.nona,
	weight = mi$weights
);

qdraw(
	ggplot(clin.sel.m, aes(x=exac_pop, y=purity, colour=cohort, alpha=log(weight+1))) +
		geom_jitter() + theme_bw() +
		facet_grid(sex ~ smoking_signature_high) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols) +
		xlab("ExAC genetic population") + ylab("Tumor purity") +
		theme(strip.text.y = element_text(angle = 0))
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "purity", "sex", "smoking", "cem", "alpha"))
);

qdraw(
	ggplot(clin.sel.m, aes(x=exac_pop, y=purity, colour=cohort, size=weight)) +
		geom_jitter(alpha=0.6) + theme_bw() +
		facet_grid(sex ~ smoking_signature_high) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols) +
		xlab("ExAC genetic population") + ylab("Tumor purity") +
		theme(strip.text.y = element_text(angle = 0))
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "purity", "sex", "smoking", "cem", "size"))
);

qdraw(
	ggplot(clin.sel.m, aes(x=exac_pop, y=purity, colour=cohort, size=weight, alpha=log(weight+1))) +
		geom_jitter() + theme_bw() +
		facet_grid(sex ~ smoking_signature_high) +
		geom_vline(xintercept=seq(1.5, 6, 1), colour="grey60") +
		theme(panel.grid.major = element_blank()) +
		scale_colour_manual(values=cols) +
		xlab("ExAC genetic population") + ylab("Tumor purity") +
		theme(strip.text.y = element_text(angle = 0))
	,
	width = 8, height = 5,
	file = insert(pdf.fname, c("exac-pop", "purity", "sex", "smoking", "cem", "size", "alpha"))
);

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

# assume that all non-zero case samples have the same weight
clin.sel.matched <- rbind(clin.sel.m.pkb[clin.sel.m.pkb$weight > 0, ], clin.sel.m.tcga.boot);

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
		ggplot(clin.sel, aes(x=cohort, y=purity, fill=cohort)) +
			geom_violin(alpha=alphav) + theme_bw() +
			scale_fill_manual(values=cols) +
			theme(legend.position="none") +
			xlab("") + ylab("Tumor purity")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "purity", "violin"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, y=purity, fill=cohort)) +
			geom_boxplot(alpha=alphav, notch=TRUE) + theme_bw() +
			scale_fill_manual(values=cols) +
			theme(legend.position="none") +
			xlab("") + ylab("Tumor purity")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "purity", "box"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, y=purity, colour=cohort)) +
			geom_jitter(width=0.2, alpha=alphav) + theme_bw() +
			scale_colour_manual(values=cols) +
			theme(legend.position="none") +
			xlab("") + ylab("Tumor purity")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "purity", "jitter"))
	);

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
			geom_col(alpha=alphav) + geom_errorbar(width=0.2) + theme_bw() +
			xlab("") + ylab("Proportion of females") +
			scale_fill_manual(values=cols) +
			theme(legend.position="none")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "sex", "female-prop"))
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
			geom_col(alpha=alphav) + geom_errorbar(width=0.2) + theme_bw() +
			xlab("") + ylab("Proportion of individuals with high smoking exposure") +
			scale_fill_manual(values=cols) +
			theme(legend.position="none")
		,
		width = 2,
		file = insert(pdf.fname, c(matching, "smoking", "prop"))
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
			scale_fill_manual(values=cols)
		,
		width = 6, height = 3,
		file = insert(pdf.fname, c(matching, "exac-pop", "prop"))
	);

	exac.ct <- with(clin.sel, table(cohort, exac_pop));
	exac.ct
	fisher.test(exac.ct)

	stage.df <- rbind(
		transmute(pkb.clin, cohort="BM-LUAD", stage=stage_clean),
		transmute(tcga.clin, cohort="TCGA-LUAD", stage=stage_clean)
	);

	stage.df$cohort <- factor(stage.df$cohort, levels=cohorts);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=exac_pop)) +
			geom_bar() + theme_bw() +
			xlab("") + ylab("Count")
		,
		width = 3,
		file = insert(pdf.fname, c(matching, "stage", "bar"))
	);

	qdraw(
		ggplot(clin.sel, aes(x=cohort, fill=exac_pop)) +
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
			scale_fill_manual(values=cols)
		,
		width = 5, height = 3,
		file = insert(pdf.fname, c(matching, "stage", "prop"))
	);

	stage.ct <- with(stage.df, table(cohort, stage));
	stage.ct
	fisher.test(stage.ct)
}

###

tcga.weights <- left_join(select(tcga.clin, clinical_id), select(clin.sel.m.tcga, clinical_id, weight, weight_norm)) %>%
	mutate(weight_norm = ifelse(is.na(weight_norm), 0, weight_norm));

pkb.weights <- left_join(select(pkb.clin, clinical_id), select(clin.sel.m.pkb, clinical_id, weight, weight_norm)) %>%
	mutate(weight_norm = ifelse(is.na(weight_norm), 0, weight_norm));

qwrite(tcga.weights, tcga.weights.fname);
qwrite(pkb.weights, pkb.weights.fname);

