library(io);
library(dplyr);
library(ggplot2);

source("~/exigens/brain-mets/common/params.R");
source("~/exigens/brain-mets/common/theme.R");

out.fname <- filename("pkb-tcga-luad-compare");
pdf.fname <- insert(out.fname, ext="pdf");

pkb <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
pkb.weights <- qread("~/exigens/brain-mets/matchit/pkb-luad_weights.tsv");
tcga <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");
tcga.weights <- qread("~/exigens/brain-mets/matchit/tcga-luad_weights.tsv");

pkb.doubled <- pkb$genome_doubling > 0;
tcga.doubled <- tcga$genome_doubling > 0;

pkb.doubled.counts <- table(pkb.doubled);
tcga.doubled.counts <- table(tcga.doubled);

pkb.doubled.counts / sum(pkb.doubled.counts)
tcga.doubled.counts / sum(tcga.doubled.counts)

ct <- matrix(
	c(table(tcga.doubled), table(pkb.doubled)),
	nrow = 2,
	byrow = TRUE
);

fisher.test(ct)

pkb <- left_join(pkb, pkb.weights, by="clinical_id");
tcga <- left_join(tcga, tcga.weights, by="clinical_id");

pheno <- full_join(data.frame(tcga, cohort="TCGA-LUAD"), data.frame(pkb, cohort="BM-LUAD"));

fit <- glm(genome_doublings > 0 ~ cohort, data = pheno, family=quasibinomial, weights=weight);
summary(fit)

with(pkb, sum(weight[genome_doublings > 0], na.rm=TRUE) / sum(weight[!is.na(genome_doublings)]))
with(tcga, sum(weight[genome_doublings > 0], na.rm=TRUE) / sum(weight[!is.na(genome_doublings)]))


alpha.level <- 0.2;
alphav <- 0.8;
cols.ch <- col.cohorts;
names(cols.ch) <- cohorts;

logistic <- function(x) {
	1 / (1 + exp(-x))
}

# get t confidence intervals
fit0 <- glm(genome_doublings > 0 ~ cohort - 1, data = pheno, family=quasibinomial, weights=weight);
summary(fit0)
coefs0 <- summary(fit0)$coefficients;
lest <- coefs0[, 1];
lse <- coefs0[, 2];


# two coefficients are
a <- qt(1 - alpha.level/2, df = fit0$df.residual);
coef.df <- data.frame(
	cohort = factor(gsub("cohort", "", rownames(coefs0)), levels=cohorts),
	mean = logistic(lest),
	lower = logistic(lest - a * lse),
	upper = logistic(lest + a * lse)
);


qdraw(
	ggplot(coef.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) + 
		geom_col(alpha=alphav) +
		geom_errorbar(width=0.2) + theme_bar() + 
		scale_fill_manual(values=cols.ch) + 
		ylab("Genome doubling frequency") + xlab("") +
		ylim(0, 1)
	,
	width = 1.5, height = 3,
	file = insert(pdf.fname, c("wgd-freq", "bar"))
);

qdraw(
	ggplot(coef.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, colour=cohort)) + 
		geom_point() +
		geom_errorbar(width=0.1) + theme_vdot() + 
		scale_colour_manual(values=cols.ch) + 
		ylab("Genome doubling frequency") + xlab("") +
		ylim(0, 1)
	,
	width = 1.2, height = 3,
	file = insert(pdf.fname, c("wgd-freq", "vdot"))
);

qdraw(
	ggplot(coef.df, aes(x=cohort, y=mean, ymin=lower, ymax=upper, colour=cohort)) + 
		geom_point() +
		geom_errorbar(width=0.1) + theme_dot() + coord_flip() +
		scale_colour_manual(values=cols.ch) + 
		expand_limits(y=0) +
		ylab("Genome doubling frequency") + xlab("") +
		ylim(0, 1)
	,
	width = 3.5, height = 1,
	file = insert(pdf.fname, c("wgd-freq", "dot"))
);

print(summary(fit))

