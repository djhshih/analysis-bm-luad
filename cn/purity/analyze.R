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

pkb.seg <- qread("~/exigens/brain-mets/gistic/brain-mets_pass_luad_absolute_bmet-only.lgr.cntf.cenf.cnvf.seg");
tcga.seg <- qread("~/exigens/brain-mets/gistic/tcga-luad_pass_absolute.lgr.cntf.cenf.cnvf.seg");

pdf.fname <- filename("cn-purity", ext="pdf");

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

hist(pkb.seg$state, breaks=100)
hist(tcga.seg$state, breaks=100)

hist(tcga.chrom.seg$state, breaks=100)

mean(tcga.chrom.seg$state < log2(1/2))
mean(tcga.chrom.seg$state > log2(3/2))

####

pkb.pheno <- filter(pkb.pheno.all, primary_histotype == "Lung adenocarcinoma") %>%
	left_join(pkb.weights);
pkb.patients <- as.character(unique(pkb.pheno$clinical_id));

#pkb.pheno.pick <- filter(pkb.pheno, sample_type == "Brain metastasis");
pkb.pheno.pick <- filter(pkb.pheno, pick, sample_type == "Brain metastasis");

pkb.clin <- filter(pkb.clin.all, clinical_id %in% pkb.patients);

# add purity of picked representative sample
pkb.clin <- left_join(pkb.clin, select(pkb.pheno.pick, sample_id, clinical_id, purity, ploidy));

#filter(pkb.pheno, clinical_id == "PB0164")

tcga.pheno <- filter(tcga.pheno.all, sample_type == "Tumor") %>%
	left_join(tcga.weights);
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


pkb.seg <- left_join(pkb.seg, select(pkb.pheno, sample=sample_id, ploidy));
pkb.seg <- mutate(pkb.seg, tcn = (2^state) * ploidy);

tcga.seg <- left_join(tcga.seg, select(tcga.pheno, sample=sample_id, ploidy));
tcga.seg <- mutate(tcga.seg, tcn = (2^state) * ploidy);

####

cortest_to_string <- function(h, statistic="Kendall tau") {
#cortest_to_string <- function(h, statistic="Pearson r") {
	paste(
		sprintf("%s = %.2f", statistic, h$estimate),
		ifelse(h$p.value < 0.01,
			"p < 0.01",
			sprintf("p = %.2f", h$p.value)
		)
		,
		sep=",  "
	)
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

variables <- c("count_hcna", "count_hamp", "count_hdel", "count_cna", "count_gain", "count_loss");
variables.broad <- c("count_cna", "count_gain", "count_loss");

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

outliers <- filter(pkb.cna, count_cna > 500)$sample;
# visual inspection of copy-number profiles revealed that the outliers
# PB0103-M and PB0251-M have noisy profiles

pkb.cna[pkb.cna$sample %in% outliers, ] <- NA;

pkb <- left_join(pkb.pheno, pkb.cna, by = c("sample_id"="sample"));

pkb.sel <- select(pkb, sample_id, purity, count_hcna, count_hamp, count_hdel, count_cna, count_gain, count_loss, weight, weight_norm);

pkb.tall <- pkb.sel %>% gather(key="variable", value="value", -sample_id, -purity, -weight, -weight_norm)

pkb.tall$variable <- factor(pkb.tall$variable, levels=variables);

pkb.cortests <- lapply(variables, function(v) {
	cor.test(pkb$purity, pkb[[v]], method="kendall", alternative="greater")
});

pkb.text <- data.frame(
	variable = variables,
	cortest = unlist(lapply(pkb.cortests, cortest_to_string))
);

qdraw(
	ggplot(pkb.tall, aes(x=purity, y=value)) +
		geom_point(alpha=0.5) + theme_line() +
		geom_text(data = pkb.text, aes(x=-Inf, y=Inf, label=cortest), hjust=-0.1, vjust=1.5) +
		stat_smooth(method = loess, span=10) +
		facet_wrap(~ variable, ncol=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = cna.labeller)
		) +
		xlim(0, 1) +
		ylab("Count of copy-number events") + xlab("Tumor purity")
	,
	height = 8,
	file = insert(pdf.fname, "pkb-bm-luad")
);

####

pkb.chrom.cna <- mutate(pkb.chrom.seg,
		gain = state > chrom.gain.cut,
		loss = state < chrom.loss.cut,
		cna = gain | loss
	) %>%
	group_by(sample) %>%
	summarize(
		count_cna = sum(cna),
		count_gain = sum(gain),
		count_loss = sum(loss)
	);

pkb.chrom <- left_join(pkb.pheno, pkb.chrom.cna, by = c("sample_id"="sample"));

pkb.chrom.sel <- select(pkb.chrom, sample_id, purity, count_cna, count_gain, count_loss, weight, weight_norm);

pkb.chrom.tall <- pkb.chrom.sel %>% gather(key="variable", value="value", -sample_id, -purity, -weight, -weight_norm)

pkb.chrom.tall$variable <- factor(pkb.chrom.tall$variable, levels=variables);


pkb.chrom.cortests <- lapply(variables.broad, function(v) {
	cor.test(pkb.chrom$purity, pkb.chrom[[v]], method="kendall", alternative="greater")
});

pkb.chrom.text <- data.frame(
	variable = variables.broad,
	cortest = unlist(lapply(pkb.chrom.cortests, cortest_to_string))
);

qdraw(
	ggplot(pkb.chrom.tall, aes(x=purity, y=value)) +
		geom_point(alpha=0.5) + theme_line() +
		geom_text(data = pkb.chrom.text, aes(x=-Inf, y=Inf, label=cortest), hjust=-0.1, vjust=1.5) +
		xlim(0, 1) +
		stat_smooth(method = loess, span=10) +
		facet_wrap(~ variable, ncol=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = cna.labeller)
		) +
		ylab("Count of copy-number events") + xlab("Tumor purity")
	,
	height = 4,
	file = insert(pdf.fname, c("pkb-bm-luad", "chrom"))
);

####

tcga.cna <- mutate(tcga.seg,
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

outliers <- filter(tcga.cna, count_cna > 500)$sample;
# no outliers

tcga.cna[tcga.cna$sample %in% outliers, ] <- NA;

tcga <- left_join(tcga.pheno, tcga.cna, by = c("sample_id"="sample"));
tcga.sel <- select(tcga, sample_id, purity, count_hcna, count_hamp, count_hdel, count_cna, count_gain, count_loss, weight, weight_norm);

tcga.tall <- tcga.sel %>% gather(key="variable", value="value", -sample_id, -purity, -weight, -weight_norm)

tcga.tall$variable <- factor(tcga.tall$variable, levels=variables);

tcga.cortests <- lapply(variables, function(v) {
	cor.test(tcga$purity, tcga[[v]], method="kendall", alternative="greater")
});

tcga.text <- data.frame(
	variable = variables,
	cortest = unlist(lapply(tcga.cortests, cortest_to_string))
);

qdraw(
	ggplot(tcga.tall, aes(x=purity, y=value)) +
		geom_point(alpha=0.5) + theme_line() +
		geom_text(data = tcga.text, aes(x=-Inf, y=Inf, label=cortest), hjust=-0.1, vjust=1.5) +
		stat_smooth(method = loess, span=10) +
		xlim(0, 1) +
		facet_wrap(~ variable, ncol=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = cna.labeller)
		) +
		ylab("Count of copy-number events") + xlab("Tumor purity")
	,
	height = 8,
	file = insert(pdf.fname, "tcga-luad")
)

####

tcga.chrom.cna <- mutate(tcga.chrom.seg,
		gain = state > chrom.gain.cut,
		loss = state < chrom.loss.cut,
		cna = gain | loss
	) %>%
	group_by(sample) %>%
	summarize(
		count_cna = sum(cna),
		count_gain = sum(gain),
		count_loss = sum(loss)
	);

tcga.chrom <- left_join(tcga.pheno, tcga.chrom.cna, by = c("sample_id"="sample"));

tcga.chrom.sel <- select(tcga.chrom, sample_id, purity, count_cna, count_gain, count_loss, weight, weight_norm);

tcga.chrom.tall <- tcga.chrom.sel %>% gather(key="variable", value="value", -sample_id, -purity, -weight, -weight_norm)

tcga.chrom.tall$variable <- factor(tcga.chrom.tall$variable, levels=variables);

tcga.chrom.cortests <- lapply(variables.broad, function(v) {
	cor.test(tcga.chrom$purity, tcga.chrom[[v]], method="kendall", alternative="greater")
});

tcga.chrom.text <- data.frame(
	variable = variables.broad,
	cortest = unlist(lapply(tcga.chrom.cortests, cortest_to_string))
);

qdraw(
	ggplot(tcga.chrom.tall, aes(x=purity, y=value)) +
		geom_point(alpha=0.5) + theme_line() +
		geom_text(data = tcga.chrom.text, aes(x=-Inf, y=Inf, label=cortest), hjust=-0.1, vjust=1.5) +
		stat_smooth(method = loess, span=10) +
		facet_wrap(~ variable, ncol=1, scales="free_y", strip.position = "top",
			labeller = labeller(variable = cna.labeller)
		) +
		ylab("Count of copy-number events") + xlab("Tumor purity") +
		xlim(0, 1) +
		ylim(0, 30)
	,
	height = 4,
	file = insert(pdf.fname, c("tcga-luad", "chrom"))
);

####

combined <- rbind(
	data.frame(tcga.sel, cohort = "TCGA-LUAD"),
	data.frame(pkb.sel, cohort = "BM-LUAD")
);

fit <- glm(count_hcna ~ cohort + purity, combined, family="quasipoisson");
summary(fit)

fit <- glm(count_hcna ~ cohort + purity, combined, family="quasipoisson", weights=combined$weight);
summary(fit)

fit <- glm(count_hcna ~ purity, combined, family="quasipoisson")
summary(fit)

fit <- glm(count_hcna ~ purity, combined, family="quasipoisson", weights=combined$weight);
summary(fit)

fit <- glm(count_cna ~ cohort + purity, combined, family="quasipoisson")
summary(fit)

fit <- glm(count_cna ~ cohort + purity, combined, family="quasipoisson", weights=combined$weight)
summary(fit)

fit <- glm(count_cna ~ purity, combined, family="quasipoisson")
summary(fit)

alphav <- 0.8;

qdraw(
	ggplot(combined, aes(x = cohort, y = count_cna + 1, fill=cohort)) +
		geom_violin(alpha=alphav) + theme_bw() + scale_y_log10() +
		scale_fill_manual(values=col.cohorts) +
		geom_jitter(alpha=0.3) +
		ylab("Count of events + 1") + xlab("") +
		ggtitle("Gains or losses")
	,
	width = 3.5,
	file = insert(pdf.fname, c("cna", "compare", "violin"))
);

wilcox.test(count_cna ~ cohort, combined)
wilcox.test(count_cna ~ cohort, combined, alternative="less")
t.test(count_cna ~ cohort, combined)


qdraw(
	ggplot(combined, aes(x = cohort, y = count_hcna + 1, fill=cohort)) +
		geom_violin(alpha=alphav) + theme_bw() + scale_y_log10() +
		scale_fill_manual(values=col.cohorts) +
		geom_jitter(alpha=0.3) +
		ylab("Count of events + 1") + xlab("") +
		ggtitle("Amplifications or deletions")
	,
	width = 3.5,
	file = insert(pdf.fname, c("hcna", "compare", "violin"))
);

wilcox.test(count_hcna ~ cohort, combined)
wilcox.test(count_hcna ~ cohort, combined, alternative="less")
t.test(count_hcna ~ cohort, combined)

####

alpha <- 0.05;
z <- qnorm(1 - alpha / 2);

fit0.hcna <- glm(count_hcna ~ cohort - 1, combined, family="quasipoisson", weights=combined$weight);
s <- summary(fit0.hcna);
hcna.means.df <- data.frame(
	cohort = factor(levels(combined$cohort), levels(combined$cohort)),
	mu = s$coefficients[,1],
	error = s$coefficients[,2]
) %>% mutate(
	y = exp(mu),
	ymin = exp(mu - z * error),
	ymax = exp(mu + z * error)
);

fit.hcna <- glm(count_hcna ~ cohort, combined, family="quasipoisson", weights=combined$weight);
p.hcna <- coef(summary(fit.hcna))[2,4];

qdraw(
	ggplot(hcna.means.df, aes(x=cohort, y=y, ymin=ymin, ymax=ymax, fill=cohort)) +
		geom_bar(stat="identity", alpha=alphav) + geom_errorbar(width=0.2) + theme_bw() +
		scale_fill_manual(values=col.cohorts) +
		annotate("text", x = 1, y = 3, label = p_to_string(p.hcna)) +
		ylab("Mean number of high-level copy-number events") + xlab("") +
		ggtitle("Amplifications or deletions")
	,
	width = 3.5,
	file = insert(pdf.fname, c("hcna", "compare", "bar"))
);

fit0.cna <- glm(count_cna ~ cohort - 1, combined, family="quasipoisson", weights=combined$weight);
s <- summary(fit0.cna);
cna.means.df <- data.frame(
	cohort = factor(levels(combined$cohort), levels(combined$cohort)),
	mu = s$coefficients[,1],
	error = s$coefficients[,2]
) %>% mutate(
	y = exp(mu),
	ymin = exp(mu - z * error),
	ymax = exp(mu + z * error)
);

fit.cna <- glm(count_cna ~ cohort, combined, family="quasipoisson", weights=combined$weight);
p.cna <- coef(summary(fit.cna))[2,4];

qdraw(
	ggplot(cna.means.df, aes(x=cohort, y=y, ymin=ymin, ymax=ymax, fill=cohort)) +
		geom_bar(stat="identity", alpha=alphav) + geom_errorbar(width=0.2) + theme_bw() +
		scale_fill_manual(values=col.cohorts) +
		annotate("text", x = 1, y = 50, label = p_to_string(p.cna)) +
		ylab("Mean number of copy-number events") + xlab("") +
		ggtitle("Gains or losses")
	,
	width = 3.5,
	file = insert(pdf.fname, c("cna", "compare", "bar"))
);

# GPLDIFF has a global difference parameter
# ABSOLUTE accounts for purity

####

combined.chrom <- rbind(
	data.frame(tcga.chrom.sel, cohort = "TCGA-LUAD"),
	data.frame(pkb.chrom.sel, cohort = "BM-LUAD")
);

fit <- glm(count_cna ~ cohort + purity, combined.chrom, family="quasipoisson")
summary(fit)

fit <- glm(count_cna ~ cohort + purity, combined.chrom, family="quasipoisson", weights=combined.chrom$weight)
summary(fit)

fit <- glm(count_cna ~ purity, combined.chrom, family="quasipoisson")
summary(fit)

alphav <- 0.8;

qdraw(
	ggplot(combined.chrom, aes(x = cohort, y = count_cna + 1, fill=cohort)) +
		geom_violin(alpha=alphav) + theme_bw() + scale_y_log10() +
		scale_fill_manual(values=col.cohorts) +
		geom_jitter(alpha=0.3) +
		ylab("Count of events + 1") + xlab("") +
		ggtitle("Chromosome gains or losses")
	,
	width = 3.5,
	file = insert(pdf.fname, c("chrom-cna", "compare", "violin"))
);

wilcox.test(count_cna ~ cohort, combined.chrom)
wilcox.test(count_cna ~ cohort, combined.chrom, alternative="less")
t.test(count_cna ~ cohort, combined.chrom)

####

alpha <- 0.05;
z <- qnorm(1 - alpha / 2);

fit0.cna <- glm(count_cna ~ cohort - 1, combined.chrom, family="quasipoisson", weights=combined.chrom$weight);
s <- summary(fit0.cna);
cna.means.df <- data.frame(
	cohort = factor(levels(combined$cohort), levels(combined$cohort)),
	mu = s$coefficients[,1],
	error = s$coefficients[,2]
) %>% mutate(
	y = exp(mu),
	ymin = exp(mu - z * error),
	ymax = exp(mu + z * error)
);

fit.cna <- glm(count_cna ~ cohort, combined.chrom, family="quasipoisson", weights=combined.chrom$weight);
p.cna <- coef(summary(fit.cna))[2,4];

qdraw(
	ggplot(cna.means.df, aes(x=cohort, y=y, ymin=ymin, ymax=ymax, fill=cohort)) +
		geom_bar(stat="identity", alpha=alphav) + geom_errorbar(width=0.2) + theme_bw() +
		scale_fill_manual(values=col.cohorts) +
		annotate("text", x = 1, y = 8, label = p_to_string(p.cna)) +
		ylab("Mean number of copy-number events") + xlab("") +
		ggtitle("Chromosome gains or losses")
	,
	width = 3.5,
	file = insert(pdf.fname, c("chrom-cna", "compare", "bar"))
);

####

graphics.off()

