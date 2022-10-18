library(io)
library(dplyr)
library(ggplot2)

source("~/exigens/brain-mets/common/params.R");
source("theme.R");

# 1-based coordinates
pkb.bm.seg <- qread("../gistic/brain-mets_pass_luad_absolute_bmet-only.lgr.cntf.cenf.cnvf.seg");
tcga.seg <- qread("../gistic/tcga-luad_pass_absolute.lgr.cntf.cenf.cnvf.seg");

pkb.bm.seg$chromosome <- as.integer(pkb.bm.seg$chromosome);
tcga.seg$chromosome <- as.integer(tcga.seg$chromosome);

chrom.lens <- qread("chrom-lens_no-cen.vtr");

# chromosome Y is not assessed
chrom.lens.sum <- sum(chrom.lens[1:23]);

cna.cut <- 0.2;


pkb.bm.cna <- filter(pkb.bm.seg, abs(state) > cna.cut) %>%
	group_by(sample) %>% summarize(prop_cna = sum(end - start + 1) / chrom.lens.sum);

tcga.cna <- filter(tcga.seg, abs(state) > cna.cut) %>%
	group_by(sample) %>% summarize(prop_cna = sum(end - start + 1) / chrom.lens.sum);

d <- rbind(
	data.frame(pkb.bm.cna, cohort="BM-LUAD"),
	data.frame(tcga.cna, cohort="TCGA-LUAD")
);

alphav <- 0.8;

ggplot(d, aes(x=cohort, y=prop_cna, fill=cohort)) +
	geom_violin(alpha=alphav) + my_theme() +
	scale_fill_manual(values=col.cohorts)

ggplot(d, aes(x=cohort, y=prop_cna, colour=cohort)) +
	geom_jitter(alpha=alphav, width=0.1) + my_theme() +
	scale_colour_manual(values=col.cohorts)

wilcox.test(prop_cna ~ cohort, data = d, alternative="greater")
t.test(prop_cna ~ cohort, data = d, alternative="greater")


amp.cut <- 2;

pkb.bm.amp <- filter(pkb.bm.seg, abs(state) > amp.cut) %>%
	group_by(sample) %>% summarize(n_amp = n());

tcga.amp <- filter(tcga.seg, abs(state) > amp.cut) %>%
	group_by(sample) %>% summarize(n_amp = n());

d.amp <- rbind(
	data.frame(pkb.bm.amp, cohort="BM-LUAD"),
	data.frame(tcga.amp, cohort="TCGA-LUAD")
);

ggplot(d.amp, aes(x=cohort, y=n_amp, fill=cohort)) +
	geom_violin(alpha=alphav) + my_theme() +
	scale_fill_manual(values=col.cohorts)

ggplot(d.amp, aes(x=cohort, y=n_amp, colour=cohort)) +
	geom_jitter(alpha=alphav, width=0.1) + my_theme() +
	scale_colour_manual(values=col.cohorts)

wilcox.test(n_amp ~ cohort, data = d.amp, alternative="greater")
t.test(n_amp ~ cohort, data = d.amp, alternative="greater")



del.cut <- -3;

pkb.bm.del <- filter(pkb.bm.seg, abs(state) < -del.cut) %>%
	group_by(sample) %>% summarize(n_del = n());

tcga.del <- filter(tcga.seg, abs(state) < -del.cut) %>%
	group_by(sample) %>% summarize(n_del = n());

d.del <- rbind(
	data.frame(pkb.bm.del, cohort="BM-LUAD"),
	data.frame(tcga.del, cohort="TCGA-LUAD")
);

ggplot(d.del, aes(x=cohort, y=n_del, fill=cohort)) +
	geom_violin(alpha=alphav) + my_theme() +
	scale_fill_manual(values=col.cohorts)

ggplot(d.del, aes(x=cohort, y=n_del, colour=cohort)) +
	geom_jitter(alpha=alphav, width=0.1) + my_theme() +
	scale_colour_manual(values=col.cohorts)

wilcox.test(n_del ~ cohort, data = d.del, alternative="greater")

t.test(n_del ~ cohort, data = d.del, alternative="greater")

