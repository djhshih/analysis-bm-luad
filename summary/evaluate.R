library(io);
library(dplyr);
library(ggplot2);

options(stringsAsFactors=FALSE);

count.snvs <- qread("count-snvs.tsv");
count.indels <- qread("count-indels.tsv");

agilent.sample.stats <- qread("../fiss/brain-mets_agilent-clean_sample-stats.tsv");
agilent.pair.stats <- qread("../fiss/brain-mets_agilent-clean_pair-stats.tsv");
ice.sample.stats <- qread("../fiss/brain-mets_ice-20160822_sample-stats.tsv");
ice.pair.stats <- qread("../fiss/brain-mets_ice-20160822_pair-stats.tsv");

clin <- qread("~/exigens/brain-mets/annot/patient-info.tsv");

pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2.tsv");

x <- pheno %>% 
	inner_join(rename(count.snvs, count_snvs = count)) %>%
	inner_join(rename(count.indels, count_indels = count)) %>%
	left_join(rbind(agilent.sample.stats, ice.sample.stats), by = c("fh_sample_id" = "sample_id")) %>%
	left_join(rbind(agilent.pair.stats, ice.pair.stats), by = c("fh_pair_id" = "pair_id")) %>%
	rename(somatic_mutation_covered_bases_capture = somatic_mutation_covered_bases_capture.x);

x <- x %>% mutate(
	covered_mb = (somatic_mutation_covered_bases_capture / 1e6),
	mutation_burden = count_snvs / covered_mb,
	mutation_burden_indel = count_indels / covered_mb
);

luad <- x %>% filter(primary_histotype == "Lung adenocarcinoma");

out <- select(x, sample_id, count_snvs, count_indels, mutation_burden, mutation_burden_indel);
qwrite(out, "burden.tsv");

##

qdraw(
	ggplot(x, aes(x=wes_platform, colour=specimen_material, y=count_snvs)) +
		geom_jitter(width=0.3) + scale_y_log10(),
	file = "counts_platform.pdf"
);

qdraw(
	ggplot(
		filter(x, primary_cancer_type == "Lung cancer"),
		aes(x=wes_platform, colour=specimen_material, y=count_snvs)) +
		geom_jitter(width=0.3) + scale_y_log10(),
	file = "counts_platform_lung.pdf"
);

wilcox.test(count_snvs ~ wes_platform, x)
t.test(count_snvs ~ wes_platform, x)

qdraw(
	ggplot(x, aes(x=specimen_material, y=count_snvs)) +
		geom_jitter(width=0.3) + scale_y_log10(),
	file = "counts_specimen.pdf"
);

wilcox.test(count_snvs ~ specimen_material, x)
anova(lm(count_snvs ~ specimen_material, x))

with(x, table(specimen_material, useNA="always"));
with(filter(x, primary_cancer_type == "Lung cancer"), table(specimen_material, useNA="always"));

with(x, table(wes_platform, useNA="always"));

qdraw(
	ggplot(x, aes(x=primary_cancer_type, colour=sample_type, y=count_snvs)) +
		geom_jitter(width=0.3) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "counts_primary-cancer-type.pdf",
	width = 10
);


##

qdraw(
	ggplot(x, aes(x=somatic_mutation_covered_bases_capture, colour=wes_platform, y=count_snvs)) +
		geom_point(alpha=0.5) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "counts_coverage_platform.pdf",
	width = 8
);

qdraw(
	ggplot(x, aes(x=somatic_mutation_covered_bases_capture, colour=wes_platform, y=count_indels)) +
		geom_point(alpha=0.5) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "counts-indels_coverage_platform.pdf",
	width = 8
);

qdraw(
	ggplot(x, aes(x=somatic_mutation_covered_bases_capture, colour=center, y=count_snvs)) +
		geom_point(alpha=0.5) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "counts_coverage_center.pdf",
	width = 8
);

qdraw(
	ggplot(x, aes(x=somatic_mutation_covered_bases_capture, colour=center, y=count_indels)) +
		geom_point(alpha=0.5) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "counts-indels_coverage_center.pdf",
	width = 8
);

##


qdraw(
	ggplot(x, aes(x=wes_platform, colour=specimen_material, y=mutation_burden)) +
		geom_jitter(width=0.3) + scale_y_log10(),
	file = "burden_platform.pdf"
);

qdraw(
	ggplot(x, aes(x=specimen_material, y=mutation_burden)) +
		geom_jitter(width=0.3) + scale_y_log10(),
	file = "burden_specimen.pdf"
);

qdraw(
	ggplot(x, aes(x=primary_cancer_type, colour=sample_type, y=mutation_burden)) +
		geom_jitter(width=0.3) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "burden_primary-cancer-type.pdf",
	width = 10
);

qdraw(
	ggplot(x, aes(x=primary_cancer_type, colour=sample_type, y=mutation_burden_indel)) +
		geom_jitter(width=0.3) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "burden-indels_primary-cancer-type.pdf",
	width = 10
);

qdraw(
	ggplot(
		mutate(x, genome_doublings=factor(genome_doublings)),
			aes(x=purity, y=mutation_burden, colour=genome_doublings)) +
		geom_point() + scale_y_log10()
	,
	file = "burden_purity.pdf",
);

qdraw(
	{
		hist(log(x$mutation_burden), breaks=20, col="grey60", border=NA, freq=FALSE);
		with(x, curve(dnorm(x, mean(log(mutation_burden)), sd(log(mutation_burden))), col="red", add=TRUE));
	},
	file = "burden_hist.pdf"
);

shapiro.test(log(x$mutation_burden))

summary(lm(log(mutation_burden) ~ purity + genome_doublings, x))
anova(lm(log(mutation_burden) ~ purity + genome_doublings, x))

##

qdraw(
	ggplot(luad, aes(x=sample_type, y=mutation_burden)) +
		geom_jitter(width=0.3) + scale_y_log10() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)),
	file = "burden_luad_sample-type.pdf",
	width = 10
);


qdraw(
	ggplot(
		mutate(luad, genome_doublings=factor(genome_doublings)),
			aes(x=purity, y=mutation_burden, colour=genome_doublings)) +
		geom_point() + scale_y_log10()
	,
	file = "burden_luad_purity.pdf",
);

##

y <- left_join(x, clin, by="clinical_id");

y.met <- filter(y, sample_type == "Brain metastasis") %>%
	rename(
		primary_cancer_type = primary_cancer_type.x,
		primary_histotype = primary_histotype.x
	);

y.met$primary_mutation_burden = unlist(lapply(
	y.met$clinical_id,
	function(id) {
		mean(y$mutation_burden[which(y$sample_type == "Primary" & y$clinical_id == id)])
	}
));

y.met$primary_count_snvs = unlist(lapply(
	y.met$clinical_id,
	function(id) {
		mean(y$count_snvs[which(y$sample_type == "Primary" & y$clinical_id == id)])
	}
));


with(y.met, wilcox.test(mutation_burden, primary_mutation_burden, paired=TRUE));
with(y.met, t.test(log(mutation_burden), log(primary_mutation_burden), paired=TRUE));

qdraw(
	ggplot(
		y.met,
			aes(x=bm_pfs_years, y=mutation_burden / primary_mutation_burden)) +
		geom_point() + scale_y_log10()
	,
	file = "burden_bm-pfs.pdf",
);

summary(lm(log(mutation_burden) - log(primary_mutation_burden) ~ bm_pfs_years, y.met))


qdraw(
	ggplot(
		y.met,
			aes(x = primary_mutation_burden, y = mutation_burden, colour=primary_cancer_type)) +
		geom_point() + scale_y_log10() + scale_x_log10() +
		geom_abline(slope=1, colour="grey40", alpha=0.5, size=1.5)
	,
	file = "burden_met-primary.pdf",
);

qdraw(
	ggplot(
		y.met,
			aes(x = primary_count_snvs, y = count_snvs, colour=primary_cancer_type)) +
		geom_point() + scale_y_log10() + scale_x_log10() +
		geom_abline(slope=1, colour="grey40", alpha=0.5, size=1.5)
	,
	file = "counts_met-primary.pdf",
);


qdraw(
	ggplot(
		y.met,
			aes(x = primary_mutation_burden, y = mutation_burden, colour = primary_cancer_type)) +
		geom_point() + scale_y_log10() + scale_x_log10() +
		geom_abline(slope=1, colour="grey40", alpha=0.5, size=1.5)
	,
	file = "burden_met-primary.pdf",
);

y.met.brca <- filter(y.met, primary_cancer_type == "Breast cancer");

summary(lm(log(mutation_burden) - log(primary_mutation_burden) ~ bm_pfs_years, y.met.brca));

summary(
	with(y.met.brca, log(mutation_burden) - log(primary_mutation_burden))
)

with(y.met.brca, wilcox.test(mutation_burden, primary_mutation_burden, paired=TRUE));
with(y.met.brca, t.test(log(mutation_burden), log(primary_mutation_burden), paired=TRUE));


y.met.luad <- filter(y.met, primary_histotype == "Lung adenocarcinoma");

qdraw(
	ggplot(
		y.met.luad,
			aes(x=bm_pfs_years, y=mutation_burden / primary_mutation_burden)) +
		geom_point() + scale_y_log10()
	,
	file = "burden_luad_bm-pfs.pdf",
);

with(y.met.luad, wilcox.test(mutation_burden, primary_mutation_burden, paired=TRUE));
with(y.met.luad, t.test(log(mutation_burden), log(primary_mutation_burden), paired=TRUE));

summary(lm(log(mutation_burden) - log(primary_mutation_burden) ~ bm_pfs_years, y.met.luad));


qdraw(
	ggplot(
		y.met.luad,
			aes(x = primary_mutation_burden, y = mutation_burden, colour = specimen_material)) +
		geom_point() + scale_y_log10() + scale_x_log10() +
		geom_abline(slope=1, colour="grey40", alpha=0.5, size=1.5)
	,
	file = "burden_luad_met-primary.pdf",
);

