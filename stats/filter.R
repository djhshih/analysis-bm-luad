library(io);
library(dplyr);
library(ggplot2);

options(stringsAsFactors=FALSE);

pheno <- qread("../annot/sample-info_wes.tsv");

out.fname <- filename("stats");

picard.only <- qread("../picard/duplicate-metrics.tsv");

pairs.agilent <- qread("../fiss/brain-mets_agilent-clean_pair-stats.tsv");
pairs.ice <- qread("../fiss/brain-mets_ice-20160822_pair-stats.tsv");

pairs <- rbind(
		mutate(pairs.agilent, wes_platform="Agilent"),
		mutate(pairs.ice, wes_platform="ICE")
	) %>% left_join(pheno, by=c("pair_id"="fh_pair_id")) %>%
	rename(wes_platform = wes_platform.x) %>%
	mutate(
		ffpe_reject_pct_capture = ffpe_reject_count_capture / (ffpe_reject_count_capture + ffpe_pass_count_capture) * 100
	);

samples.agilent <- qread("../fiss/brain-mets_agilent-clean_sample-stats.tsv");
samples.ice <- qread("../fiss/brain-mets_ice-20160822_sample-stats.tsv");

samples <- rbind(
		mutate(samples.agilent, wes_platform="Agilent"),
		mutate(samples.ice, wes_platform="ICE")
	) %>% left_join(pheno, by=c("sample_id"="fh_sample_id")) %>%
	rename(wes_platform = wes_platform.x);

wilcox.test(contamination_percentage_consensus_capture ~ wes_platform, samples)
t.test(contamination_percentage_consensus_capture ~ wes_platform, samples)


# FF samples with SNVs filtered due to FFPE filter
ff.rejected <- filter(pairs, ffpe_reject_count_capture > 0, specimen_material == "FF");

qdraw(
	ggplot(pairs, aes(x=specimen_material, colour=wes_platform, y=ffpe_reject_count_capture)) +
		geom_jitter(width=0.5, alpha=0.5)
	,
	file = insert(out.fname, tag="ffpe-reject-count", ext="pdf")
);

qdraw(
	ggplot(mutate(pairs, label = ifelse(pair_id %in% ff.rejected$pair_id, sample_id, NA)),
		aes(x=specimen_material, colour=wes_platform, y=ffpe_reject_pct_capture, label=label)) +
			geom_jitter(width=0.5, alpha=0.5) +
			geom_text(size=2, colour="grey10", nudge_x=-0.4)
	,
	file = insert(out.fname, tag="ffpe-reject-pct", ext="pdf")
);

qdraw(
	ggplot(pairs, aes(x=ifelse(is.na(specimen_material), "NA", specimen_material), colour=wes_platform, y=ffpe_Q)) +
		geom_jitter(width=0.5, alpha=0.5)
	,
	file = insert(out.fname, tag="ffpe-q", ext="pdf")
);
# higher => fewer expected artificats attributable to FFPE damage

filter(pairs, ffpe_Q < 40, specimen_material == "FF");
filter(pairs, ffpe_Q >= 40 & ffpe_Q < 65, specimen_material == "FF");

filter(pairs, ffpe_Q > 65, specimen_material == "FFPE");


qdraw(
	ggplot(pairs, aes(x=ifelse(is.na(specimen_material), "NA", specimen_material), colour=wes_platform, y=picard_oxoQ)) +
		geom_jitter(width=0.5, alpha=0.5)
	,
	file = insert(out.fname, tag="picard-oxoq", ext="pdf")
);
# higher => fewer expected artifacts attributable to oxidative damage


samples <- samples %>% mutate(
	cross_contaminated = contamination_percentage_consensus_capture > 5
);

qdraw(
	ggplot(samples, aes(x=specimen_material, colour=wes_platform, y=contamination_percentage_consensus_capture)) +
		geom_jitter(width=0.5, alpha=0.5),
	file = insert(out.fname, tag="cross-contamination", ext="pdf"),
	width = 5
);

qdraw(
	ggplot(samples, aes(x=contamination_percentage_consensus_capture)) +
		geom_histogram()
	,
	file = insert(out.fname, tag=c("cross-contamination", "hist"), ext="pdf")
);

logit <- function(x) {
	log(x) - log(1 - x)
}

qdraw(
	ggplot(samples, aes(x=logit(contamination_percentage_consensus_capture/100))) +
		geom_histogram()
	,
	file = insert(out.fname, tag=c("cross-contamination", "logit", "hist"), ext="pdf")
);

samples.ccon <- filter(samples, cross_contaminated) %>%
	select(sample_id.y, sample_id, fh_pair_id, contamination_percentage_consensus_capture, primary_cancer_type, primary_histotype) %>%
	rename(fh_sample_id = sample_id, sample_id = sample_id.y)

qwrite(samples.ccon, insert(out.fname, tag=c("cross-contamination", "platform"), ext="tsv"));

# Samples with no cross-contamination estimate (due to too few homozygous loci
# with read mismatches)
filter(samples, contamination_percentage_consensus_capture < 0)

qdraw(
	ggplot(samples, aes(x=specimen_material, colour=center, y=contamination_percentage_consensus_capture)) +
		geom_jitter(width=0.5, alpha=0.5)
	,
	file = insert(out.fname, tag=c("cross-contamination", "center"), ext="pdf")
);


pairs <- pairs %>% mutate(
	tumor_in_normal = combined_deTiN_TiN_num > 0.10
);

pairs.detin <- filter(pairs, tumor_in_normal) %>%
	select(sample_id, pair_id, combined_deTiN_TiN_num, primary_cancer_type, primary_histotype) %>%
	rename(fh_pair_id = pair_id);
# Although PB0091-L2 and PB0083-M have high tumour-in-normal contamination,
# these contamination estimates are inconsistent with the copy-number changes
# observed in the normal samples. Therefore, the estimates are suspect.

qdraw(
	ggplot(pairs, aes(x=specimen_material, colour=wes_platform, y=combined_deTiN_TiN_num * 100)) +
		geom_jitter(width=0.5, alpha=0.5)
	,
	file = insert(out.fname, tag=c("detin", "platform"), ext="pdf")
);

qdraw(
	ggplot(pairs, aes(x=specimen_material, colour=center, y=combined_deTiN_TiN_num * 100)) +
		geom_jitter(width=0.5, alpha=0.5)
	,
	file = insert(out.fname, tag=c("detin", "center"), ext="pdf")
);

qdraw(
	ggplot(pairs, aes(x=somatic_mutation_covered_bases_capture / 1e6)) +
		geom_histogram() + facet_grid(specimen_material ~ wes_platform)
	,
	file = insert(out.fname, tag=c("covered-bases"), ext="pdf")
);

# exome target regions are ~30e6 bp
# bases have to be covered by sufficient number of *non-duplicate* reads
filter(pairs, somatic_mutation_covered_bases_capture < 10e6);
# BS-528-P-3 has 95% read duplication, 2e6 unique molecules
# PB0145-P has 94% read duplication, 4e6 unique molecules
# PB0145-M looks fine: low duplication rate, 100e6 unique molecules in library
# but PB0145-N is low quality (80% read duplication, 1e6 unique molecules)

pairs$underpowered_coverage <- pairs$somatic_mutation_covered_bases_capture < 10e6;


picard <- left_join(pheno, picard.only);

qdraw(
	ggplot(picard, aes(x=duplication_rate)) +
		geom_histogram() + facet_grid(specimen_material ~ wes_platform)
	,
	file = insert(out.fname, tag=c("dup"), ext="pdf")
);

qdraw(
	ggplot(picard, aes(x=estimated_library_size)) +
		geom_histogram() + facet_grid(specimen_material ~ wes_platform)
	,
	file = insert(out.fname, tag=c("lib-size"), ext="pdf")
);

qdraw(
	ggplot(picard, aes(x=duplication_rate, y=estimated_library_size, colour=specimen_material, shape=wes_platform)) +
		geom_point()
	,
	file = insert(out.fname, tag=c("lib-size-vs-dup"), ext="pdf")
);

pairs <- left_join(pairs, picard.only) %>%
	rename(tumor_duplication_rate = duplication_rate, tumor_library_size = estimated_library_size);

idx <- match(pheno$sample_id[match(pairs$control_sample, pheno$fh_sample_id)], picard.only$sample_id);
pairs$normal_duplication_rate <- picard.only$duplication_rate[idx];
pairs$normal_library_size <- picard.only$estimated_library_size[idx];

qdraw(
	ggplot(pairs, aes(x=log10(tumor_library_size), y=log10(normal_library_size), colour=log10(somatic_mutation_covered_bases_capture))) +
		geom_point(alpha=0.5) + geom_text(aes(label=ifelse(underpowered_coverage, sample_id, NA)), angle=45, hjust=-0.15, colour="black", size=2) +
		theme_bw() +
		geom_text(aes(label=ifelse(log10(normal_library_size) < 4, sample_id, NA)), angle=45, hjust=-0.15, colour="green4", size=2)
	,
	file = insert(out.fname, tag=c("lib-size-vs-covered-bases"), ext="pdf"),
	width=10
);

