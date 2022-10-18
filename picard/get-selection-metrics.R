library(io);
library(dplyr);

options(stringsAsFactors=FALSE);

# proportion power calculation for right-tailed test
# at significant level of 0.05 and 90% power
# need to distinguish heterzygous variant at tumour purity 20%
# with base quality of 20 (probability of error = 10^-(20/10) = 0.01)
# effect size is calculated by h = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
#library(pwr);
#h <- 2*asin(sqrt(0.5 * 0.2)) - 2*asin(sqrt(0.01));
#pwr.p.test(h, sig.level=0.05, power=0.90, alternative="greater")

pheno <- qread("../annot/sample-info_wes.tsv");

extract_stats <- function(file) {
	lapply(file, extract_sample_stats)
}

extract_sample_stats <- function(file) {
	message(file)
	if (!is.na(file) && file.exists(file)) {
		x <- read.table(file, comment.char="#", blank.lines.skip=TRUE, sep="\t", header=TRUE, nrows=1);
		if (nrow(x) > 0) {
			return(c(
				hs_on_bait_bases = x$ON_BAIT_BASES,
				hs_near_bait_bases = x$NEAR_BAIT_BASES,
				hs_off_bait_bases = x$OFF_BAIT_BASES,
				hs_on_target_bases = x$ON_TARGET_BASES,
				hs_bait_coverage = x$MEAN_BAIT_COVERAGE,
				hs_target_coverage = x$MEAN_TARGET_COVERAGE,
				hs_pct_bases_on_bait = x$PCT_USABLE_BASES_ON_BAIT,
				hs_pct_bases_on_target = x$PCT_USABLE_BASES_ON_TARGET,
				hs_fold_enrichment = x$FOLD_ENRICHMENT,
				hs_pct_targets_zero_cvg = x$ZERO_CVG_TARGETS_PCT,
				hs_library_size = x$HS_LIBRARY_SIZE
			));
		}
	}

	c(
		hs_on_bait_bases = NA,
		hs_near_bait_bases = NA,
		hs_off_bait_bases = NA,
		hs_on_target_bases = NA,
		hs_bait_coverage = NA,
		hs_target_coverage = NA,
		hs_pct_bases_on_bait = NA,
		hs_pct_bases_on_target = NA,
		hs_fold_enrichment = NA,
		hs_pct_targets_zero_cvg = NA,
		hs_library_size = NA
	)
}

pheno <- mutate(pheno,
	sample_directory = sub("/[^/]*\\.bam$", "", clean_bam_file_capture),
	bam_file_stem = sub(".*/([^/]+)\\.bam$", "\\1", clean_bam_file_capture),
	selection_metrics_file = file.path(sample_directory, paste0(bam_file_stem, ".hybrid_selection_metrics"))
);

stats.list <- extract_stats(pheno$selection_metrics_file);
stats <- matrix(unlist(stats.list), nrow=length(stats.list[[1]]));
rownames(stats) <- names(stats.list[[1]]);
stats <- data.frame(sample_id=pheno$sample_id, t(stats));

qwrite(stats, filename("selection-metrics", ext="tsv"));

