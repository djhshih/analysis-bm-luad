library(io);
library(dplyr);

options(stringsAsFactors=FALSE);

pheno <- qread("../annot/sample-info_wes.tsv");

extract_duplication_stats <- function(file) {
	sapply(file, extract_single_duplication_stats)
}

extract_single_duplication_stats <- function(file) {
	if (is.na(file) || !file.exists(file)) {
		list(duplication_rate=NA, estimated_library_size=NA)
	} else {
		x <- scan(file, what=character(), comment.char="#", blank.lines.skip=TRUE, nlines=9, sep="\n");
		sep <- "\t";
		ss <- strsplit(x, sep);
		y <- as.list(ss[[2]]);
		names(y) <- ss[[1]];

		list(duplication_rate=as.numeric(y$PERCENT_DUPLICATION), estimated_library_size=as.numeric(y$ESTIMATED_LIBRARY_SIZE))
	}
}

pheno <- mutate(pheno,
	sample_directory = sub("/[^/]*\\.bam$", "", clean_bam_file_capture),
	bam_file_stem = sub(".*/([^/]+)\\.bam$", "\\1", clean_bam_file_capture),
	duplicate_metrics_file = file.path(sample_directory, paste0(bam_file_stem, ".duplicate_metrics"))
);

n <- nrow(pheno);
duplication.rates <- numeric(n);
estimated.library.sizes <- integer(n);
for (i in 1:n) {
	x <- pheno$duplicate_metrics_file[i];
	message(x);
	y <- extract_single_duplication_stats(x);
	duplication.rates[i] <- y$duplication_rate;
	estimated.library.sizes[i] <- y$estimated_library_size;
}

stats <- data.frame(
	sample_id = pheno$sample_id,
	duplication_rate = duplication.rates,
	estimated_library_size = estimated.library.sizes
);

qwrite(stats, filename("duplicate-metrics", ext="tsv"));

