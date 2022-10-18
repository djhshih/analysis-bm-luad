library(io);
library(dplyr);

x <- qread("brain-mets_agilent-master_bams.tsv");

fix_bad_bam_path <- function(x) {
	x %>%
		sub("tier3b/", "", .) %>%
		sub("picard_aggregation2", "picard_aggregation", .)
}

y <- x %>%
	filter(grepl("tier3b/", clean_bam_file_capture) | grepl("picard_aggregation2/", clean_bam_file_capture));

y <- y %>%
	mutate(clean_bam_file_capture=fix_bad_bam_path(clean_bam_file_capture));

qwrite(y, "brain-mets_agilent-master_bams_to-fix.tsv");

