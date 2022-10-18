library(io);
library(dplyr);

x <- qread("../sample-info_wes_stage2.tsv")

qwrite(as.character(na.omit(x$clean_bam_file_capture)), "sample-bams.vtr")

y <- select(x, sample_id=fh_sample_id, clean_bam_file_capture);
qwrite(y, "sample-bams.tsv");

agilent.clean <- filter(x, wes_platform == "Agilent", !is.na(clean_bam_file_capture));

y2 <- agilent.clean %>%
	select(sample_id=fh_sample_id, clean_bam_file_capture);
qwrite(y2, "sample-bams_agilent.tsv");

y3 <- agilent.clean %>%
	filter(!is.na(fh_pair_id)) %>%
	transmute(pair_set_id = "clean_tumor_normal_pair_set", pair_id = fh_pair_id);
qwrite(y3, "clean_tumor-normal_pair-set.tsv");

