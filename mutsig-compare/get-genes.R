library(io);
library(dplyr);

x <- qread("~/exigens/brain-mets/mutsig-compare/pkb-luad_brain-mets_black-f_coding.pset.mmaf", type="tsv");
pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2.tsv");

genes.all <- qread("~/exigens/brain-mets/mutsig-compare/pkb-luad_sig-genes.txt", type="vtr");
exclude <- c("TP53", "KRAS", "EGFR", "CDKN2A", "STK11", "KEAP1", "NF1", "ARID1A", "APC");

genes <- setdiff(genes.all, exclude);
length(genes)

y <- filter(x, Hugo_Symbol %in% genes) %>% arrange(Hugo_Symbol, Start_Position) %>%
	mutate(sample_id = pheno$sample_id[match(Tumor_Sample_Barcode, pheno$fh_sample_id)]) %>%
	select(-Tumor_Sample_Barcode);

samples <- as.character(unique(y$sample_id));

# merge data from absolute
absolute.path <- "~/share/exigens/brain-mets/absolute/ABSOLUTE_results/brain-mets_pass_luad_absolute-1-4/reviewed/SEG_MAF";

sample.abs.paths <- file.path(absolute.path, sprintf("%s_ABS_MAF.txt", samples));
# NB  I/O intensive
absolutes <- lapply(sample.abs.paths, qread, type="tsv");

absolutes.sel <- mapply(
	function(samp, d) {
		transmute(d, Chromosome, Start_Position = Start_position, sample_id = samp, pr_germline = Pr_germline, pr_wt0 = Pr_wt0, ccf_hat, ccf_lower = ccf_CI95_low, ccf_upper = ccf_CI95_high, detection_power, n_ref_count, n_alt_count, t_ref_count = ref, t_alt_count = alt, ref_context)
	},
	samples, absolutes, SIMPLIFY=FALSE
);

absolute.df <- do.call(rbind, absolutes.sel);

y2 <- left_join(y, absolute.df, by = c("Chromosome", "Start_Position", "sample_id"));

# ensure that no mutation got duplicated during joining
stopifnot(nrow(y) == nrow(y2));

qwrite(y2, filename("pkb-luad", tag=c("brain-mets", "black-f", "coding", "sig-genes"), ext="tsv"));

