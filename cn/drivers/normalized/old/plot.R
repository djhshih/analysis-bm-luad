library(io);
library(dplyr);
library(ggplot2);
library(reshape2);

cnames <- c("chromosome", "start", "stop", "name", "value");

prepare <- function(d) {
	names(d) <- cnames;
	d$position <- ceiling((d$start + d$stop) / 2);
	d
}

combine_cn <- function(xs) {
	y <- data.frame(xs[[1]], id = names(xs)[1]);
	for (i in 2:length(xs)) {
		y <- rbind(y, data.frame(xs[[i]], id = names(xs)[i]));
	}
	y
}

lratio_to_linear <- function(x) {
	2^(x + 1)
}

mark_region <- function(x, region.start, region.end) {
	mutate(x, marked = position >= region.start & position <= region.end)
}

region_plot <- function(cn.l, patient, chrom, window.lim, region.lim) {
	window.start <- window.lim[1];
	window.end <- window.lim[2];
	region.start <- region.lim[1];
	region.end <- region.lim[2];

	patient.d <- combine_cn(cn.l) %>% 
		filter(chromosome == chrom, start >= window.start, stop <= window.end) %>%
		mark_region(region.start, region.end);

	ggplot(patient.d, aes(x = position / 1e6, y = lratio_to_linear(value), colour=marked)) +
		geom_hline(yintercept=2, linetype="dotted") +
		geom_point(alpha=0.3) + theme_bw() + 
		scale_colour_manual(values=c("black", "red")) +
		facet_wrap(~ id, ncol=1, strip.position="top") +
		xlab("Genomic position (Mbp)") + ylab("Raw copy-number") +
		ggtitle(sprintf("Patient %s on chr%s", patient, chrom)) +
		theme(
			legend.position = "none",
			strip.background = element_blank(),
			strip.text = element_text(hjust = 0)
		)
}

pdf.fname <- filename("cnplot", ext="pdf");

patient <- "LS-020";
cn.l <- list(
	"Normal" = qread("LS-020-N.tn.tsv") %>% prepare(),
	"Primary" = qread("LS-020-P.tn.tsv") %>% prepare(),
	"BM-1" = qread("LS-020-M-1.tn.tsv") %>% prepare(),
	"BM-2" = qread("LS-020-M-2.tn.tsv") %>% prepare(),
	"BM-3" = qread("LS-020-M-3.tn.tsv") %>% prepare()
);

# YAP1/MMP13
qdraw(
	region_plot(cn.l, patient, "11", c(90e6, 110e6), c(101981192, 102826463)),
	width = 6,	
	height = 10,
	file = insert(pdf.fname, c(patient, "yap1-mmp13"))
);

# CKDN2A/CDKN2B
qdraw(
	region_plot(cn.l, patient, "9", c(10e6, 30e6), c(21967751, 22009312)) + 
		scale_colour_manual(values=c("black", "royalblue")) + ylim(0, 3)
	,
	width = 6,
	height = 10,
	file = insert(pdf.fname, c(patient, "cdkn2ab"))
);


patient <- "PB0321";
cn.l <- list(
	"Normal" = qread("PB0321-N.tn.tsv") %>% prepare(),
	"Primary" = qread("PB0321-P.tn.tsv") %>% prepare(),
	"BM" = qread("PB0321-M.tn.tsv") %>% prepare()
);

# YAP1/MMP13
qdraw(
	region_plot(cn.l, patient, "11", c(90e6, 110e6), c(101981192, 102826463)),
	width = 6,	
	height = 6,
	file = insert(pdf.fname, c(patient, "yap1-mmp13"))
);


patient <- "PB0034";
cn.l <- list(
	"Normal" = qread("PB0034-N.tn.tsv") %>% prepare(),
	"Primary" = qread("PB0034-P.tn.tsv") %>% prepare(),
	"BM" = qread("PB0034-M.tn.tsv") %>% prepare()
);

qdraw(
	region_plot(cn.l, patient, "8", c(120e6, 140e6), c(128748315, 128753680)),
	width = 6,	
	height = 6,
	file = insert(pdf.fname, c(patient, "myc"))
);

