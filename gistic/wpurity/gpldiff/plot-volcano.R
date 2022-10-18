library(ggplot2);
library(ggrepel);
library(io);

results.fname <- "pkb-tcga-luad_gpldiff.tsv";
out.fname <- filename("pkb-tcga-luad", tag=c("wpurity", "gpldiff"));

#library(gpldiff)
#fits.fname <- "pkb-tcga-luad_gpldiff.rds";
#fits <- qread(fits.fname);
#results <- summary.gistic_gpldiffs(fits);

results <- qread(results.fname);

idx <- which(results$ldiff > 0.5 & results$fdr < 0.05 & results$gistic_q < 0.025);
print(results[idx, ])

# FIXME fragile!
# label significant results
results$cytoband <- NA;
results$cytoband[idx] <- c("11q22.2", NA, "8q24.21", NA,"9p21.3");

rev_log10_trans <- function() {
	scales::trans_new(
		name = "rev_log10",
		transform = function(x) -log10(x),
		inverse = function(x) 10^-x,
		breaks = function(x) c(0.1, 0.05, 0.01, 0.005, 0.001),
		domain = c(1e-100, Inf)
	)
}

qdraw(
	{
		ggplot(results, aes(x=diff, y=fdr, colour=type, label=cytoband)) +
			geom_vline(xintercept=0.5, colour="grey40", linetype="dashed") +
			geom_hline(yintercept=0.05, colour="grey40", linetype="dashed") +
			geom_point() + geom_label_repel() + theme_bw() +
			theme(panel.grid=element_blank()) +
			scale_colour_manual(values=c("#DE2D26", "#3182BD")) +
			scale_y_continuous(trans=rev_log10_trans()) +
			xlab("GISTIC score difference") + ylab("False discovery rate")
	},
	height = 4,
	file = insert(out.fname, tag="volcano", ext="pdf")
);

