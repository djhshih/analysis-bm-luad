library(io);
library(ggplot2);
library(ggrepel);

source("~/exigens/brain-mets/common/params.R");

#case.fname <- "mutsig2cv-sig-genes_pkb-luad_brain-mets_black-f_coding.tsv";
#control.fname <- "mutsig2cv-sig-genes_tcga-luad_pass_black-f_coding.tsv";

case.fname <- "mutsig2cv-sig-genes_pkb-luad_brain-mets_black-f_coding_orxnxf.tsv";
control.fname <- "mutsig2cv-sig-genes_tcga-luad_pass_black-f_coding_orxnxf.tsv";

out.fname <- filename("pkb-luad");

case <- qread(case.fname, quote="");
control <- qread(control.fname, quote="");

case[case$q < 0.25, ]

genes <- intersect(case$gene, control$gene);
G <- length(genes);

case <- case[match(genes, case$gene), ];
control <- control[match(genes, control$gene), ];

#plot(-log10(control$q), -log10(case$q))
#plot(-log10(control$p), -log10(case$p))

control$prank <- order(control$p, runif(G));
case$prank <- order(case$p, runif(G));

#show.label <- pmin(control$q, case$q) < 0.25;
#show.label <- pmin(control$rank, case$rank) < 20;
show.label <- case$q < 0.25 | control$q < 0.05;

max.rank <- max(case$rank, control$rank);
probs <- (max.rank - pmin(case$rank, control$rank)) / max.rank;
show.point <- show.label | rbinom(G, 1, prob=probs^100);

d <- data.frame(
	x = -log10(control$rank),
	y = -log10(case$rank),
	q = -log10(pmin(control$q, case$q)),
	significant = factor(case$q < 0.25, label=c(expression("q >= 0.25"), expression("q < 0.25"))),
	gene = ifelse(show.label, genes, NA),
	show = show.point
);

labs <- paste0("-log10(rank) in ", cohorts);

qdraw(
	ggplot(d, aes(x=x, y=y, label=gene)) +
		stat_density2d(aes(fill=..density..^0.2), geom="raster", contour=FALSE) +
		scale_fill_gradient(low="white", high="royalblue2") +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=x, y=y), data=d[d$show, ], size=1) +
		#geom_point(aes(x=x, y=y, size=q), data=d[d$show, ]) +
		#scale_size(range=c(0.1, 2)) +
		geom_label_repel(aes(colour=significant), fontface="italic", size=4.1) +
		scale_colour_manual(values=c("dodgerblue4", "firebrick")) +
		labs(x = labs[1], y = labs[2], fill="density", size="-log(q)", colour="BM-LUAD") +
		theme_bw() +
		theme(
			axis.text=element_text(size=12),
			axis.title=element_text(size=14),
			legend.text=element_text(size=12),
			legend.title=element_text(size=12)
		)
	,
	height = 6,
	width = 7.2,
	file = insert(out.fname, c("compare", "mutsig2cv-rank"), ext="pdf")
);

sig.genes <- as.character(case$gene[case$q < 0.25]);
qwrite(sig.genes, insert(out.fname, c("sig-genes"), ext="txt"));

