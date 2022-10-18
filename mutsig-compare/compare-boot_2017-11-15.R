library(io);
library(ggplot2);
library(ggrepel);
library(reshape2);
library(dplyr);

case.fname <- "mutsig2cv-sig-genes_pkb-luad_brain-mets_black-f_coding.tsv";
control.fname <- "mutsig2cv-sig-genes_tcga-luad_pass_black-f_coding.tsv";
boot.indir <- "bootstrap";

out.fname <- filename("pkb-luad");

case <- qread(case.fname, quote="");
control <- qread(control.fname, quote="");
boot.orig <- qread(boot.indir, quote="");

genes <- union(case$gene, control$gene);
G <- length(genes);
B <- length(boot.orig);
sig.idx <- which(case$q < 0.25);
sig.genes <- genes[sig.idx];

case <- case[match(genes, case$gene), ];
control <- control[match(genes, control$gene), ];


# rearrnage boot to match gene order in case
boot <- lapply(
	boot.orig,
	function(x) {
		x[match(genes, x$gene), ]
	}
);

# store realized samples of case and control
case.r <- case;
control.r <- control;



rank.diff <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$rank - case$rank
	}
)));
rownames(rank.diff) <- genes;

hist(rank.diff, breaks=100)

logp.diff <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		-log10(case$p) - (-log10(x$p))
	}
)));
rownames(logp.diff) <- genes;

hist(logp.diff, breaks=100)

logp.diff.mean <- rowMeans(logp.diff, na.rm=TRUE);

hist(logp.diff.mean, breaks=100)
summary(logp.diff.mean)

logp.diff.perm <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		-log10(case$p) - sample(-log10(x$p))
	}
)));
rownames(logp.diff.perm) <- genes;

hist(logp.diff.perm, breaks=100)

logp.diff.perm.mean <- rowMeans(logp.diff.perm, na.rm=TRUE);
logp.diff.perm.sd <- apply(logp.diff.perm, 1, sd, na.rm=TRUE);

logp.diff.z <- (logp.diff.mean - logp.diff.perm.mean) / logp.diff.perm.sd;
logp.diff.p <- 1 - pnorm(logp.diff.z);

hist(logp.diff.z, breaks=100)

hist(10^logp.boot, breaks=100)

logp.boot <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		-log10(x$p)
	}
)));
rownames(logp.boot) <- genes;

hist(10^-logp.boot, breaks=100)

p.boot <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$p
	}
)));
rownames(p.boot) <- genes;
hist(p.boot[p.boot < 1], breaks=100)

logp.boot.d <- melt(logp.boot);

colnames(logp.boot.d) <- c("gene", "type", "statistic");
logp.boot.d$type <- "bootstrap";

logp.case.d <- data.frame(
	gene = case$gene,
	type = "realization",
	statistic = -log10(case$p)
);
logp.case.d.sel <- logp.case.d[sig.idx, ];

logp.d <- rbind(logp.case.d, logp.boot.d);
logp.d.sel <- filter(logp.d, gene %in% sig.genes);
logp.d.sel$gene <- factor(logp.d.sel$gene, levels=rev(genes));

qdraw(
	ggplot(logp.d.sel, aes(x=gene, y=statistic)) +
		geom_jitter(colour="grey60") + theme_bw() + coord_flip() +
		geom_point(data = filter(logp.d.sel, type == "realization"), colour="firebrick", size=2) +
		ylab("MutSig2CV -log(p)")
	,
	width = 6,
	file = insert(out.fname, tag=c("bootstrap", "logp", "jitter"), ext="pdf")
);

qdraw(
	ggplot(logp.d.sel, aes(x=statistic)) +
		facet_grid(gene ~ ., scale="free_y") +
		geom_density(fill="grey60", bw=0.5) + theme_bw() +
		geom_vline(aes(xintercept=statistic), data = filter(logp.d.sel, type == "realization"), colour="firebrick") +
		xlab("MutSig2CV -log(p)")
	,
	height = 15,
	file = insert(out.fname, tag=c("bootstrap", "logp", "density"), ext="pdf")
);


p.le <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$p <= case$p
	}
)));
rownames(p.le) <- genes;

p.le.p <- rowMeans(p.le, na.rm=TRUE);
sort(pmax(p.le.p[sig.idx], 1/B))
sort(p.adjust(pmax(p.le.p[sig.idx], 1/B), "BH"))

hist(p.le.p)



rank.le <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$rank <= case$rank
	}
)));
rownames(rank.le) <- genes;

boot.rank.le.p <- rowMeans(rank.le, na.rm=TRUE);

hist(boot.rank.le.p)

sort(p.adjust(pmax(boot.rank.le.p[sig.idx], 1/B), "BH"))

rank.ge <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$rank >= case$rank
	}
)));
rownames(rank.ge) <- genes;

boot.rank.ge.p <- rowMeans(rank.ge, na.rm=TRUE);

hist(boot.rank.ge.p)

sort(p.adjust(pmax(boot.rank.ge.p[control$q < 0.1], 1/B), "BH"))






# illustrate difference due to bootstrapping

compare_plot <- function(case, control) {
	#show.label <- pmin(control$q, case$q) < 0.25;
	#show.label <- pmin(control$rank, case$rank) < 20;
	show.label <- case$q < 0.25 | control$q < 0.05;

	max.rank <- max(case$rank, control$rank, na.rm=TRUE);
	probs <- (max.rank - pmin(case$rank, control$rank)) / max.rank;
	show.point <- show.label | rbinom(G, 1, probs^100);

	d <- data.frame(
		x = -log10(control$rank),
		y = -log10(case$rank),
		q = -log10(pmin(control$q, case$q)),
		significant = case$q < 0.25,
		gene = ifelse(show.label, genes, NA),
		show = show.point
	);

	ggplot(d, aes(x=x, y=y, label=gene)) +
		stat_density2d(aes(fill=..density..^0.2), geom="raster", contour=FALSE) +
		scale_fill_gradient(low="white", high="royalblue3") +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=x, y=y, size=q), data=d[d$show, ]) +
		scale_size(range=c(0.1, 2)) +
		geom_label_repel(aes(colour=significant), fontface = "italic") +
		scale_colour_manual(values=c("dodgerblue2", "firebrick")) +
		labs(x = "-log(rank) in TCGA LUAD", y = "-log(rank) in BM LUAD", fill="density", size="-log(q)", colour="BM q < 0.25") +
		theme_bw()
}

qdraw(compare_plot(boot[[1]], control), height=7, width=7.5);
