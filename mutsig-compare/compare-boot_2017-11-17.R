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

genes <- intersect(case$gene, control$gene);
G <- length(genes);
B <- length(boot.orig);
sig.idx <- which(case$q < 0.25);
control.sig.idx <- which(control$q < 0.1);
sig.genes <- genes[sig.idx];

# rearrnage boot to match gene order in case, control, and bootstrap samples

case <- case[match(genes, case$gene), ];
control <- control[match(genes, control$gene), ];

boot <- lapply(
	boot.orig,
	function(x) {
		x[match(genes, x$gene), ]
	}
);

# store realized samples of case and control
case.r <- case;
control.r <- control;



p.boot <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$p
	}
)));
rownames(p.boot) <- genes;

hist(p.boot[p.boot < 1], breaks=100)


logp.boot <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		-log10(x$p)
	}
)));
rownames(logp.boot) <- genes;

logp.boot.d <- melt(logp.boot);

colnames(logp.boot.d) <- c("gene", "type", "statistic");
logp.boot.d$type <- "bootstrap";

logp.case.d <- data.frame(
	gene = case$gene,
	type = "realization",
	statistic = -log10(case$p)
);

logp.case.d.sel <- logp.case.d[sig.idx, ];
logp.boot.d.sel <- filter(logp.boot.d, gene %in% sig.genes);

logp.boot.d.sel$gene <- factor(logp.boot.d.sel$gene, levels=rev(genes));
qdraw(
	ggplot(logp.boot.d.sel, aes(x=gene, y=statistic)) +
		geom_jitter(colour="grey60", alpha=0.3, width=0.2) + theme_bw() + coord_flip() +
		geom_point(data=logp.case.d.sel, colour="firebrick", size=2) +
		ylab("MutSig2CV -log(p)")
	,
	width = 4,
	height = 5,
	file = insert(out.fname, tag=c("bootstrap", "logp", "jitter"), ext="pdf")
);

logp.boot.d.sel$gene <- factor(logp.boot.d.sel$gene, levels=genes);
qdraw(
	ggplot(logp.boot.d.sel, aes(x=statistic)) +
		facet_wrap(~gene, scale="free_y") +
		geom_density(fill="grey60", bw=0.5) + theme_bw() +
		geom_vline(aes(xintercept=statistic), data=logp.case.d.sel, colour="firebrick") +
		xlab("MutSig2CV -log(p)")
	,
	width = 6,
	height = 5,
	file = insert(out.fname, tag=c("bootstrap", "logp", "density"), ext="pdf")
);


# calculate bootstrap p-values for more significant MutSig2Cv p-values

p.le <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$p <= case$p
	}
)));
rownames(p.le) <- genes;

p.le.p <- rowMeans(p.le, na.rm=TRUE);
hist(p.le.p, breaks=100);
sort(pmax(p.le.p[sig.idx], 1/B))
p.le.sel.bh <- sort(p.adjust(pmax(p.le.p[sig.idx], 1/B), "BH"))


p_qqplot <-  function(p, title="Quantile-quantile plot of p-values", spartan=F) {
	n <- length(p);
	o <- -log10(sort(p, decreasing=F))
	e <- -log10( 1:n/n )
	# you could use base graphics
	#plot(e,o,pch=19,cex=0.25, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
	#lines(e,e,col="red")
	#You'll need ggplot2 installed to do the rest

	g <- qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + 
		geom_abline(intercept=0,slope=1, col="red") +
		labs(title=title) +
		scale_x_continuous(name=expression(Expected~~-log[10](italic(p)))) +
		scale_y_continuous(name=expression(Observed~~-log[10](italic(p))));
	
	if (spartan) {
		g <- g + theme(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
	}
	g
}

qdraw(
	qplot(p.le.p, geom="histogram", bins=100) + theme_bw(),
	file = insert(out.fname, c("mutsigcv2-p", "bootstrap-p-le", "hist"), ext="pdf")
);
qdraw(
	qplot(p.le.p[p.le.p < 1], geom="histogram", bins=100) + theme_bw(),
	file = insert(out.fname, c("mutsigcv2-p", "bootstrap-p-le", "lt-1", "hist"), ext="pdf")
);
qdraw(
	p_qqplot(p.le.p[p.le.p < 1]) + theme_bw()
	,
	file = insert(out.fname, c("mutsigcv2-p", "bootstrap-p-le", "qq"), ext="pdf")
);


inflation_factor <- function(p, quantiles=0.5) {
	chisq <- qchisq(1 - p, 1);
	lambda <- quantile(chisq, quantiles) / qchisq(quantiles, 1);
	lambda
}

p.le.p.lambda <- inflation_factor(p.le.p[p.le.p < 1]);
p.le.p.lambdas <- inflation_factor(p.le.p[p.le.p < 1], seq(0, 1, 0.1));

#p.le.p.deflated <- pmin(p.le.p * p.le.p.lambda, 1);
p.le.p.deflated <- pmin(p.le.p * p.le.p.lambdas[3], 1);

qplot(p.le.p.deflated[p.le.p.deflated < 1], geom="histogram", bins=11) + theme_bw()

p_qqplot(p.le.p.deflated[p.le.p.deflated < 1]) + theme_bw()

# inflation factor is *not* constant across
# the full range of p-values... 
# so a constant correction is insufficient!

# calculate bootstrap p-values for less significant MutSig2Cv p-values

p.ge <- as.matrix(as.data.frame(lapply(
	boot,
	function(x) {
		x$p >= case$p
	}
)));
rownames(p.ge) <- genes;

p.ge.p <- rowMeans(p.ge, na.rm=TRUE);
hist(p.ge.p, breaks=100);
sort(pmax(p.ge.p[control.sig.idx], 1/B))
p.ge.sel.bh <- sort(p.adjust(pmax(p.ge.p[control.sig.idx], 1/B), "BH"));

qdraw(
	qplot(p.ge.p, geom="histogram", bins=100) + theme_bw(),
	file = insert(out.fname, c("mutsigcv2-p", "bootstrap-p-ge", "hist"), ext="pdf")
);

qdraw(
	qplot(p.ge.p[p.ge.p < 1], geom="histogram", bins=100) + theme_bw(),
	file = insert(out.fname, c("mutsigcv2-p", "bootstrap-p-ge", "lt-1", "hist"), ext="pdf")
);

qdraw(
	p_qqplot(p.ge.p[p.ge.p < 1]) + theme_bw(),
	file = insert(out.fname, c("mutsigcv2-p", "bootstrap-p-ge", "qq"), ext="pdf")
);

p.le.boot <- lapply(
	1:B,
	function(b) {
		rowMeans(as.matrix(as.data.frame(lapply(
			setdiff(1:B, b),
			function(b2) {
				boot[[b2]]$p <= boot[[b]]$p
			}
		))), na.rm=TRUE)
	}
);

p.le.p.sel <- pmax(sort(p.le.p[sig.idx]), 1/B);
# mathematically, 1/B should be 1/(B-1), because the p-values of the bootstrap
# samples are based on B - 1 other samples
# but 1/B is more conservative
# ideally, we should add an additional boostrap sample for estimate comparison
# p-values in the bootstrap
p.le.boot.sel <- lapply(p.le.boot, function(x) pmax(sort(x[sig.idx]), 1/B));

# TODO optimize
#' @param p  numeric vector of sorted p-values
#' @param boot.p  list of vectors of sorted p-values; one vector for each
#'                bootstrap samples
#' @return numeric vector of false discovery rates at each p-value threshold
boot_fdr <- function(p, boot.p) {
	unlist(lapply(
		p,
		function(p.cut) {
			# R = number of rejections
			# V = number of false discoveries
			# estimate E[R] using observed number of rejections
			r <- sum(p <= p.cut);
			# estimate E[V] using number of rejections observed in bootstrap samples
			v <- mean(unlist(lapply(boot.p, function(x) sum(x <= p.cut))));
			# estimate FDR = E[V/R] \approx E[V] / E[R]
			#v / r
			if (v > r) 1 else v / r
		}
	))
}

p.le.sel.q <- boot_fdr(p.le.p.sel, p.le.boot.sel);

cor(p.le.sel.q, p.le.sel.bh)
plot(p.le.sel.q, p.le.sel.bh)
abline(a=0, b=1)

#p.le.p.sorted <- pmax(sort(p.le.p), 1/B);
#p.le.boot.sorted <- lapply(p.le.boot, function(x) pmax(sort(x), 1/B));

#p.le.q <- boot_fdr(p.le.p.sorted, p.le.boot.sorted);
#p.le.q[names(p.le.q) %in% sig.genes]
#p.le.bh <- p.adjust(p.le.p.sorted, "BH");

#plot(p.le.p.sorted, p.le.bh)

#cor(p.le.q, p.le.bh)
#plot(p.le.q, p.le.bh)
#abline(a=0, b=1)

#hist(p.le.bh, breaks=100)
#hist(p.le.q, breaks=100)


p.ge.boot <- lapply(
	1:B,
	function(b) {
		rowMeans(as.matrix(as.data.frame(lapply(
			setdiff(1:B, b),
			function(b2) {
				boot[[b2]]$p >= boot[[b]]$p
			}
		))), na.rm=TRUE)
	}
);

p.ge.p.sel <- pmax(sort(p.ge.p[control.sig.idx]), 1/B);
p.ge.boot.sel <- lapply(p.ge.boot, function(x) pmax(sort(x[control.sig.idx]), 1/B));

p.ge.sel.q <- boot_fdr(p.ge.p.sel, p.ge.boot.sel);

cor(p.ge.sel.bh, p.ge.sel.q)
plot(p.ge.sel.bh, p.ge.sel.q)


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

qdraw(
	qdraw(compare_plot(boot[[1]], control), height=7, width=7.5),
	file = insert(out.fname, c("compare-boot", "b1", "mutsig2cv-rank"), ext="pdf")
);

qdraw(
	qdraw(compare_plot(boot[[2]], control), height=7, width=7.5),
	file = insert(out.fname, c("compare-boot", "b2", "mutsig2cv-rank"), ext="pdf")
);

qdraw(
	qdraw(compare_plot(boot[[3]], control), height=7, width=7.5),
	file = insert(out.fname, c("compare-boot", "b3", "mutsig2cv-rank"), ext="pdf")
);

