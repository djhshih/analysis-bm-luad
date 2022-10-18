library(io);
library(ggplot2);
library(ggrepel);
library(reshape2);
library(dplyr);

#case.fname <- "mutsig2cv-sig-genes_pkb-luad_brain-mets_black-f_coding.tsv";
#control.fname <- "mutsig2cv-sig-genes_tcga-luad_pass_black-f_coding.tsv";

case.fname <- "mutsig2cv-sig-genes_pkb-luad_brain-mets_black-f_coding_orxnxf.tsv";
control.fname <- "mutsig2cv-sig-genes_tcga-luad_pass_black-f_coding_orxnxf.tsv";

#boot.fname <- "mutsig2cv_bootstrap_p.rds";
#boot.type <- "boot";

#boot.fname <- "mutsig2cv_subsample_p.rds";
#boot.type <- "subsample";

#boot.fname <- "mutsig2cv_subsamplew_p.rds";
#boot.type <- "subsamplew";

boot.fname <- "mutsig2cv_subsamplew-orxnxf_p.rds";
boot.type <- "subsamplew-orxnxf";



out.fname <- filename("pkb-luad");

case <- qread(case.fname, quote="");
control <- qread(control.fname, quote="");

boot.orig <- qread(boot.fname);


genes <- intersect(case$gene, control$gene);
G <- length(genes);
B <- ncol(boot.orig$p);

qqplot.width <- 5;
qqplot.height <- 5;

compare.fdr <- 0.025;

# rearrange boot to match gene order in case, control, and bootstrap samples

case <- case[match(genes, case$gene), ];
control <- control[match(genes, control$gene), ];

sig.idx <- which(case$q < 0.25);
control.sig.idx <- which(control$q < 0.05);
sum(case$q > 0.25 & control$q < 0.05)
sig.genes <- genes[sig.idx];

boot.idx <- match(genes, boot.orig$genes);

boot <- boot.orig;
boot$genes <- as.character(boot$genes[boot.idx]);
boot$p <- boot$p[boot.idx, ];

ranks <- apply(boot$p, 2,
	function(p) {
		match(1:length(p), order(p))
	}
);
rownames(ranks) <- genes;

# store realized samples of case and control
case.r <- case;
control.r <- control;

p.boot <- boot$p;
rownames(p.boot) <- genes;

# MutSig2CV p-values are strange
# because many genes cannot be evaluated
# and their p-values are assigned values of 1
# or close to 1 (final p-value is a combination of
# 3 p-values)
hist(p.boot[p.boot < 1], breaks=100)


logp.boot.d <- melt(-log10(p.boot));

colnames(logp.boot.d) <- c("gene", "boot", "statistic");
logp.boot.d$type <- boot.type;

logp.case.d <- data.frame(
	gene = case$gene,
	type = "realization",
	statistic = -log10(case$p)
);

logp.case.d.sel <- logp.case.d[sig.idx, ];
logp.boot.d.sel <- dplyr::filter(logp.boot.d, gene %in% sig.genes);

logp.boot.d.sel$gene <- factor(logp.boot.d.sel$gene, levels=rev(genes));
qdraw(
	ggplot(filter(logp.boot.d.sel, boot < 50), aes(x=gene, y=statistic)) +
		geom_jitter(colour="grey60", alpha=0.5, width=0.2) + theme_bw() + coord_flip() +
		geom_point(data=logp.case.d.sel, colour="firebrick", size=2) +
		ylab("MutSig2CV -log(p)") + xlab("") +
		theme(
			axis.text.y = element_text(face="bold.italic", colour="black")
		)
	,
	width = 4,
	height = 5,
	file = insert(out.fname, tag=c(boot.type, "logp", "jitter"), ext="pdf")
);

qdraw(
	ggplot(logp.boot.d.sel, aes(x=gene, y=statistic)) +
		geom_violin(fill="grey60", alpha=0.5, adjust=3) + theme_bw() + coord_flip() +
		geom_point(data=logp.case.d.sel, colour="firebrick", size=2) +
		ylab("MutSig2CV -log(p)") + xlab("") +
		theme(
			axis.text.y = element_text(face="bold.italic", colour="black")
		)
	,
	width = 4,
	height = 5,
	file = insert(out.fname, tag=c(boot.type, "logp", "violin"), ext="pdf")
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
	file = insert(out.fname, tag=c(boot.type, "logp", "density"), ext="pdf")
);


# calculate bootstrap p-values for more significant MutSig2Cv p-values

hist(case$p, breaks=1000)
hist(boot$p, breaks=1000)

mean(case$p == 1)
mean(boot$p == 1, na.rm=TRUE)

# when MutSigCV p-value is 1, it means that significance cannot be assessed
# therefore, draw a random number from the uniform distribution as
null.p.idx <- case$p > 0.99;
n.null.p <- sum(null.p.idx);

set.seed(1234);

p.case <- matrix(rep(case$p, B), nrow=length(case$p));
p.case[null.p.idx] <- runif(sum(null.p.idx));

# same issue with subsample p-values
boot.null.p.idx <- is.na(p.boot) | p.boot > 0.99;
p.boot[boot.null.p.idx] <- runif(sum(boot.null.p.idx));


hist(p.case, breaks=1000)
hist(p.boot, breaks=1000)

mean(p.case == 1)
mean(p.boot == 1)


#p.le <- p.boot <= case$p;
p.le <- p.boot <= p.case;

p.le.p <- rowMeans(p.le, na.rm=TRUE);
p.le.p.orig <- p.le.p;
sort(pmax(p.le.p[sig.idx], 1/B))
p.le.sel.bh <- sort(p.adjust(pmax(p.le.p[sig.idx], 1/B), "BH"))
p.le.sel.bh

sort(p.le.p[sig.idx])
p.le.sel.bh <- sort(p.adjust(p.le.p[sig.idx], "BH"))
p.le.sel.bh

p.boot[p.le.p == 1, ]
case$p[p.le.p == 1]

summary(case$p[p.le.p == 1])


# Calculate p-values using GPD

#library(fExtremes);
#x <- gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
#fit <- gpdFit(x, u = min(x), type = "pwm")

hist(-log10(p.boot))

lp.boot <- -log10(p.boot);
lp.case <- -log10(p.case);

lp.ge <- lp.boot >= lp.case;

hist(lp.boot)
hist(lp.case)


library(gPdtest)

gfit <- gpd.fit(lp.case, method="combined");



pgp <- function(x, shape, scale, lower.tail = TRUE, log.p = FALSE) {
	if (shape == 0) {
		pexp(x, rate = 1 / scale, lower.tail, log.p)
	} else {
		v <- (1 + (shape / scale) * x);
		if (shape < 0) {
			# x must be < -scale/shape
			# if not, v will be negative (invalid)
			# therefore, 
			v <- pmax(v, .Machine$double.eps);
		}
		if (lower.tail) {
			p <- 1 - v^(-1/shape);
		} else {
			p <- v^(-1/shape);
		}
		if (log.p) {
			log(p)	
		} else {
			p
		}
	}
}

z <- rgp(1000, shape = gfit[1], scale = gfit[2]);
hist(z)
plot(ecdf(z))
lines(curve(pgp(x, gfit[1], gfit[2]), from=0, to=4, add=TRUE), col="red")



plot(curve(pexp(x, 1/gfit[2]), from=0, to=10), type="l")
plot(curve(dexp(x, 1/gfit[2]), from=0, to=10), type="l")


p_gpd <- function(x0, x, u=NULL, top=250, method="combined") {
	y <- sort(x, decreasing=TRUE);
	if (is.null(u)) {
		u <- (y[top] + y[top+1]) / 2;
		z <- y[1:top] - u;
	} else {
		z <- y[y > u] - u;
		top <- length(z);
	}
	
	params <- gPdtest::gpd.fit(z, method=method);
	print(params)
	if (is.finite(params[1]) && params[1] < -0.5) {
		params[1] <- -0.5;
	}
	N <- length(x);

	(top / N) * (1 - pgp(x0 - u, shape=params[1], scale=params[2]))
}

library(MASS)
params <- fitdistr(z, "exponential");

p_exp <- function(x0, x, u=NULL, top=200, alternative=c("greater", "less")) {
	n <- length(x);
	alternative <- match.arg(alternative);
	if (alternative == "greater") {
		y <- sort(x, decreasing=TRUE)
		if (is.null(u)) {
			u <- (y[top] + y[top+1]) / 2;
			z <- y[1:top] - u;
		} else {
			z <- y[y > u] - u;
			top <- length(z);
		}

		rate <- coef(fitdistr(z, "exponential"));
		(top/n) * pexp(x0 - u, rate=rate, lower.tail=FALSE)
	} else {
		y <- sort(x, decreasing=FALSE)
		if (is.null(u)) {
			u <- (y[top] + y[top+1]) / 2;
			z <- u - y[1:top];
		} else {
			z <- u - y[y < u];
			top <- length(z);
		}


		rate <- coef(fitdistr(z, "exponential"));

		#hist(z, freq=FALSE);
		#lines(curve(dexp(x, rate=rate), add=TRUE), col="red");

		(top/n) * pexp(u - x0, rate=rate, lower.tail=FALSE)
	}
}




#p_gpd(lp.case[1,1], lp.boot[1,])
#p.le.p[1]

#p.le.p.orig <- p.le.p;

case[case$gene == "OR8H2", ]
lp.case[case$gene == "OR8H2", 1]
lp.boot[case$gene == "OR8H2", ]

# approximate the tail p-value
approx.p.idx <- which(p.le.p.orig < 30 / B);
p.le.p[approx.p.idx] <- unlist(lapply(approx.p.idx,
	function(i) {
		#p_gpd(lp.case[i, 1], lp.boot[i, ])
		p_exp(lp.case[i, 1], lp.boot[i, ])
	}
));

r <- 2;

p_gpd(lp.case[approx.p.idx[r], 1], lp.boot[approx.p.idx[r], ])
lp.case[approx.p.idx[r]]
p.le.p.orig[approx.p.idx[r]]
p.le.p[approx.p.idx[r]]

z <- lp.boot[approx.p.idx[r], ];

fit.amle <- gPdtest::gpd.fit(z, method="amle");
fit.amle
fit.combined <- gPdtest::gpd.fit(z, method="combined");
fit.combined
fit.exp <- fitdistr(z, "exponential");
fit.exp

(300 / B) * (1 - pgp(lp.case[approx.p.idx[1]], shape=fit.combined[1], scale=fit.combined[2]))
pgp(lp.case[approx.p.idx[1]], shape=fit.combined[1], scale=fit.combined[2], lower.tail=FALSE)

(300 / B) * (1 - pgp(lp.case[approx.p.idx[1]], shape=fit.amle[1], scale=fit.amle[2]))

# parameter estimation by the combined method is more accurate
hist(z, breaks=50)
plot(ecdf(z))
lines(curve(pgp(x, fit.amle[1], fit.amle[2]), from=0, to=4, add=TRUE), col="blue")
lines(curve(pgp(x, fit.combined[1], fit.combined[2]), from=0, to=4, add=TRUE), col="red")
lines(curve(pexp(x, coef(fit.exp)), from=0, to=4, add=TRUE), col="green")

sum(is.na(p.le.p))

plot(p.le.p.orig, p.le.p)
abline(a=0, b=1, col="red")
plot(p.le.p.orig, p.le.p, xlim=c(0, 0.1))
abline(a=0, b=1, col="red")

hist(p.le.p.orig, breaks=500)
hist(p.le.p, breaks=500)
hist(p.le.p, breaks=500, xlim=c(0, 0.3))

inflation_factor <- function(p, quantiles=0.5) {
	chisq <- qchisq(1 - p, 1);
	lambda <- quantile(chisq, quantiles) / qchisq(quantiles, 1);
	lambda
}



p_qqplot <- function(p, label.idx=NULL, fdr=0.1, title="Quantile-quantile plot of p-values", spartan=FALSE, colour="black", force=2, subsample=1.0) {
	n <- length(p);
	idx <- order(p, decreasing=FALSE);
	o <- -log10(p[idx]);
	e <- -log10((1:n)/n);
	# you could use base graphics
	#plot(e,o,pch=19,cex=0.25, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
	#lines(e,e,col="red")
	#You'll need ggplot2 installed to do the rest

	#p.le.p.lambda <- inflation_factor(p.le.p[p.le.p < 1]);
	#p.le.p.lambdas <- inflation_factor(p.le.p[p.le.p < 1], seq(0, 1, 0.1));
	lambda <- inflation_factor(p);

	if (subsample < 1.0) {
		h <- hist(o, breaks=min(floor(n/10), 500), plot=FALSE);
		prob <- h$density / sum(h$density);
		weight <- 1 / prob;
		binned <- as.integer(cut(o, h$breaks));
		binned[is.na(binned)] <- 1;
		ind <- sample.int(n, size=ceiling(n * subsample), prob=weight[binned], replace=FALSE);
		keep.idx <- 1:n %in% ind;
	} else {
		keep.idx <- rep(TRUE, n);
	}

	omax <- max(e, o) * 1.1;
	d <- data.frame(e, o, label=names(o), significant=FALSE);

	if (is.null(label.idx)) {
		if (n > 10) {
			d$label <- NA;
		}
	} else {
		labels.sel <- names(p)[label.idx];
		d$label[! d$label %in% labels.sel ] <- NA;

		p.sel <- p[label.idx];
		q <- p.adjust(p.sel, "BH");
		sig.idx <- q[idx] < fdr;
		d$significant[sig.idx] <- TRUE;
	}

	d.sel <- d[keep.idx | !is.na(d$label), ];

	g <- ggplot(d.sel, aes(x=e, y=o, label=label)) +
		#geom_smooth(color="grey60") +
		geom_point(size=1) +
		#xlim(0, omax) + ylim(0, omax) + 
		geom_abline(intercept=0,slope=1, col="red") +
		geom_label_repel(aes(fill=significant), colour=colour, size=4, fontface="italic", force=force) +
		#geom_text_repel(aes(colour=significant)) +
		labs(title=title) +
		xlab(expression(Expected~~-log[10](italic(p)))) +
		ylab(expression(Observed~~-log[10](italic(p)))) +
		annotate("text", x=0, y=max(o), hjust=-0.1, vjust=0.1, label=sprintf("lambda = %s", format(lambda, digits=2)), size=5) +
		theme_bw() +
		theme(
			axis.text=element_text(size=12),
			axis.title=element_text(size=14),
			legend.text=element_text(size=12),
			legend.title=element_text(size=12)
		)
	
	if (spartan) {
		g <- g + theme(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
	}
	g
}


qdraw(
	qplot(p.le.p, geom="histogram", bins=50) + theme_bw()
	,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-le", "hist"), ext="pdf")
);

qdraw(
	qplot(p.le.p[p.le.p < 1], geom="histogram", bins=50) + theme_bw()
	,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-le", "lt-1", "hist"), ext="pdf")
);

p.le.p.q <- p.adjust(p.le.p[sig.idx], "BH");
sort(p.le.p.q)

qwrite(names(p.le.p.q[p.le.p.q < compare.fdr]), insert(out.fname, c("mutsig2cv-compare", "sig"), ext="vtr"))

qdraw(
	p_qqplot(p.le.p, sig.idx, fdr=compare.fdr, colour="firebrick", subsample=0.01,
		title="BM-LUAD MutSig2CV q < 0.25") + 
		scale_fill_manual(values=c("white", "#FFD6E1")) +
		labs(colour=expression(italic(q) < compare.fdr)) +
		xlim(0, 5.5) + ylim(0, 5.5) +
		theme(legend.position="none")
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-le-approx", "qq"), ext="pdf")
);

qdraw(
	p_qqplot(p.le.p.orig, sig.idx, fdr=compare.fdr, colour="firebrick", subsample=0.01) + 
		scale_fill_manual(values=c("white", "pink")) +
		labs(colour=expression(italic(q) < compare.fdr)) +
		xlim(0, 5.5) + ylim(0, 5.5) +
		theme(legend.position="none")
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-le", "qq"), ext="pdf")
);

qdraw(
	p_qqplot(p.le.p, subsample=0.02) + theme(legend.position="none") +
		xlim(0, 5.5) + ylim(0, 5.5)
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-le-approx", "unlabeled", "qq"), ext="pdf")
);

qdraw(
	p_qqplot(p.le.p.orig) +  xlim(0, 5.5) + ylim(0, 5.5) + theme(legend.position="none")
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-le", "unlabeled", "qq"), ext="pdf")
);

#stop()

inflation_factor <- function(p, quantiles=0.5) {
	chisq <- qchisq(1 - p, 1);
	lambda <- quantile(chisq, quantiles) / qchisq(quantiles, 1);
	lambda
}

# START OBSOLETE ################

#p.le.p.lambda <- inflation_factor(p.le.p[p.le.p < 1]);
#p.le.p.lambdas <- inflation_factor(p.le.p[p.le.p < 1], seq(0, 1, 0.1));
p.le.p.lambda <- inflation_factor(p.le.p);
p.le.p.lambdas <- inflation_factor(p.le.p, seq(0, 1, 0.1));

#p.le.p.deflated <- pmin(p.le.p * p.le.p.lambda, 1);
p.le.p.deflated <- pmin(p.le.p * p.le.p.lambdas[3], 1);

qplot(p.le.p.deflated[p.le.p.deflated < 1], geom="histogram", bins=11)

p_qqplot(p.le.p.deflated[p.le.p.deflated < 1])

# inflation factor is *not* constant across
# the full range of p-values... 
# so a constant correction is insufficient!

# END OBSOLETE ################


# calculate bootstrap p-values for less significant MutSig2Cv p-values

p.ge <- p.boot >= p.case;

p.ge.p <- rowMeans(p.ge, na.rm=TRUE);
p.ge.p.orig <- p.ge.p;
hist(p.ge.p, breaks=100);
sort(pmax(p.ge.p[control.sig.idx], 1/B))
p.ge.sel.bh <- sort(p.adjust(pmax(p.ge.p[control.sig.idx], 1/B), "BH"));

approx.p.ge.idx <- which(p.ge.p.orig < 30 / B);
p.ge.p[approx.p.ge.idx] <- unlist(lapply(approx.p.ge.idx,
	function(i) {
		#p_gpd(lp.case[i, 1], lp.boot[i, ])
		p_exp(lp.case[i, 1], lp.boot[i, ], top=150, alternative="less")
	}
));

head(sort(p.ge.p))
p.ge.sel.bh <- sort(p.adjust(p.ge.p, "BH"));
head(p.ge.sel.bh)

qdraw(
	qplot(p.ge.p, geom="histogram", bins=100)
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-ge", "hist"), ext="pdf")
);

qdraw(
	qplot(p.ge.p[p.ge.p < 1], geom="histogram", bins=100)
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-ge", "lt-1", "hist"), ext="pdf")
);

qdraw(
	p_qqplot(p.ge.p, setdiff(control.sig.idx, sig.idx), fdr=compare.fdr, colour="dodgerblue4",
		force=3, title="BM-LUAD MutSig2CV q >= 0.25", subsample=0.01) + 
		scale_fill_manual(values=c("white", "lightblue")) +
		labs(colour=expression(italic(q) < compare.fdr)) +
		xlim(0, 2.5) + ylim(0, 2.5) +
		theme(legend.position = "none")
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-ge-approx", "qq"), ext="pdf")
);

qdraw(
	p_qqplot(p.ge.p, subsample=0.02) + theme(legend.position="none") +
		xlim(0, 2.5) + ylim(0, 2.5) +
		theme(legend.position = "none")
	,
	width = qqplot.width, height = qqplot.height,
	file = insert(out.fname, c("mutsigcv2-p", boot.type, "p-ge", "unlabeled", "qq"), ext="pdf")
);


#p.le.boot <- lapply(
#	1:B,
#	function(b) {
#		rowMeans(as.matrix(as.data.frame(lapply(
#			setdiff(1:B, b),
#			function(b2) {
#				#boot[[b2]]$p <= boot[[b]]$p
#				boot$p[, b2] <= boot$p[, b]
#			}
#		))), na.rm=TRUE)
#	}
#);

#p.le.p.sel <- pmax(sort(p.le.p[sig.idx]), 1/B);
# mathematically, 1/B should be 1/(B-1), because the p-values of the bootstrap
# samples are based on B - 1 other samples
# but 1/B is more conservative
# ideally, we should add an additional boostrap sample for estimate comparison
# p-values in the bootstrap
#p.le.boot.sel <- lapply(p.le.boot, function(x) pmax(sort(x[sig.idx]), 1/B));

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

#p.le.sel.q <- boot_fdr(p.le.p.sel, p.le.boot.sel);

#cor(p.le.sel.q, p.le.sel.bh)
#plot(p.le.sel.q, p.le.sel.bh)
#abline(a=0, b=1)

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


#p.ge.boot <- lapply(
#	1:B,
#	function(b) {
#		rowMeans(as.matrix(as.data.frame(lapply(
#			setdiff(1:B, b),
#			function(b2) {
#				boot$p[, b2] >= boot$p[, b]
#			}
#		))), na.rm=TRUE)
#	}
#);

#p.ge.p.sel <- pmax(sort(p.ge.p[control.sig.idx]), 1/B);
#p.ge.boot.sel <- lapply(p.ge.boot, function(x) pmax(sort(x[control.sig.idx]), 1/B));

#p.ge.sel.q <- boot_fdr(p.ge.p.sel, p.ge.boot.sel);

#cor(p.ge.sel.bh, p.ge.sel.q)
#plot(p.ge.sel.bh, p.ge.sel.q)


# illustrate difference due to bootstrapping

compare_plot <- function(case, control, genes, cohorts) {
	#show.label <- pmin(control$q, case$q) < 0.25;
	#show.label <- pmin(control$rank, case$rank) < 20;
	show.label <- case$q < 0.25 | control$q < 0.05;

	max.rank <- max(nrow(case), nrow(control));
	case$rank <- match(1:nrow(case), order(case$p));
	control$rank <- match(1:nrow(control), order(control$p));

	probs <- (max.rank - pmin(case$rank, control$rank)) / max.rank;
	show.point <- show.label | rbinom(G, 1, probs^100);

	d <- data.frame(
		x = -log10(control$rank),
		y = -log10(case$rank),
		q = -log10(pmin(control$q, case$q)),
		significant = factor(case$q < 0.25, label=c("q >= 0.25", "q < 0.25")),
		gene = ifelse(show.label, genes, NA),
		show = show.point
	);

	labs <- paste0("-log10(rank) in ", cohorts);

	ggplot(d, aes(x=x, y=y, label=gene)) +
		stat_density2d(aes(fill=..density..^0.2), geom="raster", contour=FALSE) +
		scale_fill_gradient(low="white", high="royalblue2") +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=x, y=y), data=d[d$show, ], size=1) +
		#geom_point(aes(x=x, y=y, size=q), data=d[d$show, ]) +
		#scale_size(range=c(0.1, 2)) +
		geom_label_repel(aes(colour=significant), fontface = "italic", size=4.1) +
		scale_colour_manual(values=c("dodgerblue4", "firebrick")) +
		labs(x = labs[1], y = labs[2], fill="density", size="-log(q)", colour="BM-LUAD") +
		theme_bw() +
		theme(
			axis.text=element_text(size=12),
			axis.title=element_text(size=14),
			legend.text=element_text(size=12),
			legend.title=element_text(size=12)
		)
}

make_mutsig2cv <- function(p) {
	data.frame(
		p = p,
		q = p.adjust(p, "BH")
	)
}

cohorts.re <- c("Control cohort", "Control cohort subsample");
qdraw(
	compare_plot(make_mutsig2cv(boot$p[,1]), control, genes, cohorts=cohorts.re),
	height=6, width=7.2,
	file = insert(out.fname, c("compare-boot", boot.type, "b1", "mutsig2cv-rank"), ext="pdf")
);

qdraw(
	compare_plot(make_mutsig2cv(boot$p[,2]), control, genes, cohorts=cohorts.re),
	height=6, width=7.2,
	file = insert(out.fname, c("compare-boot", boot.type, "b2", "mutsig2cv-rank"), ext="pdf")
);

qdraw(
	compare_plot(make_mutsig2cv(boot$p[,3]), control, genes, cohorts=cohorts.re),
	height=6, width=7.2,
	file = insert(out.fname, c("compare-boot", boot.type, "b3", "mutsig2cv-rank"), ext="pdf")
);

qdraw(
	compare_plot(make_mutsig2cv(boot$p[,4]), control, genes, cohorts=cohorts.re),
	height=6, width=7.2,
	file = insert(out.fname, c("compare-boot", boot.type, "b4", "mutsig2cv-rank"), ext="pdf")
);


####

hist(lp.boot[5,])
p_exp(10, lp.boot[5, ], alternative="less")

hist(lp.boot[1, ])
p.ge.p.orig[1]
p_exp(16, lp.boot[1, ], top=400, alternative="less")

hist(lp.boot[10, ], breaks=30)
p_exp(16, lp.boot[10, ], top=100, alternative="less")

head(approx.p.ge.idx)
hist(lp.boot[2304, ], breaks=30)
p_exp(2, lp.boot[2304, ], top=150, alternative="less")

hist(lp.boot[2797, ], breaks=30)
p_exp(2, lp.boot[2797, ], top=150, alternative="less")

hist(lp.boot[5145, ], breaks=30)
p_exp(2, lp.boot[5145, ], top=150, alternative="less")

