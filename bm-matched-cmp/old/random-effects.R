library(io);
library(reshape2);
library(tidyr);
library(rstan);
library(dplyr);

# Analyze matched BM and primary with candidate CNA drivers
# using a random effects model (hierarchical Bayesian binomial model)

out.fname <- filename("pkb-luad", tag="bm-matched-cmp");
pdf.fname <- insert(out.fname, ext="pdf");

compart <- qread("cna-compart-specificity.tsv");
genes <- as.character(unique(compart$gene));

# TODO do transformation automatically instead of hard-coding values
compart.data <- list(
	J = 4,
	x = c(4, 2, 4, 3),
	n = c(5, 4, 4, 4)
);

fit <- stan("hier-binom.stan", data = compart.data, chains=4);
print(fit)

# @param fit  \code{stanfit} object
get_phi_df <- function(fit) {
	theta <- rstan::extract(fit, "theta")[[1]];
	phi <- log(theta / (1 - theta));
	phi.q <- apply(phi, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)));
	colnames(phi.q) <- genes;
	rownames(phi.q) <- c("ymin", "y", "ymax");
	phi.m <- melt(phi.q, varnames=c("variable", "gene"));
	phi.df <- spread(phi.m, "variable", "value") %>%
		mutate(significant = (ymin > 0 & ymax > 0) | (ymin < 0 & ymax < 0));

	phi.df
}

phi.df <- get_phi_df(fit);

qdraw(
	ggplot(phi.df, aes(x=gene, ymin=ymin, y=y, ymax=ymax, colour=significant)) +
		geom_hline(yintercept=0, colour="grey70", lwd=1) +
		geom_errorbar(width=0.5) + geom_point() + 
		scale_colour_manual(values=c("black", "firebrick")) +
		theme_bw() + coord_flip() +
		xlab("") + ylab("log odds ratio of BM vs. Primary exclusivity") +
		theme(legend.position="none")
	,
	width = 4, height = 2,
	file = insert(pdf.fname, c("bm-exclusivity"))
);

####

late <- qread("cna-bm-timing.tsv");

# TODO
early.data <- list(
	J = 4,
	x = c(3, 3, 5, 13),
	n = c(7, 5, 9, 16)
);

fit2 <- stan("hier-binom.stan", data = early.data, chains=4);
print(fit2);

phi.df <- get_phi_df(fit2);

qdraw(
	ggplot(phi.df, aes(x=gene, ymin=ymin, y=y, ymax=ymax, colour=significant)) +
		geom_hline(yintercept=0, colour="grey70", lwd=1) +
		geom_errorbar(width=0.5) + geom_point() + 
		scale_colour_manual(values=c("black", "firebrick")) +
		theme_bw() + coord_flip() +
		xlab("") + ylab("log odds ratio of early vs. late BM mutation") +
		theme(legend.position="none")
	,
	width = 4, height = 2,
	file = insert(pdf.fname, c("bm-early"))
);

