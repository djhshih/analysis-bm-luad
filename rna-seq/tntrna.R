library(io);
library(tntrna);
library(dplyr);

x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");
pheno <- qread("../annot/sample-info_wes_stage2.tsv");

out.fname <- filename("luad-bm_tntrna");

mapping <- match(x$pheno$tumor_id, pheno$fh_sample_id);

# visually check that the histology annotations are consistent
cbind(
	as.character(x$pheno$histology_group),
	as.character(pheno$primary_histotype[mapping])
)

cbind(
	as.character(x$pheno$tumor_id),
	as.character(pheno$fh_sample_id[mapping])
)

# update sample names with standardized annotation
colnames(x$expr) <- pheno$sample_id[mapping];

x$pheno <- mutate(x$pheno,
	sample_id = pheno$sample_id[mapping],
	sample_type = pheno$sample_type[mapping],
	primary_histotype = pheno$primary_histotype[mapping]
);

hist(x$expr, breaks=100)

# a truncated proportion of the expression distribution is approximatly normal
# there *is* a better way to do this with actually truncated normal models
z <- x$expr[x$expr > 1 & x$expr < 5.5];
hist(z, breaks=100, freq=FALSE);
curve(dnorm(x, mean(z), 1.3*sd(z))*1.15, add=TRUE);

mu <- mean(z);
sigma <- 1.35 * sd(z);
hist(x$expr, breaks=100, freq=FALSE, col="skyblue", border=0);
curve(dnorm(x, mu, sigma), add=TRUE, col="royalblue", lwd=2);


####

# use all Non-small cell lung cancers; sample size is too small when analyzing
# Lung adenocarcinoma alone
met.idx <- with(x$pheno,
	#primary_histotype == "Lung adenocarcinoma" &
	sample_type == "Brain metastasis"
);

expr.met <- x$expr[, met.idx]
purity.met <- x$pheno$purity[met.idx];

# use all Non-small cell lung cancers; sample size is too small when analyzing
# Lung adenocarcinoma alone
prim.idx <- with(x$pheno,
	#primary_histotype == "Lung adenocarcinoma" &
	sample_type == "Primary"
);

expr.prim <- x$expr[, prim.idx];
purity.prim <- x$pheno$purity[prim.idx];

fit_gene <- function(y, purity) {
	data <- list(
		N = length(y),
		y = y,
		p = purity,
		nu = 1,
		mu = 5,
		tau = 2,
		kappa = 2,
		omega = 0.3,
		psi = 1
	);

	tntrna_gem_normexp(data, init = list(theta = c(3, 3), sigma=2))
}

plot_gene <- function(y, purity, fit) {
	plot(purity, y, xlim=c(0, 1))
	abline(a=fit$parameters$theta[1], b=diff(fit$parameters$theta), col="green4")
}

fits.prim <- apply(expr.prim, 1, fit_gene, purity = purity.prim);
fits.met <- apply(expr.met, 1, fit_gene, purity = purity.met);
fits.all <- apply(x$expr, 1, fit_gene, purity = x$pheno$purity);

qwrite(fits.prim, insert(out.fname, tag=c("fits", "prim"), ext="rds"));
qwrite(fits.met, insert(out.fname, tag=c("fits", "met"), ext="rds"));

# k = 1: non-tumour
# k = 2: tumour
# FIXME these checks would not be needed if the fitting algorithm
#       were not so brittle...
get_tntexpr <- function(x, k) {
	z <- x$parameters$theta[k];
	if (z > 10 && !x$converged) {
		# something had probably gone wrong during model fitting
		# replace with NA
		z <- NA;
	} else if (z > 20) {
		# such values are simply too big to be plausible, even if the algorithm
		# has apparently converged
		z <- NA;
	}
	z
}

# visually inspect that the estimated expression levels are reasonable

texpr.prim <- unlist(lapply(fits.prim, get_tntexpr, k = 2));
hist(texpr.prim, breaks=100)

ntexpr.prim <- unlist(lapply(fits.prim, get_tntexpr, k = 1));
hist(ntexpr.prim, breaks=100)

texpr.met <- unlist(lapply(fits.met, get_tntexpr, k = 2));
hist(texpr.met, breaks=100)
hist(texpr.met[texpr.met < 10], breaks=100)

ntexpr.met <- unlist(lapply(fits.met, get_tntexpr, k = 1));
hist(ntexpr.met, breaks=100)

texpr.all <- unlist(lapply(fits.all, get_tntexpr, k = 2));
hist(texpr.all, breaks=100)

ntexpr.all <- unlist(lapply(fits.all, get_tntexpr, k = 1));
hist(ntexpr.all, breaks=100)

smoothScatter(ntexpr.met, texpr.met)


sigma.prim <- unlist(lapply(fits.prim, function(x) x$parameters$sigma));

# replace missing values with the median level
missing.prim.idx <- is.na(texpr.prim);
texpr.prim[missing.prim.idx] <- apply(expr.prim[missing.prim.idx, ], 1, median)
missing.met.idx <- is.na(texpr.met);
texpr.met[missing.met.idx] <- apply(expr.met[missing.met.idx, ], 1, median)


qwrite(texpr.prim, insert(out.fname, tag=c("texpr", "prim"), ext="rds"));
qwrite(texpr.met, insert(out.fname, tag=c("texpr", "met"), ext="rds"));
qwrite(texpr.all, insert(out.fname, tag=c("texpr", "all"), ext="rds"));


# interesting observations...

purity_plot(x$expr["MYC", ], as.character(x$pheno$met_type), x$pheno$purity);
plot_gene(expr.prim["MYC", ], purity.prim, fits.prim[["MYC"]])
plot_gene(expr.met["MYC", ], purity.met, fits.met[["MYC"]])

# MYCN is not in x$expr matrix

purity_plot(x$expr["MYCL", ], as.character(x$pheno$met_type), x$pheno$purity);

purity_plot(x$expr["EGFR", ], as.character(x$pheno$met_type), x$pheno$purity);
# Ah... the negative correlation is likely due to sampling?!
# Check this correlation in the TCAG LUAD cohort!
# EGFR mutants might interact with non-EGFR mutants...?
# GBMs with EGFR activating alterations have limited mutant EGFR expression
# "The distribution of EGFR mutations around the ‘catalytic kinase domain’ is
# distinct in lung adenocarcinoma and contrasts the mutations in glioblastomas
# that are located in the extracellular portion of EGFR."  (Siegelin, 2014)
plot_gene(expr.prim["EGFR", ], purity.prim, fits.prim[["EGFR"]])
plot_gene(expr.met["EGFR", ], purity.met, fits.met[["EGFR"]])

purity_plot(x$expr["KRAS", ], as.character(x$pheno$met_type), x$pheno$purity);

purity_plot(x$expr["KEAP1", ], as.character(x$pheno$met_type), x$pheno$purity);
fits.prim["EGFR"]
fits.prim["KRAS"]
fits.prim["KEAP1"]
fits.prim["TP53"]
fits.prim["STK11"]
fits.prim["NF1"]
fits.prim["BRAF"]
fits.prim["CDKN2A"]
fits.prim["PIK3CA"]
# Not even CDKN2A is found in the matrix... RNA degradation or transient expression??!
x$expr[grep("CDKN", rownames(x$expr)), ]

purity_plot(x$expr["COL6A2", ], as.character(x$pheno$met_type), x$pheno$purity);
plot_gene(x$expr["COL6A2", ], x$pheno$purity, fits.all[["COL6A2"]])
plot_gene(x$expr["COL5A2", ], x$pheno$purity, fits.all[["COL6A2"]])

purity_plot(x$expr["COL5A2", ], as.character(x$pheno$met_type), x$pheno$purity);
plot_gene(x$expr["COL5A2", ], x$pheno$purity, fits.all[["COL5A2"]])
plot_gene(x$expr["COL5A2", ], x$pheno$purity, fits.all[["COL5A2"]])

# LAMA1 not in matrix
purity_plot(x$expr["LAMA2", ], as.character(x$pheno$met_type), x$pheno$purity);
purity_plot(x$expr["FN1", ], as.character(x$pheno$met_type), x$pheno$purity);
purity_plot(x$expr["FBN2", ], as.character(x$pheno$met_type), x$pheno$purity);
purity_plot(x$expr["SERPINE2", ], as.character(x$pheno$met_type), x$pheno$purity);
purity_plot(x$expr["PLOD2", ], as.character(x$pheno$met_type), x$pheno$purity);

purity_plot(x$expr["NF1", ], as.character(x$pheno$met_type), x$pheno$purity);
purity_plot(x$expr["SPTAN1", ], as.character(x$pheno$met_type), x$pheno$purity);


# Need to check the TCGA LUAD cohort...

purity_plot(x$expr["TP53", ], as.character(x$pheno$met_type), x$pheno$purity);
plot_gene(expr.prim["TP53", ], purity.prim, fits.prim[["TP53"]])

# too few data points for regression!
# sampling error too high, and shrinkage prior dominates

expressed.cut <- 1;
# use which to remove NA's
expressed.prim <- names(texpr.prim)[which(texpr.prim > expressed.cut)];
expressed.met <- names(texpr.met)[which(texpr.met > expressed.cut)];
expressed.all <- names(texpr.all)[which(texpr.all > expressed.cut)];

genes <- c("EGFR", "KRAS", "KEAP1", "TP53", "STK11", "NF1", "BRAF", "CDKN2A", "PIK3CA");

genes %in% expressed.prim
genes %in% expressed.met
genes %in% expressed.all

qwrite(expressed.prim, insert(out.fname, tag=c("genes", "expressed", "prim"), ext="vtr"));
qwrite(expressed.met, insert(out.fname, tag=c("genes", "expressed", "met"), ext="vtr"));
qwrite(expressed.all, insert(out.fname, tag=c("genes", "expressed", "all"), ext="vtr"));
