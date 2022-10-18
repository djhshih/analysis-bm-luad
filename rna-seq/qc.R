library(io);
library(ggplot2);
library(reshape2);
library(preprocessCore);
library(mmalign);
library(lme4);

rpkm <- qread("luad-bm_rna-seq.tsv", type="mtx");

pheno.all <- qread("../annot/sample-info.tsv");

qc <- qread("Brastianos_BMET_RNAQC.txt", type="tsv");
qc$rnaseq_id <- qc$"Collaborator Sample ID";
qc <- qc[match(colnames(rpkm), qc$rnaseq_id), ];

fresh.frozen.samples <- c("PB0052-MT", "0067-MT", "PB0225-MT", "PB0300-MT", "PB0321-MT");

pheno <- pheno.all[match(colnames(rpkm), pheno.all$rnaseq_id), ];
stopifnot(all(pheno$rnaseq_id == colnames(rpkm)));

pheno$dv200 <- qc$DV200;

pheno$fresh_frozen <- pheno$rnaseq_id %in% fresh.frozen.samples;


out.fname <- filename("luad-bm", ext="rds");


# Independent common functions

rpkm_to_tpm <- function(x) {
	x / sum(x) * 1e6
}

detectable <- function(X, cutoff=0) {
	apply(X, 1, mean) > cutoff
}

prop_zero <- function(X) {
	apply(X, 2, function(x) mean(x == 0))
}


# Set up expression matrices

tpm.m.colnames <- c("ensembl", "sample", "log2_tpm");

tpm <- apply(rpkm, 2, rpkm_to_tpm);

pheno$sum_rpkm <- apply(rpkm, 2, sum);


# Assess baseline density distribution before any normalization or filtering
# Plot only a subsample for efficiency

sample.idx <- sample(1:nrow(tpm), 10000);
tpm.sub <- tpm[sample.idx, ];

tpm.sub.l.m <- melt(log2(tpm.sub+1));
colnames(tpm.sub.l.m) <- tpm.m.colnames;
g <- ggplot(tpm.sub.l.m, aes(x=log2_tpm, fill=sample)) + geom_density(alpha=0.2) + ylim(0, 1);
g


# Identity detectable genes and 
# filter out any genes not detectable across samples

idx.detectable <- detectable(tpm, 2.5);
sum(idx.detectable);

tpm.e <- tpm[idx.detectable, ];


# Assess 75% quantile of each sample

sample.q <- apply(tpm.e, 2, quantile, probs=0.75);

ggplot(melt(sample.q), aes(value)) + geom_histogram(bins=30);

low.quartile <- sample.q == 0;
sum(low.quartile);
low.quartile.samples <- colnames(tpm.e)[low.quartile];


# Assess proportion of zero measurements of each sample

sample.pzero <- prop_zero(tpm.e);
pheno$prop_zero <- sample.pzero;

ggplot(melt(sample.pzero), aes(value)) + geom_histogram(binwidth=0.05);


# Assess correlation between DV200 (proportion of RNA > 200 bp as reported by
# BioAnalyzer), and measures derived from RNA-seq libraries

qc.pdf.fname <- filename("qc", ext="pdf");

g <- ggplot(pheno, aes(x=dv200, y=prop_zero)) +
	geom_point();
qdraw(g, insert(qc.pdf.fname, tag=c("dv200-vs-propzero")));

g <- ggplot(pheno, aes(x=dv200, y=sum_rpkm)) +
	geom_point();
qdraw(g, insert(qc.pdf.fname, tag=c("dv200-vs-sumrpkm")));

g <- ggplot(pheno, aes(x=prop_zero, y=sum_rpkm)) +
	geom_point();
qdraw(g, insert(qc.pdf.fname, tag=c("propzero-vs-sumrpkm")));


# Perhaps more reads will help


# Filter samples

#high.pzero <- sample.pzero > 0.8;
#high.pzero <- sample.pzero > 0.5;
#high.pzero <- sample.pzero > 0.25;
high.pzero <- sample.pzero > 0.1;

sum(high.pzero);

high.pzero.samples <- colnames(tpm.e)[high.pzero];


# Assess density distributions of measurements

# In original scale

tpm.m <- melt(tpm.e);
colnames(tpm.m) <- c("ensembl", "sample", "tpm");
tpm.m$high_pzero <- FALSE;
tpm.m$high_pzero[tpm.m$sample %in% low.quartile.samples] <- TRUE;

g <- ggplot(tpm.m, aes(x=tpm, fill=sample)) + geom_density(alpha=0.2) + xlim(0, 1e2) +
	facet_grid(high_pzero ~ .);
g

# In log scale

tpm.l.m <- melt(log2(tpm.e+1));
colnames(tpm.l.m) <- tpm.m.colnames;
tpm.l.m$high_pzero <- FALSE;
tpm.l.m$high_pzero[tpm.l.m$sample %in% high.pzero.samples] <- TRUE;
tpm.l.m$fresh_frozen <- FALSE;
tpm.l.m$fresh_frozen[tpm.l.m$sample %in% fresh.frozen.samples] <- TRUE;

g <- ggplot(tpm.l.m, aes(x=log2_tpm, fill=sample)) + geom_density(alpha=0.2) +
	facet_grid(high_pzero ~ ., scales="free_y");
g
g <- ggplot(tpm.l.m, aes(x=log2_tpm, fill=sample)) + geom_density(alpha=0.2) +
	facet_grid(fresh_frozen ~ ., scales="free_y");
g


# Filter samples based on pzero

tpm.f <- tpm.e[, !high.pzero];
pheno.f <- pheno[!high.pzero, ];


# Assess box plot after sample filtering

tpm.f.l <- log2(tpm.f+1);

g <- ggplot(melt(tpm.f.l), aes(x=Var2, y=value)) + geom_boxplot() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1));
g

# Assess PCA plots to identify further outliers

tpm.pca <- pca(tpm.f.l);
g <- pca_plot(tpm.pca, pheno.f, aes(colour=met_type, shape=histology_group));
g

# Identify outliers from PCA

pca.outlier <- tpm.pca$Z[1, ] < -200 | tpm.pca$Z[2,] < -100;
colnames(tpm.pca$Z)[pca.outlier]

# Filter out outliers from PCA

tpm.f.pf <- tpm.f[, !pca.outlier];
pheno.f.pf <- pheno.f[!pca.outlier, ];


# Re-assess PCA plots after sample filtering

tpm.f.pf.l <- log2(tpm.f.pf+1);
tpm.pca <- pca(tpm.f.pf.l);

g <- pca_plot(tpm.pca, pheno.f.pf, aes(colour=met_type, shape=histology_group));
g
g <- pca_plot(tpm.pca, pheno.f.pf, aes(colour=met_type, shape=histology_group), dims=2:3);
g
g <- pca_plot(tpm.pca, pheno.f.pf, aes(colour=met_type, shape=histology_group), dims=3:4);
g


# Evaluate the patient-specific effect

summary(lmer(tpm.pca$Z[1,] ~ met_type + (1|patient_id), pheno.f.pf));


# Quantile normalization

tpm.f.qn <- normalize.quantiles(tpm.f.pf);
dimnames(tpm.f.qn) <- dimnames(tpm.f.pf);

# Assess box plot after quantile normalization

tpm.f.qn.l <- log2(tpm.f.qn+1);

tpm.f.qn.l.m <- melt(tpm.f.qn.l);
colnames(tpm.f.qn.l.m) <- tpm.m.colnames;

g <- ggplot(tpm.f.qn.l.m, aes(x=sample, y=log2_tpm)) + geom_boxplot() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1));
g

# Assess density plot after quantile normalization

g <- ggplot(tpm.f.qn.l.m, aes(x=log2_tpm, fill=sample)) + geom_density(alpha=0.2);
g


# Assess PCA plots

tpm.qn.pca <- pca(tpm.f.qn.l);
g <- pca_plot(tpm.qn.pca, pheno.f.pf, aes(colour=met_type, shape=histology_group));
g
g <- pca_plot(tpm.qn.pca, pheno.f.pf, aes(colour=met_type, shape=histology_group), dims=2:3);
g
g <- pca_plot(tpm.qn.pca, pheno.f.pf, aes(colour=met_type, shape=histology_group), dims=3:4);
g

g <- pca_plot(tpm.qn.pca, pheno.f.pf, aes(colour=patient_id, shape=met_type));
g
g <- pca_plot(tpm.qn.pca, pheno.f.pf, aes(colour=patient_id, shape=met_type), dims=2:3);
g
g <- pca_plot(tpm.qn.pca, pheno.f.pf, aes(colour=patient_id, shape=met_type), dims=3:4);
g

dv200_scaled <- pheno$dv200 / max(pheno$dv200, na.rm=TRUE);
g <- pca_plot(tpm.qn.pca, pheno.f.pf,
	aes(colour=patient_id, shape=met_type, alpha=dv200_scaled));
g

# Assess RNA quality with PCA component

g <- ggplot(pheno.f.pf, aes(x=dv200, y=tpm.qn.pca$Z[1, ])) + geom_point();
g
g <- ggplot(pheno.f.pf, aes(x=dv200, y=tpm.qn.pca$Z[2, ])) + geom_point();
g
g <- ggplot(pheno.f.pf, aes(x=dv200, y=tpm.qn.pca$Z[3, ])) + geom_point();
g

# Assess correlation of RNA with each PC

ks <- 1:nrow(tpm.qn.pca$Z);
dv200.cor <- unlist(lapply(ks, function(k) cor(tpm.qn.pca$Z[k, ], pheno.f.pf$dv200)));
dv200.cor.p <- unlist(lapply(ks, function(k) cor.test(tpm.qn.pca$Z[k, ], pheno.f.pf$dv200)$p.value));
g <- ggplot(, aes(x=ks, y=dv200.cor, fill=-log10(dv200.cor.p))) + geom_bar(stat="identity");
g


# Assess association of a PC with other covariates

k <- 1;
#k <- 2;

fit.null <- lm(tpm.qn.pca$Z[k,] ~ 1, pheno.f.pf);
fit.dv200 <- lm(tpm.qn.pca$Z[k,] ~ dv200, pheno.f.pf);
fit.pt <- lm(tpm.qn.pca$Z[k,] ~ patient_id, pheno.f.pf);
fit.met <- lm(tpm.qn.pca$Z[k,] ~ met_type, pheno.f.pf);
fit.met.pt <- lm(tpm.qn.pca$Z[k,] ~ met_type + patient_id + met_type:patient_id, pheno.f.pf);
fit.met.dv200 <- lm(tpm.qn.pca$Z[k,] ~ met_type + dv200, pheno.f.pf);

anova(fit.null, fit.dv200, test="Chisq");
anova(fit.null, fit.pt, test="Chisq");
anova(fit.null, fit.met, test="Chisq");
anova(fit.pt, fit.met.pt, test="Chisq");
anova(fit.dv200, fit.met.dv200, test="Chisq");

summary(lmer(tpm.qn.pca$Z[k,] ~ met_type + (1|patient_id), pheno.f.pf));
summary(lmer(tpm.qn.pca$Z[k,] ~ (1|patient_id) + (1|met_type), pheno.f.pf));

summary(lmer(tpm.qn.pca$Z[k,] ~ dv200 + (1|patient_id), pheno.f.pf));


# Hierarchical clustering

tpm.qn.hc <- hclust(as.dist(1 - cor(tpm.f.qn.l)), method="average");
plot(tpm.qn.hc);


# Output normalized and filtered expression data

#qwrite(tpm.f.pf.l, insert(out.fname, c("tpm", "f")));
#qwrite(tpm.f.qn.l, insert(out.fname, c("tpm", "f", "qn")));
#qwrite(rpkm, insert(out.fname, c("rpkm")));
qwrite(pheno, insert(out.fname, "qc"));

