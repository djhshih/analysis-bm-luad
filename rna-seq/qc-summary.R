library(io);

expr.norm <- qread("luad-bm_tpm_f_qn_gene.rds");
expr.raw <- qread("luad-bm_rpkm_gene.rds");

pheno <- qread("../annot/sample-info.tsv");

# standardize all sample names

map_sample_ids <- function(x, pheno, from, to) {
	idx <- match(x, pheno[, from]);
	stopifnot(all(!is.na(idx)));
	as.character(pheno[idx, to])
}

samples.raw <- map_sample_ids(colnames(expr.raw), pheno, "rnaseq_id", "pair_id");
samples.norm <- map_sample_ids(colnames(expr.norm), pheno, "rnaseq_id", "pair_id");

s <- data.frame(
	pair_id = samples.raw,
	pzero_pass = samples.raw %in% samples.norm
);

qwrite(s, "luad-bm_qc-summary.tsv");
