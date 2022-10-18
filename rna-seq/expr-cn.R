# Assess corrleation between copy-number and expression of specific genes

library(io);
library(ggplot2);
library(reshape2);

expr <- qread("luad-bm_tpm_f_qn_gene.rds");
cn <- qread("../cn/combined_gene_corrected_CN.mtx");
mut <- qread("../exome-seq/luad-bm_sig-snvs_exome-seq_mutsigcv.tsv", quote="");
pheno <- qread("../annot/sample-info.tsv");

# standardize all sample names

map_sample_ids <- function(x, pheno, from, to) {
	idx <- match(x, pheno[, from]);
	stopifnot(all(!is.na(idx)));
	pheno[idx, to]
}

colnames(expr) <- map_sample_ids(colnames(expr), pheno, "rnaseq_id", "pair_id");

# cn has some unannoated samples, but its IDs are already pair_id
#colnames(cn) <- map_sample_ids(colnames(cn), pheno, "pair_id", "pair_id");

mut$pair_id <- map_sample_ids(mut$Tumor_Sample_Barcode, pheno, "tumor_id", "pair_id");

samples <- intersect(colnames(expr), colnames(cn));

# samples with expr data but no cn data
setdiff(colnames(expr), colnames(cn));

expr.n <- expr[, match(samples, colnames(expr))];
cn.n <- cn[, match(samples, colnames(cn))];
pheno.n <- pheno[match(samples, pheno$pair_id), ];

# annotate pheno.n with mutation information
# assume wildtype if no mutation in reported in gene!
make_mut <- function() {
	mut.n <- data.frame(pair_id=pheno.n$pair_id);
	idx <- match(samples, mut$pair_id)	
	idx.valid <- !is.na(idx);
	pair.ids <- mut$pair_id[idx[idx.valid]]	
	genes <- as.character(mut$Hugo_Symbol[idx[idx.valid]]);
	mut.types <- as.character(mut$Variant_Classification[idx[idx.valid]]);
	for (i in 1:length(genes)) {
		mut.n[, genes[i]] <- "Wildtype";
		mut.n[match(pair.ids[i], mut.n$pair_id), genes[i]] <- mut.types[i];
	}

	#x <- t(mut.n[, 2:ncol(mut.n)]);
	#colnames(x) <- mut.n$pair_id;
	mut.n
}

mut.n <- make_mut();

make_gene <- function(gene) {
	data.frame(
		snv=mut.n[, gene],
		expr=expr.n[gene, ],
		cn=cn.n[gene, ],
		met_type=pheno.n$met_type
	)
}

genes <- colnames(mut.n)[-1];
for (gene in genes) {
	print(gene);
	try({
		g <- ggplot(make_gene(gene), aes(x=log2(cn), y=expr, shape=snv, colour=met_type)) + 
			geom_point() + ggtitle(gene);
		qdraw(g, filename("cn-expr-snv", tag=gene, ext="pdf"));
	});
}

make_gene_expr_cn <- function(gene) {
	data.frame(
		expr=expr.n[gene, ],
		cn=cn.n[gene, ],
		met_type=pheno.n$met_type
	)
}

other.genes <- c(setdiff(levels(mut$Hugo_Symbol), genes), "YAP1", "MYC");
yet.other.genes <- c("RUNX1", "CRYAB", "EDNRA");
my.other.genes <- c("SALL1", "FKBP5", "GFAP", "PMP2", "SLC1A3", "GAP43");
for (gene in c(other.genes, yet.other.genes, my.other.genes)) {
	print(gene);
	try({
		g <- ggplot(make_gene_expr_cn(gene), aes(x=log2(cn), y=expr, colour=met_type)) + 
			geom_point() + ggtitle(gene);
		qdraw(g, filename("cn-expr", tag=gene, ext="pdf"));
	});
}

cryab <- make_gene_expr_cn("CRYAB");
t.test(expr ~ met_type, cryab);

pik <- c(grep("PIK", rownames(expr), value=TRUE), "PTEN");

make_genes_expr_cn <- function(genes) {
	data.frame(
		expr=c(expr.n[genes, ]),
		cn=c(cn.n[genes, ]),
		gene=rep(genes, times=ncol(expr.n)),
		met_type=rep(pheno.n$met_type, each=length(genes))
	)
}

g <- ggplot(make_genes_expr_cn(pik), aes(x=log2(cn), y=expr, colour=met_type)) + 
	geom_point() + facet_wrap(~gene);
qdraw(g, filename("cn-expr", tag="pik", ext="pdf"), width=8, height=6);

mmp <- c(intersect(
	grep("^MMP", rownames(expr), value=TRUE),
	grep("^MMP", rownames(cn), value=TRUE)
));

g <- ggplot(make_genes_expr_cn(c(mmp, "YAP1")), aes(x=log2(cn), y=expr, colour=met_type)) + 
	geom_point() + facet_wrap(~gene);
qdraw(g, filename("cn-expr", tag="mmp", ext="pdf"), width=8, height=6);

g <- ggplot(make_genes_expr_cn(mmp), aes(x=met_type, y=expr, colour=gene)) + geom_jitter();
qdraw(g, filename("cn-expr", tag=c("mmp", "jitter"), ext="pdf"));

output <- list(cn=cn.n, expr=expr.n, pheno=pheno.n, mut=mut.n);
qwrite(output, filename("luad-bm", tag=c("cn-expr-snv-pheno"), ext="rds"));

