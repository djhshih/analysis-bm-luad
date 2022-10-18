library(io);
library(data.table);
library(mmalign);

x <- qread("2016-01-21/luad-bm_tpm_f_qn_20160121T172905.rds");
annot <- qread("~/data/biomart/release-83/ensembl-id_hgnc-symbol_complete.tsv");

out.fname <- filename("luad-bm", ext="rds");

pheno.all <- qread("../annot/sample-info.tsv");
pheno <- pheno.all[match(colnames(x), pheno.all$rnaseq_id), ];

ensembl <- gsub("\\.\\d+_\\d+$", "", rownames(x));
symbol <- annot$hgnc_symbol[match(ensembl, annot$ensembl_gene_id)];

idx.mappable <- !is.na(symbol);

sum(!idx.mappable);

x.sub <- x[idx.mappable, ];

x.dt <- data.table(x);
x.dt.agg <- x.dt[, lapply(.SD, max, na.rm=TRUE), by=list(gene=symbol)];

x.mat <- as.matrix(x.dt.agg[, 2:ncol(x.dt.agg), with=FALSE]);
rownames(x.mat) <- x.dt.agg$gene;

x.pca <- pca(x.mat);
pca_plot(x.pca, pheno, aes(colour=patient_id, shape=met_type));
pca_plot(x.pca, pheno, aes(colour=patient_id, shape=met_type), dims=2:3);
pca_plot(x.pca, pheno, aes(colour=patient_id, shape=met_type), dims=3:4);

qwrite(x.mat, insert(out.fname, c("tpm", "f", "qn", "gene")));

