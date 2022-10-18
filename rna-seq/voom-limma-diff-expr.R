library(limma);
library(io);
library(edgeR);

x <- qread("2016-01-21/luad-bm_tpm_f_20160121T172905.rds");
pheno.all <- qread("../annot/sample-info.tsv");
pheno <- pheno.all[match(colnames(x), pheno.all$rnaseq_id), ];

met.type <- factor(pheno$met_type, levels=c("Primary", "BM"));
design <- model.matrix(~ met.type);

# Estimate voom parameters and using it during linear fitting for identifying
# differentially expressed genes

# It is not ideal that the input data is TPM instead of counts...

dge <- DGEList(counts=2^x-1);
dge <- calcNormFactors(dge);
v <- voom(dge, design, plot=TRUE, normalize="quantile");

fit <- lmFit(v, design);
eb <- eBayes(fit);
sig.df <- topTable(eb, number=Inf);
sig.df <- cbind(gencode=rownames(sig.df), sig.df);
qwrite(sig.df, filename("luad-bm", tag=c("voom-qn", "limma", "met-vs-prim"), ext="tsv"));

