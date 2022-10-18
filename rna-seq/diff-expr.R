# Identify differentially expressed genes between primary and met

library(io);
library(limma);
library(lme4);

expr <- qread("luad-bm_tpm_f_qn_gene.rds");
pheno.all <- qread("../annot/sample-info.tsv");
pheno <- pheno.all[match(colnames(expr), pheno.all$rnaseq_id), ];

met.type <- factor(pheno$met_type, levels=c("Primary", "BM"));
design <- model.matrix(~ met.type);

# Remove that pesky row with an invalid gene name
# (Where did it come from?)
expr <- expr[!is.na(rownames(expr)), ];

fit <- lmFit(expr, design);
eb <- eBayes(fit);
sig.df <- topTable(eb, number=Inf);
sig.df <- cbind(symbol=rownames(sig.df), sig.df);

out.fname <- filename("luad-bm", tag=c("limma", "met-vs-prim"));

qwrite(sig.df, insert(out.fname, ext="tsv"));
qwrite(sig.df, insert(out.fname, ext="rds"));


rnk <- data.frame(
	gene = rownames(sig.df),
	t = sig.df$t
);
write.table(rnk, file=tag(out.fname, ext="rnk"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");

# Promising hits: SALL1 (EMT), FKBP5, CRYAB

# Notes:
# The apparent upregulated expression of GFAP, PMP2, SLC1A3, GAP43 in the metastasis
# suggests significant stroma contamination
# It may also suggest tumour-cell mediated remodeling of stroma

# Downregulation of FCMR may mean that the metastasis is not infiltrated by
# immune cells

genes <- rownames(sig.df)[sig.df$adj.P.Val < 0.25];

# The association of metastasis of the top genes cannot be explained by 
# the patient covariate: their upregulation is generalizable across patients.

for (gene in genes) {
	print(gene);
	fit <- lmer(expr[gene, ] ~ met.type + (1|pheno$patient_id));
	print(summary(fit));
	cat('\n');
}

