library(io);
library(RColorBrewer);
library(magrittr);
library(dplyr);
library(ggplot2);
library(kmplot);
library(survival);


x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");
maf <- qread("../exome-seq/unpaired_Bmet_LUAD_and_carcinoma.maf", type="tsv", quote="");


d <- table(maf$pair_id);
d <- data.frame(pair_id=names(d), nmuts=as.numeric(d));
median(d$nmuts);

g <- ggplot(d, aes(x=pair_id, y=nmuts)) + geom_bar(stat="identity") + theme_bw();
g

