library(io);
library(dplyr);

source("df_helper.R");

options(stringsAsFactors=FALSE);

pheno.fname <- as.filename("sample-info_wes_stage2.tsv");
out.fname <- filename(pheno.fname$fstem, tag=pheno.fname$tag, ext=pheno.fname$ext);
out.fname$tag[length(out.fname$tag)] <- "stage3";

pheno <- qread(pheno.fname);
burden <- qread("../summary/burden.tsv");

out <- left_join(pheno, burden, by="sample_id");

qwrite(out, out.fname);

