library(io);
library(dplyr);

in.fname <- as.filename("../brain-mets_pass_luad_absolute_bmet-only.lgr.cntf.cenf.cnvf.seg");
out.fname <- set_fpath(insert(in.fname, ext="wt"), ".");


weights.df <- qread("~/exigens/brain-mets/matchit/wstage/pkb-luad_weights.tsv");
pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes.tsv");
weights.df <- left_join(weights.df, select(pheno, clinical_id, sample_id));

seg <- qread(in.fname);

idx <- !is.na(weights.df$weight) & weights.df$weight > 0;
samples <- as.character(weights.df$sample_id[idx]);

seg <- seg[seg$sample %in% samples, ];

weights <- weights.df$weight[match(seg$sample, weights.df$sample_id)];
#weights[is.na(weights)] <- 0;

seg$state <- seg$state * weights;

print(summary(seg$state))

qwrite(seg, out.fname);
