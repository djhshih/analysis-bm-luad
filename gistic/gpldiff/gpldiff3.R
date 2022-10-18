library(gpldiff)
library(io);

options(mc.cores=8);

load_all("~/projects/r/gpldiff");

# Version used in initial BM-LUAD manuscript of 2018-06

control.fname <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_absolute_cnvf_wt/scores.gistic";
case.fname <- "../luad_absolute_bmet-only_cnvf/scores.gistic";
out.fname <- filename("pkb-tcga-luad", tag="gpldiff");

hparams <- list(
	nu2 = 5^2,
	lambda2 = 1^2,
	alpha = 10,
	beta = 1,
	tau2 = 1^2
);

params <- NULL;

case <- read_gistic(case.fname);
control <- read_gistic(control.fname);

fits <- compare_gistics(case.fname, control.fname, params, hparams);
qwrite(fits, insert(out.fname, ext="rds"));
#fits <- qread("pkb-tcga-luad_gpldiff.rds");

results <- summary.gistic_gpldiffs(fits);

results.f <- results[results$diff > 0.6 & results$fdr < 0.05, ];

qwrite(results, insert(out.fname, ext="tsv"));

