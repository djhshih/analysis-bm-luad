library(io)

load_all("~/projects/r/gpldiff");

options(mc.cores=8);

control.path <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/wpurity/tcga-luad_absolute_cnvf_wt";
control.fname <- file.path(control.path, "scores.gistic");
control.peaks.fname <- file.path(control.path, "all_lesions.conf_95.txt");
case.fname <- "../luad_absolute_bmet-only_cnvf/scores.gistic";
case.peaks.fname <- "../luad_absolute_bmet-only_cnvf/all_lesions.conf_95.txt";
out.fname <- filename("pkb-tcga-luad", tag="gpldiff");

hparams <- list(
	nu2 = 5^2,
	lambda2 = 1^2,
	alpha = 0.1,
	beta = 0.1,
	tau2 = 1^2
);

params <- NULL;

case <- read_gistic(case.fname);
control <- read_gistic(control.fname);

case.peaks <- read_gistic_peaks(case.peaks.fname);
control.peaks <- read_gistic_peaks(control.peaks.fname);

fits <- compare_gistics(case.fname, control.fname, params, hparams);
qwrite(fits, insert(out.fname, ext="rds"));
#fits <- qread("pkb-tcga-luad_gpldiff.rds");

results <- summary.gistic_gpldiffs(fits);
results <- summary_append_gistic_peaks(results, case.peaks)
results.f <- results[which(results$diff > 0.5 & results$fdr < 0.05 & results$gistic_q < 0.05), ];
print(results.f)

results.down <- summary.gistic_gpldiffs(fits, direction=-1);
results.down <- summary_append_gistic_peaks(results.down, control.peaks)
results.down.f <- results.down[which(results.down$diff < -0.5 & results.down$fdr < 0.05 & results.down$gistic_q < 0.05), ];
print(results.down.f)

qwrite(results, insert(out.fname, ext="tsv"));
qwrite(results.f, insert(out.fname, tag="filt", ext="tsv"));
qwrite(results.down, insert(out.fname, tag="down", ext="tsv"));
qwrite(results.down.f, insert(out.fname, tag=c("down", "filt"), ext="tsv"));

