library(gpldiff)
library(io);

load_all("~/projects/r/gpldiff");

#d <- qread("gistic-chr8q-amp_bm-luad_tcga-luad.rds");
#d <- qread("gistic-chr11q-amp_bm-luad_tcga-luad.rds");
#d <- qread("gistic-chr4q-amp_bm-luad_tcga-luad.rds");
#d <- qread("gistic-chr14q-amp_bm-luad_tcga-luad.rds");
#d <- qread("gistic-chr5p-amp_bm-luad_tcga-luad.rds");
d <- qread("gistic-chr19q-amp_bm-luad_tcga-luad.rds");
# g has to be in {-0.5, 0.5}

hparams <- list(
	nu2 = 10,
	lambda2 = 1^2,
	alpha = 10,
	beta = 1,
	tau2 = 0.1^2
);

#params <- list(
#	mu = 0,
#	sigma2 = 0.1,
#	f = rep(0, d$J)
#);
params <- NULL;

fit <- gpldiff(d, NULL, hparams, adapt=FALSE);

str(fit);

qdraw(
	{
		plot(fit, d);
	},
	height = 6,
	width = 10,
	#file = "gpldiff_chr8q-amp_bm-luad_tcga-luad.pdf"	
	#file = "gpldiff_chr11q-amp_bm-luad_tcga-luad.pdf"	
	#file = "gpldiff_chr4q-amp_bm-luad_tcga-luad.pdf"	
	#file = "gpldiff_chr14q-amp_bm-luad_tcga-luad.pdf"	
	#file = "gpldiff_chr5p-amp_bm-luad_tcga-luad.pdf"	
	file = "gpldiff_chr19q-amp_bm-luad_tcga-luad.pdf"	
);


qdraw(
	{
		hist(fit$params$f, breaks=50);
	},
	#file = "gpldiff_f-hist_chr8q-amp_bm-luad_tcga-luad.pdf"
	#file = "gpldiff_f-hist_chr11q-amp_bm-luad_tcga-luad.pdf"
	#file = "gpldiff_f-hist_chr4q-amp_bm-luad_tcga-luad.pdf"
	#file = "gpldiff_f-hist_chr14q-amp_bm-luad_tcga-luad.pdf"
	#file = "gpldiff_f-hist_chr5p-amp_bm-luad_tcga-luad.pdf"
	file = "gpldiff_f-hist_chr19q-amp_bm-luad_tcga-luad.pdf"
);

