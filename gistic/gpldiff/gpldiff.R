library(gpldiff)
library(io);

load_all("~/projects/r/gpldiff");

#d <- qread("gistic-chr8q-amp_bm-luad_tcga-luad.rds");
#d <- qread("gistic-chr8q-amp_bm-luad_tcga-luad_ends.rds");
d <- qread("gistic-chr8-amp_bm-luad_tcga-luad.rds");

#idx <- d$x > 50;
#d$J <- sum(idx);
#d$x <- d$x[idx];
#d$g <- d$g[idx];
#d$y <- d$y[idx];

hparams <- list(
	nu2 = 5^2,
	lambda2 = 1^2,
	alpha = 10,
	beta = 1,
	tau2 = 1^2
);

fit <- gpldiff(d, NULL, hparams, adapt=FALSE);

str(fit)

prob <- summary(fit);
lodds <- log10(prob) - log10(1 - prob);
idx <- which(lodds > 2)
d$x[idx]

qdraw(
	{
		plot(fit, d);
	},
	height= 6,
	width = 10,
	#file = "gpldiff_chr8q-amp_bm-luad_tcga-luad_ends.pdf"	
	#file = "gpldiff_chr8q-amp_bm-luad_tcga-luad.pdf"	
	file = "gpldiff_chr8-amp_bm-luad_tcga-luad.pdf"	
);

idx <- which(d$x > 12.8 & d$x < 12.9) + -5:5;

data.frame(
	d$x[idx],
	d$y[idx],
	d$g[idx],
	fit$params$f[idx, ]
)

