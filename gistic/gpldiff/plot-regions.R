library(gpldiff)
library(io);

options(mc.cores=8);

load_all("~/projects/r/gpldiff");

control.fname <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_absolute_cnvf_wt/scores.gistic";
case.fname <- "../luad_absolute_bmet-only_cnvf/scores.gistic";
fits.fname <- "pkb-tcga-luad_gpldiff.rds";
out.fname <- filename("pkb-tcga-luad", tag="gpldiff");

fits <- qread(fits.fname);

####

# chromosome plots

ns <- strsplit(names(fits), ".", fixed=TRUE);

mapply(
	function(fit, type, chrom) {
		qdraw(
			{
				plot(fit$model, fit$data);
			},
			height = 6, width = 10,
			file = insert(out.fname, tag=c(sprintf("chr%s", chrom), tolower(type)), ext="pdf")
		)
	},
	fits,
	unlist(lapply(ns, function(x) x[1])),
	unlist(lapply(ns, function(x) x[2]))
);

####

# chr8:126445706-130762261
fit.myc <- subset_gpldiff(fits[["Amp.8"]], start = 126445706 / 1e6 - 1, end = 130762261 / 1e6 + 1);

qdraw(
	{
		plot(fit.myc$model, fit.myc$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("myc", "amp"), ext="pdf")
);

# chr11:101863542-102585362
fit.yap1 <- subset_gpldiff(fits[["Amp.11"]], start = 101863542 / 1e6 - 3, end = 102585362 / 1e6 + 3);

qdraw(
	{
		plot(fit.yap1$model, fit.yap1$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("yap1", "amp"), ext="pdf")
);

# chr14:36334948-39510012
fit.nkx21 <- subset_gpldiff(fits[["Amp.14"]], start = 36334948 / 1e6 - 3, end = 39510012 / 1e6 + 3);

qdraw(
	{
		plot(fit.nkx21$model, fit.nkx21$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("nkx2-1", "amp"), ext="pdf")
);

# chr19:38189869-40883989
fit <- subset_gpldiff(fits[["Amp.19"]], start = 38189869 / 1e6 - 3, end = 40883989 / 1e6 + 3);

qdraw(
	{
		plot(fit$model, fit$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("ccne1", "amp"), ext="pdf")
);

# chr1:149897699-153391702 (1q21)
fit <- subset_gpldiff(fits[["Amp.1"]], start = 149897699 / 1e6 - 3, end = 153391702 / 1e6 + 3);

qdraw(
	{
		plot(fit$model, fit$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("1q21", "amp"), ext="pdf")
);

# chr12:66935663-69645890
fit <- subset_gpldiff(fits[["Amp.12"]], start = 66935663 / 1e6 - 3, end = 69645890 / 1e6 + 3);

qdraw(
	{
		plot(fit$model, fit$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("mdm2", "amp"), ext="pdf")
);

# chr9:21227888-23693464
fit <- subset_gpldiff(fits[["Del.9"]], start = 21227888 / 1e6 - 5, end = 23693464 / 1e6 + 5);

qdraw(
	{
		plot(fit$model, fit$data)
	},
	width = 2,
	file = insert(out.fname, tag=c("cdkn2a-cdkn2b", "del"), ext="pdf")
);

