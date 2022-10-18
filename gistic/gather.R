library(io);
library(dplyr);
library(reshape2);

absolute.dir <- "~/share/exigens/brain-mets/absolute/ABSOLUTE_results/brain-mets_pass_luad_absolute-1-4/reviewed";

segs <- qread(file.path(absolute.dir, "SEG_MAF"), pattern="seg.tsv");

# make seg data.frame from ABSOLUTE seg output
make_seg <- function(x) {
	clonal <- !is.na(x$subclonal.a1) & x$subclonal.a1 == 0 &
		!is.na(x$subclonal.a2) & x$subclonal.a2 == 0;
	tcn <- x$rescaled_total_cn;
	# round clonal estimates to the nearest integer, except for amplified regions
	# it is difficult to estimate clonality of highly amplified events due to
	# ambiguouity between two models: many cells harbour modest level
	# amplication and few cells harbour high level amplification
	idx <- tcn < 8 & clonal;
	tcn[idx] <- round(tcn[idx]);

	with(x,
		data.frame(
			sample = sample,
			chromosome = Chromosome,
			start = Start.bp,
			end = End.bp,
			n_probes = n_probes,
			rescaled_total_cn = tcn
		)
	)
}

seg <- do.call(rbind, lapply(segs, make_seg));
rownames(seg) <- NULL;

# rescaled_total_cn = copy ratio / delta - b
# where basebline copy ratio b and step size delta are calculated based on
# inferred tau (tumour ploidy) and alpha (purity)
# it is an estimate of the absolute copy-number with measurement noise and
# possibility of subclonal copy-number events

#seg <- qread(file.path(absolute.dir, "brain-mets_pass_luad_absolute-1-4_rescaled_total_cn.IGV.seg.txt"), type="tsv");

absolute <- qread(file.path(absolute.dir, "absolute_brain-mets_pass_luad.tsv"));

samples.failed <- as.character(absolute$sample[
	absolute[["call status"]] == "low purity" | absolute[["call status"]] == "FAILED"
]);

seg.f <- seg[! seg$sample %in% samples.failed, ];

# divide the absolute copy-number by the ploidy
ploidy <- absolute$ploidy[match(seg.f$sample, absolute$sample)];
# lower bound by 0.01 to avoid NaN after log2 transformation
seg.f$rescaled_total_cn <- log2(pmax(seg.f$rescaled_total_cn / ploidy, 0.01));

qwrite(seg.f, "brain-mets_pass_luad_absolute.lgr.seg");


pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv");

# argmax_x y
# @param idx  subset index
argmax <- function(x, y, idx=NULL) {
	if (is.null(idx)) {
		idx <- 1:length(x);
	}
	i <- which.max(y[idx]);
	if (length(i) == 0) {
		NA
	} else {
		x[idx][i]
	}
}

# select the sample with the highest purity for each sample
patients <- lapply(
	split(pheno, pheno$clinical_id),
	function(x) {
		list(
			normal = as.character(argmax(x$sample_id, x$contamination_percentage_consensus_capture, x$sample_type == "Normal")),
			primary = as.character(argmax(x$sample_id, x$purity, x$sample_type == "Primary")),
			brain_metastasis = as.character(argmax(x$sample_id, x$purity, x$sample_type == "Brain metastasis"))
		)
	}
);


# the brain metastasis of these two patients failed ABSOLUTE...
# should they be excluded from the cohort?
filter(pheno, clinical_id == "PB0164")
# PB0164-M failed because of problem with AllelicCapseg hypersegmentation
# it does have clear detectable CNAs

samples <- melt(patients);
colnames(samples) <- c("sample_id", "sample_type", "clinical_id");

omit_na <- function(x) {
	x[!is.na(x)]
}

samples.bmet <- omit_na(as.character(filter(samples, sample_type == "brain_metastasis")$sample_id));
samples.prim <- omit_na(as.character(filter(samples, sample_type == "primary")$sample_id));
#samples.normal <- omit_na(as.character(filter(samples, sample_type == "normal")$sample_id));


seg.f.bmet <- filter(seg.f, sample %in% samples.bmet);
seg.f.prim <- filter(seg.f, sample %in% samples.prim);


qwrite(seg.f.bmet, "brain-mets_pass_luad_absolute_bmet-only.lgr.seg");
qwrite(seg.f.prim, "brain-mets_pass_luad_absolute_prim-only.lgr.seg");

