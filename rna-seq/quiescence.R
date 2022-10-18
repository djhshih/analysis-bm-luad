library(io);
library(RColorBrewer);
library(magrittr);
library(limma);
library(dplyr);
library(ggplot2);
library(kmplot);
library(survival);

source("../annot/df_helper.R")

x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");

pat.pheno <- qread("../annot/patient-info.tsv");
#pat.pheno$rnaseq <- pat.pheno$patient_id %in% x$pheno$patient_id;
#qwrite(pat.pheno, "../annot/Clinical_Data_Absolute_Outputs_03_01_2016.tsv");

censor_survival <- function(x, max.time) {
	idx <- x$os_year > max.time;
	x$os_years[idx] <- max.time;
	x$death[idx] <- 0;
	x
}

remove_rare_levels <- function(x, n) {
	d <- table(x);
	factor(x, levels=names(d)[d >= n])
}


pat.pheno %<>% mutate(
	gender = factor(gender, levels=c("M", "F")),
	histotype = normalize_empty(remove_rare_levels(histotype, 10)),
	smoking_hx_gt5 = factor(smoking_hx > 5)
);
pat.pheno$death[is.na(pat.pheno$os_years)] <- NA;

pat.pheno.c <- censor_survival(pat.pheno, 10);

# FIXME kmplot's at risk is shifted when genders contains a "" level

surv <- with(pat.pheno.c, Surv(os_years, death) ~ gender);
kmplot(surv);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ histotype);
kmplot(surv);

#surv <- with(pat.pheno.c, Surv(os_years, death) ~ histotype + gender);
#kmplot(surv);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ factor(smoker));
kmplot(surv);

with(pat.pheno, table(smoker, histotype));

#surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoker + histology);
#kmplot(surv);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoking_hx);
coxph(surv);

hist(pat.pheno.c$smoking_hx, breaks=12);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoking_hx_gt5);
survdiff(surv);
kmplot(surv);
plot(survfit(surv));


# Brain metastases

met.pheno <- x$pheno[x$pheno$met_type == "BM", ]
table(met.pheno$patient_id)

# ignore the first metastasis of patient PB-LS-020
met.pheno <- met.pheno[met.pheno$pair_id != "PB-LS-020-TM-NT-SM-7M18Z-SM-7M18X", ];

pheno.met.idx <- x$pheno$met_type == "BM" & x$pheno$pair_id != "PB-LS-020-TM-NT-SM-7M18Z-SM-7M18X";

met.idx <- match(met.pheno$patient_id, pat.pheno$patient_id);
met.pheno$treated <- pat.pheno$tx_prior_to_bm[met.idx];
met.pheno$cytotoxic_agents <- pat.pheno$number_of_cytotoxic_agents[met.idx];
met.pheno$t_stage <- pat.pheno$stage_of_dx_of_primary_t_stage[met.idx];
met.pheno$n_treatments <- pat.pheno$number_of_tx_regimens_prior_to_bm[met.idx];
met.pheno$os_years <- pat.pheno.c$os_years[met.idx];
met.pheno$death <- pat.pheno.c$death[met.idx];
met.pheno$bm_pfs_years <- pat.pheno.c$bm_pfs_years[met.idx];
met.pheno$bm <- 1;

met.pheno$quiescent <- factor(met.pheno$quiescence > 0);
met.pheno$genome_doubled <- factor(met.pheno$genome_doublings > 0);

ggplot(met.pheno, aes(x=treated, y=quiescence)) + geom_jitter(width=0.1);
ggplot(met.pheno, aes(x=cytotoxic_agents, y=quiescence)) + geom_point();
ggplot(met.pheno, aes(x=t_stage, y=quiescence)) + geom_point();
ggplot(met.pheno, aes(x=n_treatments, y=quiescence)) + geom_point();

surv <- with(met.pheno, Surv(os_years, death) ~ quiescence);
coxph(surv);

surv <- with(met.pheno, Surv(os_years, death) ~ quiescence);
coxph(surv);

surv <- with(met.pheno, Surv(bm_pfs_years, bm) ~ quiescence);
coxph(surv);

surv <- with(met.pheno, Surv(os_years, death) ~ quiescent);
survdiff(surv)
kmplot(surv);

surv <- with(met.pheno, Surv(bm_pfs_years, bm) ~ quiescent);
kmplot(surv);


hist(x$expr["PLOD2", ], breaks=12);
plod2 <- factor(x$expr["PLOD2", ] > median(x$expr["PLOD2", ]));
plod2.met <- plod2[pheno.met.idx];



expr.met <- x$expr[, pheno.met.idx];

surv <- with(met.pheno, Surv(os_years, death) ~ plod2.met);
kmplot(surv);
surv <- with(met.pheno, Surv(bm_pfs_years, bm) ~ plod2.met);
kmplot(surv);

with(met.pheno, table(plod2.met, bm_pfs_years < 0))
with(met.pheno, fisher.test(plod2.met, bm_pfs_years < 0))

stripchart(expr.met["PLOD2", ] ~ bm_pfs_years < 0, met.pheno)
wilcox.test(expr.met["PLOD2", ] ~ bm_pfs_years < 0, met.pheno)


surv <- with(met.pheno, Surv(os_years, death) ~ genome_doubled);
kmplot(surv);
surv <- with(met.pheno, Surv(bm_pfs_years, bm) ~ genome_doubled);
kmplot(surv);

with(met.pheno, table(genome_doubled, bm_pfs_years < 0))
with(met.pheno, fisher.test(genome_doubled, bm_pfs_years < 0))


# Primary tumours

pri.pheno <- x$pheno[x$pheno$met_type == "Primary", ]
table(met.pheno$patient_id)

pheno.pri.idx <- x$pheno$met_type == "Primary";

pri.idx <- match(pri.pheno$patient_id, pat.pheno$patient_id);
pri.pheno$t_stage <- pat.pheno[["Stage of diagnosis of primary tumor (T stage)"]][pri.idx];
pri.pheno$os_years <- pat.pheno.c$os_years[pri.idx];
pri.pheno$death <- pat.pheno.c$death[pri.idx];
pri.pheno$bm_pfs_years <- pat.pheno.c$bm_pfs_years[pri.idx];
pri.pheno$bm <- 1;

pri.pheno$genome_doubled <- factor(pri.pheno$genome_doublings > 0);
pri.pheno$quiescent <- factor(pri.pheno$quiescence > 0);

surv <- with(pri.pheno, Surv(os_years, death) ~ quiescent);
survdiff(surv);
#kmplot(surv);

surv <- with(pri.pheno, Surv(bm_pfs_years, bm) ~ quiescent);
survdiff(surv);
kmplot(surv);
with(pri.pheno, table(quiescent, bm_pfs_years < 0))
with(pri.pheno, fisher.test(quiescent, bm_pfs_years < 0))

surv <- with(pri.pheno, Surv(os_years, death) ~ genome_doubled);
survdiff(surv);
#kmplot(surv);

surv <- with(pri.pheno, Surv(bm_pfs_years, bm) ~ genome_doubled);
survdiff(surv);
kmplot(surv);

with(pri.pheno, table(genome_doubled, bm_pfs_years < 0))
with(pri.pheno, fisher.test(genome_doubled, bm_pfs_years < 0))

plod2.pri <- plod2[pheno.pri.idx];
surv <- with(pri.pheno, Surv(os_years, death) ~ plod2.pri);
kmplot(surv);
surv <- with(pri.pheno, Surv(bm_pfs_years, bm) ~ plod2.pri);
kmplot(surv);

###!
with(pri.pheno, table(pri_plod_high=plod2.pri, bm_at_diag=bm_pfs_years < 0))
with(pri.pheno, fisher.test(plod2.pri, bm_pfs_years < 0))

stripchart(expr.pri["PLOD2", ] ~ bm_pfs_years < 0, pri.pheno)
plot(expr.pri["PLOD2", ] ~ as.factor(bm_pfs_years < 0), pri.pheno)
wilcox.test(expr.pri["PLOD2", ] ~ as.factor(bm_pfs_years < 0), pri.pheno)
t.test(expr.pri["PLOD2", ] ~ as.factor(bm_pfs_years < 0), pri.pheno)


expr.pri <- x$expr[, pheno.pri.idx];

plot(expr.pri["MYC", ] ~ as.factor(bm_pfs_years < 0), pri.pheno)
wilcox.test(expr.pri["MYC", ] ~ as.factor(bm_pfs_years < 0), pri.pheno)


hist(x$expr["MYC", ])
myc <- factor(x$expr["MYC", ] > median(x$expr["MYC", ]));
myc.pri <- myc[pheno.pri.idx];

with(pri.pheno, table(myc.pri, bm_pfs_years < 0))
with(pri.pheno, fisher.test(myc.pri, bm_pfs_years < 0))


is.met <- x$pheno$met_type == "BM";
plot(x$pheno$quiescence, x$expr["GFAP", ]);
plot(x$pheno$quiescence[is.met], x$expr["GFAP", is.met]);


#pheno <- x$pheno %>%
#	select(patient_id, tumor_id, met_type) %>%
#	mutate(has_maf = tumor_id %in% maf$Tumor_Sample_Barcode);

#qwrite(as.character(pheno$tumor_id[!pheno$has_maf]), file="missing_maf.vtr")

mafs <- lapply(
	x$pheno$pair_id,
	function(pair_id) {
		qread(sprintf("../exome-seq/maf/%s.union_patient_calls.annotated.maf", pair_id))
	}
);
names(mafs) <- x$pheno$tumor_id;


cannonical_driver_events <- list(
	KRAS = c("p.G12A", "p.G12C", "p.G12D", "p.G12F", "p.G12R", "p.G12S", "p.G12V", "p.G12Y", "p.Q61L"),
	BRAF = c("p.G466A", "p.G469L", "p.G466V", "p.G469V", "p.N581S", "p.D594H", "p.D594N",  "p.V600E"),
	EGFR = c("p.ELREA701del", "p.I706T", "p.G719A", "p.L858R", "p.728_729insH" ),
	ERBB2 = c("p.S310F", "p.774_775insAYVM", "p.776_776G>VC", "p.V777L", "H.amp"),
	MAP2K1 = c("p.K57M", "p.K57N", "p.C121S"),
	MET = c("p.D963_splice", "p.Y1003*", "p.D1010_splice", "H.amp"),
	HRAS = c("p.Q61L"),
	NRAS = c("p.Q61L")
);

get_samples_with_event_common_maf <- function(gene, protein_change, maf) {
	as.character((
		maf %>%
			filter(Hugo_Symbol == gene) %>%
			filter(Protein_Change %in% protein_change) %>%
			select(Tumor_Sample_Barcode)
	)[, 1])
}

make_mutations_common_maf <- function(samples, genes, protein_changes, maf) {
	names(genes) <- genes;
	d <- as.matrix(as.data.frame(mapply(
		function(gene, protein_change) {
			samples %in% get_samples_with_event_common_maf(gene, protein_change, maf)
		},
		genes,
		protein_changes
	)));
	
	# determine covered samples
	has.maf = samples %in% maf$Tumor_Sample_Barcode;

	# blank out non-covered samples
	d[!has.maf, ] <- NA;

	rownames(d) <- samples;

	d
}

get_samples_with_event_common_maf <- function(gene, protein_change, maf) {
	as.character((
		maf %>%
			filter(Hugo_Symbol == gene) %>%
			filter(Protein_Change %in% protein_change) %>%
			select(Tumor_Sample_Barcode)
	)[, 1])
}

make_mutations <- function(samples, genes, protein_changes, mafs) {
	names(genes) <- genes;
	d <- as.matrix(as.data.frame(mapply(
		function(gene, protein_change) {
			samples_mutated <- unlist(lapply(mafs,
				function(maf) {
					get_samples_with_event_common_maf(gene, protein_change, maf)
				}
			));
			samples %in% samples_mutated
		},
		genes,
		protein_changes
	)));

	rownames(d) <- samples;

	d
}

muts <- make_mutations(x$pheno$tumor_id, names(cannonical_driver_events), cannonical_driver_events, mafs);

apply(muts, 2, sum);

plot(quiescence ~ KRAS, cbind(x$pheno, muts));
summary(lm(quiescence ~ KRAS, cbind(x$pheno, muts)));
t.test(quiescence ~ KRAS, cbind(x$pheno, muts));

plot(quiescence ~ EGFR, cbind(x$pheno, muts));
summary(lm(quiescence ~ EGFR, cbind(x$pheno, muts)));
t.test(quiescence ~ EGFR, cbind(x$pheno, muts));

muts.f <- apply(muts, 2, function(x) factor(as.numeric(x), levels=c(0, 1), labels=c("non-mutated", "mutated")));
muts.gene <- unlist(apply(muts, 1,
	function(x) {
		idx <- which(x);
		if (length(idx) > 0) {
			colnames(muts)[which(x)]
		} else {
			"none"
		}
	}
));

rtk.path.mut <- factor(apply(muts, 1, any), labels=c("non-mutated", "mutated"));
plot(quiescence ~ rtk.path.mut, x$pheno);
d <- data.frame(
	quiescence = x$pheno$quiescence,
	rtk_pathway = rtk.path.mut,
	driver = muts.gene
);
g <- ggplot(d, aes(x=rtk_pathway, y=quiescence, colour=driver)) + geom_jitter(width=0.1);
qdraw(g, file="rtk-path-mut-vs-quiescence.pdf", width=4);


summary(lm(quiescence ~ rtk.path.mut, x$pheno));
t.test(quiescence ~ rtk.path.mut, x$pheno);

###!

with(pri.pheno, table(rtk.path.mut[pheno.pri.idx], quiescent))
with(pri.pheno, fisher.test(rtk.path.mut[pheno.pri.idx], quiescent))

with(met.pheno, table(rtk.path.mut[pheno.met.idx], quiescent))
with(met.pheno, fisher.test(rtk.path.mut[pheno.met.idx], quiescent))

with(pri.pheno, table(rtk.path.mut[pheno.pri.idx], bm_pfs_years < 0))
with(pri.pheno, fisher.test(rtk.path.mut[pheno.pri.idx], bm_pfs_years < 0))

with(met.pheno, table(rtk.path.mut[pheno.met.idx], bm_pfs_years < 0))
with(met.pheno, fisher.test(rtk.path.mut[pheno.met.idx], bm_pfs_years < 0))

with(pri.pheno, table(muts[pheno.pri.idx, "KRAS"], bm_pfs_years < 0))
with(pri.pheno, fisher.test(muts[pheno.pri.idx, "KRAS"], bm_pfs_years < 0))

with(met.pheno, table(muts[pheno.met.idx, "KRAS"], bm_pfs_years < 0))
