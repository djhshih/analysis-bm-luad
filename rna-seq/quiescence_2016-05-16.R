library(io);
library(RColorBrewer);
library(magrittr);
library(limma);
library(dplyr);
library(ggplot2);
library(kmplot);
library(survival);

x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");

pat.pheno <- qread("../annot/Clinical_Data_Absolute_Outputs_01_20_2016.txt", type="tsv");
pat.pheno$rnaseq <- pat.pheno$patient_id %in% x$pheno$patient_id;
qwrite(pat.pheno, "../annot/Clinical_Data_Absolute_Outputs_03_01_2016.tsv");




date.src.format <- "%m/%d/%y";

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


pat.pheno <- mutate(
	pat.pheno,
	death = as.numeric(!grepl("alive|visit", tolower(pat.pheno[["Date of Death"]]))),
	date_of_death = as.Date(
		gsub("[A-Za-z ]", "", as.character(pat.pheno[["Date of Death"]])),
		date.src.format
	),
	date_of_diagnosis = as.Date(
		as.character(pat.pheno[["Date of Diagnosis of Primary Tumor"]]),
		date.src.format
	),
	os = date_of_death - date_of_diagnosis,
	os_years = as.numeric(os) / 356,
	histology = remove_rare_levels(pat.pheno[["Histology group"]], 10),
	gender = factor(Gender, levels=c("M", "F")),
	smoker = factor(pat.pheno[["Smoking history"]]),
	smoking_pack_year = pat.pheno[["Smoking Pack-Year History"]],
	smoking_pack_year_gt40 = factor(smoking_pack_year > 40)
);
pat.pheno$death[is.na(pat.pheno$os_years)] <- NA;

pat.pheno.c <- censor_survival(pat.pheno, 10);

# FIXME kmplot's at risk is shifted when genders contains a "" level

surv <- with(pat.pheno.c, Surv(os_years, death) ~ gender);
kmplot(surv);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ histology);
kmplot(surv);

#surv <- with(pat.pheno.c, Surv(os_years, death) ~ histology + gender);
#kmplot(surv);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoker);
kmplot(surv);

with(pat.pheno, table(smoker, histology));

#surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoker + histology);
#kmplot(surv);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoking_pack_year);
coxph(surv);

hist(pat.pheno.c$smoking_pack_year, breaks=12);

surv <- with(pat.pheno.c, Surv(os_years, death) ~ smoking_pack_year_gt40);
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
met.pheno$treated <- pat.pheno[["treatment prior to BM"]][met.idx];
met.pheno$cytotoxic_agents <- pat.pheno[["Number of cytotoxic agents"]][met.idx];
met.pheno$t_stage <- pat.pheno[["Stage of diagnosis of primary tumor (T stage)"]][met.idx];
met.pheno$n_treatments <- pat.pheno[["Number of treatment regimens prior to brain metastasis"]][met.idx];
met.pheno$os_years <- pat.pheno.c$os_years[met.idx];
met.pheno$death <- pat.pheno.c$death[met.idx];

ggplot(met.pheno, aes(x=treated, y=quiescence)) + geom_jitter(width=0.1);
ggplot(met.pheno, aes(x=cytotoxic_agents, y=quiescence)) + geom_point();
ggplot(met.pheno, aes(x=t_stage, y=quiescence)) + geom_point();
ggplot(met.pheno, aes(x=n_treatments, y=quiescence)) + geom_point();

surv <- with(met.pheno, Surv(os_years, death) ~ quiescence);
coxph(surv);

met.pheno$quiescent <- factor(met.pheno$quiescence > 0);
surv <- with(met.pheno, Surv(os_years, death) ~ quiescent);
survdiff(surv)
kmplot(surv);


hist(x$expr["PLOD2", ], breaks=12);
plod2 <- factor(x$expr["PLOD2", ] > median(x$expr["PLOD2", ]));
plod2.met <- plod2[pheno.met.idx];
surv <- with(met.pheno, Surv(os_years, death) ~ plod2.met);
kmplot(surv);

met.pheno$genome_doubled <- factor(met.pheno$genome_doublings > 0);
surv <- with(met.pheno, Surv(os_years, death) ~ genome_doubled);
kmplot(surv);


# Primary tumours

pri.pheno <- x$pheno[x$pheno$met_type == "Primary", ]
table(met.pheno$patient_id)

pri.idx <- match(pri.pheno$patient_id, pat.pheno$patient_id);
pri.pheno$t_stage <- pat.pheno[["Stage of diagnosis of primary tumor (T stage)"]][pri.idx];
pri.pheno$os_years <- pat.pheno.c$os_years[pri.idx];
pri.pheno$death <- pat.pheno.c$death[pri.idx];

pri.pheno$genome_doubled <- factor(pri.pheno$genome_doublings > 0);
pri.pheno$quiescent <- factor(pri.pheno$quiescence > 0);

surv <- with(pri.pheno, Surv(os_years, death) ~ quiescent);
survdiff(surv);
#kmplot(surv);

surv <- with(pri.pheno, Surv(os_years, death) ~ genome_doubled);
survdiff(surv);
#kmplot(surv);


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
			samples %in% get_samples_with_event(gene, protein_change, maf)
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
					get_samples_with_event(gene, protein_change, maf)
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

