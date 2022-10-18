library(io);
library(dplyr);
library(ggplot2);

options(stringsAsFactors=FALSE);

cohorts <- c("TCGA", "Present");

pdf.fname <- filename("pkb-tcga-luad-compare", ext="pdf");

patients.pkb <- as.character(unique(qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv")$clinical_id));

pheno.pkb.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3.tsv");
pheno.pkb <- filter(pheno.pkb.all, clinical_id %in% patients);

pheno.tcga <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");

#mut.genes <- as.character(qread("~/exigens/brain-mets/comut/mutsigcv_sig-genes.vtr"));
mut.genes <- c("RPTN", "PABPC1", "FAHD2B", "SIRPA", "PKP1", "TPRA1", "CYBB", "NSD1", "SPTAN1", "MMP13", "SERPINE1");
names(mut.genes) <- mut.genes;

cn.genes <- as.character(qread("~/exigens/brain-mets/comut/gistic_sig-genes_manual.vtr"));
names(cn.genes) <- cn.genes;

del.genes <- c("CDKN2A", "CDKN2B");
names(del.genes) <- del.genes;

amp.genes <- setdiff(cn.genes, del.genes);
names(amp.genes) <- amp.genes;

out.fname <- filename("pkb-luad", date=NA);

samples.prims <- pheno.pkb$sample_id[pheno.pkb$sample_type == "Primary"];
samples.bmets <- pheno.pkb$sample_id[pheno.pkb$sample_type == "Brain metastasis"];

cn <- qread("~/exigens/brain-mets/comut/brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);
cn.prims <- cn[, colnames(cn) %in% samples.prims];
cn.bmets <- cn[, colnames(cn) %in% samples.bmets];

mut.prims <- qread("~/exigens/brain-mets/comut/pkb-luad_prims_black-f_20171031T113017.pset.mmaf", type="tsv", quote="\"");
mut.bmets <- qread("~/exigens/brain-mets/comut/pkb-luad_brain-mets_black-f_20171031T112959.pset.mmaf", type="tsv", quote="\"");

mut.tcga <- qread("~/exigens/brain-mets/comut/tcga-luad/tcga-luad_pass_black-f_20171113T105434.pset.mmaf", type="tsv", quote="\"");
cn.tcga <- qread("~/exigens/brain-mets/comut/tcga-luad/tcga-luad_pass_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn.tcga <- as.matrix(cn.tcga);

weights.tcga <- qread("~/exigens/brain-mets/matchit/tcga-luad_weights.tsv");
weights.pkb <- qread("~/exigens/brain-mets/matchit/pkb-luad_weights.tsv");

pheno.tcga <- left_join(pheno.tcga, weights.tcga, by="clinical_id");
pheno.pkb <- left_join(pheno.pkb, weights.pkb, by="clinical_id");

patients.tcga <- as.character(unique(pheno.tcga$clinical_id));

weights.tcga <- weights.tcga[match(patients.tcga, weights.tcga$clinical_id), ];
weights.pkb <- weights.pkb[match(patients.pkb, weights.pkb$clinical_id), ];

N <- length(patients.pkb);
N.tcga <- length(patients.tcga);

idx.tcga <- match(mut.tcga$Tumor_Sample_Barcode, pheno.tcga$fh_sample_id);
mut.tcga$sample_id <- pheno.tcga$sample_id[idx.tcga];
mut.tcga$clinical_id <- pheno.tcga$clinical_id[idx.tcga];
mut.tcga$weight <- pheno.tcga$weight[idx.tcga];

decompose_samples <- function(x) {
	unlist(lapply(strsplit(x, ":"), function(z) z[1]));
}

idx.bmets <- match(decompose_samples(mut.bmets$Tumor_Sample_Barcode), pheno.pkb$fh_sample_id);
mut.bmets$sample_id <- pheno.pkb$sample_id[idx.bmets];
mut.bmets$clinical_id <- pheno.pkb$clinical_id[idx.bmets];
mut.bmets$weight <- pheno.pkb$weight[idx.bmets];


lof.classes <- c("Missense_Mutation", "Nonsense_mutation", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation", "Start_Codon_Del", "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins");

hi.classes <- lof.classes;

gof.classes <- c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins");

cols.ch <- c(Control="royalblue3", Case="darkorange");
cohorts <- c("TCGA", "Present");
names(cols.ch) <- cohorts; 

####

# samples from the same patient have already been summarized together
# however, a patient may have two mutations in the same gene
# therefore, the below frequency estimates are biased

sum(filter(mut.tcga, Hugo_Symbol == "KRAS", Variant_Classification %in% gof.classes)$weight)
sum(filter(mut.bmets, Hugo_Symbol == "KRAS", Variant_Classification %in% gof.classes)$weight)

sum(filter(mut.tcga, Hugo_Symbol == "EGFR", Variant_Classification %in% gof.classes)$weight)
sum(filter(mut.bmets, Hugo_Symbol == "EGFR", Variant_Classification %in% gof.classes)$weight)


filter(mut.bmets, Hugo_Symbol %in% gene)



compare_weighted_mut_freq <- function(mut.case, mut.control, weights.case, weights.control, gene, variant.classes, test=c("glm", "fisher"), alpha.level=0.20) {

	test <- match.arg(test);

	d.case <- filter(mut.case, Hugo_Symbol %in% gene, Variant_Classification %in% variant.classes);
	d.control <- filter(mut.control, Hugo_Symbol %in% gene, Variant_Classification %in% variant.classes);

	# use the weights directly output by cem: do not normalize
	s <- rbind(
		data.frame(
			cohort = "Control",
			mutation = weights.control$clinical_id %in% d.control$clinical_id,
			weight = weights.control$weight
		),
		data.frame(
			cohort = "Case",
			mutation = weights.case$clinical_id %in% d.case$clinical_id,
			weight = weights.case$weight
		)
	);
	s$cohort <- factor(s$cohort, levels=c("Control", "Case"));

	if (test == "fisher") {
		ct <- with(s, table(mutation, cohort));
		fit <- fisher.test(ct);

		list(
			model = fit,
			p.value = fit$p.value,
			estimate = fit$estimate
		)
	} else {
		# calculate significance
		# use quasibinomial to allow observations to be real numbers
		fit <- glm(mutation ~ cohort, weights = weight, data = s, family=quasibinomial);
		coefs <- summary(fit)$coefficients;

		# get t confidence intervals
		fit0 <- glm(mutation ~ cohort - 1, weights = weight, data = s, family=quasibinomial);
		coefs0 <- summary(fit0)$coefficients;
		lest <- coefs0[, 1];
		lse <- coefs0[, 2];

		coef.df <- data.frame(
			mean = exp(lest),
			lower = exp(lest - qnorm(1 - alpha.level/2) * lse),
			upper = exp(lest + qnorm(1 - alpha.level/2) * lse)
		)

		list(
			model = fit,
			p.value = coefs[2,4],
			estimate = coefs[2,1:2],
			coefficients = coef.df
		)
	}
}

mut.genes.res <- lapply(
	mut.genes,
	function(gene) {
		compare_weighted_mut_freq(mut.bmets, mut.tcga, weights.pkb, weights.tcga, gene, hi.classes)
	}
);

mut.genes.ps <- unlist(lapply(mut.genes.res, function(x) x$p.value));
mut.genes.qs <- p.adjust(mut.genes.ps, "BH");
which(mut.genes.qs < 0.25)

mut.genes.or <- unlist(lapply(mut.genes.res, function(x) exp(x$estimate[1])));

mut.genes.ps <- unlist(lapply(mut.genes.res, function(x) x$p.value));

mut.genes.coefs <- do.call(rbind,
	mapply(function(x, gene) {
		data.frame(
			gene = gene,
			cohort = gsub("cohort", "", rownames(x$coefficients)),
			x$coefficients
		)
	}, mut.genes.res, mut.genes, SIMPLIFY=FALSE));
rownames(mut.genes.coefs) <- NULL;
mut.genes.coefs$cohort <- factor(mut.genes.coefs$cohort, levels=c("Control", "Case"), labels=c("TCGA", "Present"));

alphav <- 0.8;

qdraw(
	ggplot(mut.genes.coefs, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) + 
		geom_col(alpha=alphav) +
		geom_errorbar(width=0.2) + theme_bw() + 
		#facet_wrap(~ gene, scales="free_y") +
		facet_wrap(~ gene) + ylim(0, 0.2) +
		scale_fill_manual(values=cols.ch) + 
		theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
		ylab("Somatic SNV frequency") + xlab("")
	,
	width = 6, height = 5,
	file = insert(pdf.fname, c("snv-freq", "bar"))
);

#filter(mut.bmets, Hugo_Symbol == "OR8H2")

####


# TODO support multiple genes
compare_weighted_cna_freq <- function(cn.case, cn.control, pheno.case, pheno.control, weights.case, weights.control, gene, cna.type=c("amp", "del", "either"), test=c("glm", "fisher"), alpha.level=0.20) { 
	test <- match.arg(test);
	cna.type <- match.arg(cna.type);

	amp.cut <- 8;
	del.cut <- 0.5;

	if (cna.type == "amp") {
		b.case <- cn.case[gene, ] > amp.cut;
		b.control <- cn.control[gene, ] > amp.cut;
	} else if (cna.type == "del") {
		b.case <- cn.case[gene, ] < del.cut;
		b.control <- cn.control[gene, ] < del.cut;
	} else {
		b.case <- cn.case[gene, ] > amp.cut | cn.case[gene, ] < del.cut;
		b.control <- cn.control[gene, ] > amp.cut | cn.control[gene, ] < del.cut;
	}

	patients.case <- pheno.case$clinical_id[match(names(b.case)[b.case], pheno.case$sample_id)];
	patients.control <- pheno.control$clinical_id[match(names(b.control)[b.control], pheno.control$sample_id)];

	# use the weights directly output by cem: do not normalize
	s <- rbind(
		data.frame(
			cohort = "Control",
			mutation = weights.control$clinical_id %in% patients.control,
			weight = weights.control$weight
		),
		data.frame(
			cohort = "Case",
			mutation = weights.case$clinical_id %in% patients.case,
			weight = weights.case$weight
		)
	);
	s$cohort <- factor(s$cohort, levels=c("Control", "Case"));

	if (test == "fisher") {
		ct <- with(s, table(mutation, cohort));
		fit <- fisher.test(ct);

		list(
			model = fit,
			p.value = fit$p.value,
			estimate = fit$estimate
		)
	} else {
		# calculate significance
		# use quasibinomial to allow observations to be real numbers
		fit <- glm(mutation ~ cohort, weights = weight, data = s, family=quasibinomial);
		coefs <- summary(fit)$coefficients;

		# get t confidence intervals
		fit0 <- glm(mutation ~ cohort - 1, weights = weight, data = s, family=quasibinomial);
		coefs0 <- summary(fit0)$coefficients;
		lest <- coefs0[, 1];
		lse <- coefs0[, 2];

		coef.df <- data.frame(
			mean = exp(lest),
			lower = exp(lest - qnorm(1 - alpha.level/2) * lse),
			upper = exp(lest + qnorm(1 - alpha.level/2) * lse)
		)

		list(
			model = fit,
			p.value = coefs[2,4],
			estimate = coefs[2,1:2],
			coefficients = coef.df
		)
	}
}


amp.genes.res <- lapply(
	amp.genes,
	function(gene) {
		compare_weighted_cna_freq(cn.bmets, cn.tcga, pheno.pkb, pheno.tcga, weights.pkb, weights.tcga, gene, "amp", alpha.level=0.20);
	}
);

amp.genes.ps <- unlist(lapply(amp.genes.res, function(x) x$p.value));
amp.genes.qs <- p.adjust(amp.genes.ps, "BH");
which(amp.genes.qs < 0.25)

amp.genes.or <- unlist(lapply(amp.genes.res, function(x) exp(x$estimate[1])));

amp.genes.coefs <- do.call(rbind,
	mapply(function(x, gene) {
		data.frame(
			gene = gene,
			cohort = gsub("cohort", "", rownames(x$coefficients)),
			x$coefficients
		)
	}, amp.genes.res, amp.genes, SIMPLIFY=FALSE));
rownames(amp.genes.coefs) <- NULL;
amp.genes.coefs$cohort <- factor(amp.genes.coefs$cohort, levels=c("Control", "Case"), labels=c("TCGA", "Present"));


qdraw(
	ggplot(amp.genes.coefs, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) + 
		geom_col(alpha=alphav) +
		geom_errorbar(width=0.2) + theme_bw() + 
		#facet_wrap(~ gene, scales="free_y") +
		facet_wrap(~ gene) + ylim(0, 0.25) +
		scale_fill_manual(values=cols.ch) + 
		theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
		ylab("Somatic CNA frequency") + xlab("")
	,
	width = 6, height = 5,
	file = insert(pdf.fname, c("amp-freq", "bar"))
);


del.genes.res <- lapply(
	del.genes,
	function(gene) {
		compare_weighted_cna_freq(cn.bmets, cn.tcga, pheno.pkb, pheno.tcga, weights.pkb, weights.tcga, gene, "del", alpha.level=0.20);
	}
);

del.genes.ps <- unlist(lapply(del.genes.res, function(x) x$p.value));
del.genes.qs <- p.adjust(del.genes.ps, "BH");
which(del.genes.qs < 0.25)

del.genes.or <- unlist(lapply(del.genes.res, function(x) exp(x$estimate[1])));

del.genes.coefs <- do.call(rbind,
	mapply(function(x, gene) {
		data.frame(
			gene = gene,
			cohort = gsub("cohort", "", rownames(x$coefficients)),
			x$coefficients
		)
	}, del.genes.res, del.genes, SIMPLIFY=FALSE));
rownames(del.genes.coefs) <- NULL;
del.genes.coefs$cohort <- factor(del.genes.coefs$cohort, levels=c("Control", "Case"), labels=c("TCGA", "Present"));


qdraw(
	ggplot(del.genes.coefs, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) + 
		geom_col(alpha=alphav) +
		geom_errorbar(width=0.2) + theme_bw() + 
		#facet_wrap(~ gene, scales="free_y") +
		facet_wrap(~ gene) + ylim(0, 0.55) +
		scale_fill_manual(values=cols.ch) + 
		theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
		ylab("Somatic CNA frequency") + xlab("")
	,
	width = 6, height = 5,
	file = insert(pdf.fname, c("del-freq", "bar"))
);

#hist(cn.bmets["CDKN2A",], breaks=100)
#which(cn.bmets["CDKN2A", ] > 0.1 & cn.bmets["CDKN2A", ] < 0.5)
#cn.bmets["CDKN2A", "PB0046-M"]


cna.genes <- c("MYC", "YAP1", "MMP13", "CDKN2A", "CDKN2B");

cna.genes.coefs <- rbind(amp.genes.coefs, del.genes.coefs) %>% filter(gene %in% cna.genes);
cna.genes.coefs$gene  <- factor(cna.genes.coefs$gene, levels=cna.genes);

qdraw(
	ggplot(cna.genes.coefs, aes(x=cohort, y=mean, ymin=lower, ymax=upper, fill=cohort)) + 
		geom_col(alpha=alphav) +
		geom_errorbar(width=0.2) + theme_bw() + 
		facet_grid(. ~ gene) + 
		scale_fill_manual(values=cols.ch) + 
		theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
		ylab("Somatic CNA frequency") + xlab("")
	,
	width = 6, height = 3,
	file = insert(pdf.fname, c("cna-freq", "bar"))
);

