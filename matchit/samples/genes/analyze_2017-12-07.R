library(io);
library(dplyr);

options(stringsAsFactors=FALSE);

cohorts <- c("TCGA", "Present");

patients <- as.character(unique(qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv")$clinical_id));
patients.pkb <- patients;

pheno.pkb.all <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3.tsv");
pheno.pkb <- filter(pheno.pkb.all, clinical_id %in% patients);

pheno.tcga <- qread("~/exigens/tcga/tcga-luad/annot/sample-info_uniq_pass_tcga-luad_stage3.tsv");

mut.genes <- as.character(qread("~/exigens/brain-mets/comut/mutsigcv_sig-genes.vtr"));
cn.genes <- as.character(qread("~/exigens/brain-mets/comut/gistic_sig-genes_manual.vtr"));

out.fname <- filename("pkb-luad", date=NA);

cn <- qread("~/exigens/brain-mets/comut/brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);


mut.prims <- qread("~/exigens/brain-mets/comut/pkb-luad_prims_black-f_20171031T113017.pset.mmaf", type="tsv", quote="\"");
mut.bmets <- qread("~/exigens/brain-mets/comut/pkb-luad_brain-mets_black-f_20171031T112959.pset.mmaf", type="tsv", quote="\"");

mut.tcga <- qread("~/exigens/brain-mets/comut/tcga-luad/tcga-luad_pass_black-f_20171113T105434.pset.mmaf", type="tsv", quote="\"");

weights.tcga <- qread("~/exigens/brain-mets/matchit/tcga-luad_weights.tsv");
weights.pkb <- qread("~/exigens/brain-mets/matchit/pkb-luad_weights.tsv");

pheno.tcga <- left_join(pheno.tcga, weights.tcga, by="clinical_id");
pheno.pkb <- left_join(pheno.pkb, weights.pkb, by="clinical_id");

patients.tcga <- as.character(unique(pheno.tcga$clinical_id));

weights.tcga <- weights.tcga[match(patients.tcga, weights.tcga$clinical_id), ];
weights.pkb <- weights.pkb[match(patients.pkb, weights.pkb$clinical_id), ];

G <- length(genes);
N <- length(patients);
N.tcga <- length(patients.tcga);

idx.tcga <- match(mut.tcga$Tumor_Sample_Barcode, pheno.tcga$fh_sample_id);
mut.tcga$sample_id <- pheno.tcga$sample_id[idx.tcga];
mut.tcga$clinical_id <- pheno.tcga$clinical_id[idx.tcga];
mut.tcga$weight <- pheno.tcga$weight[idx.tcga];

decompose_samples <- function(x) {
	unlist(lapply(strsplit(x, ":"), function(z) z[1]));
}

idx.bmets <- match(decompose_samples(mut.bmets$Tumor_Sample_Barcode), pheno$fh_sample_id);
mut.bmets$sample_id <- pheno$sample_id[idx.bmets];
mut.bmets$clinical_id <- pheno$clinical_id[idx.bmets];
mut.bmets$weight <- pheno$weight[idx.bmets];


lof.classes <- c("Missense_Mutation", "Nonsense_mutation", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation", "Start_Codon_Del", "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins");

hi.classes <- lof.classes;

gof.classes <- c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins");


####

# samples from the same patient have already been summarized together
# however, a patient may have two mutations in the same gene
# therefore, the below frequency estimates are biased

sum(filter(mut.tcga, Hugo_Symbol == "KRAS", Variant_Classification %in% gof.classes)$weight)
sum(filter(mut.bmets, Hugo_Symbol == "KRAS", Variant_Classification %in% gof.classes)$weight)

sum(filter(mut.tcga, Hugo_Symbol == "EGFR", Variant_Classification %in% gof.classes)$weight)
sum(filter(mut.bmets, Hugo_Symbol == "EGFR", Variant_Classification %in% gof.classes)$weight)


filter(mut.bmets, Hugo_Symbol %in% gene)

#gene <- "EGFR";
#gene <- "KRAS";
#gene <- "SPTAN1";
#gene <- "NSD1";
gene <- "RPTN";



d.tcga <- filter(mut.tcga, Hugo_Symbol %in% gene, Variant_Classification %in% hi.classes);
d.bmets <- filter(mut.bmets, Hugo_Symbol %in% gene, Variant_Classification %in% hi.classes);

s <- rbind(
	data.frame(
		cohort = "TCGA",
		mutation = patients.tcga %in% d.tcga$clinical_id,
		weight_raw = weights.tcga$weight,
		weight = weights.tcga$weight_norm
	),
	data.frame(
		cohort = "Present",
		mutation = patients %in% d.bmets$clinical_id,
		weight_raw = weights.pkb$weight,
		weight = weights.pkb$weight_norm
	)
);
s$cohort <- factor(s$cohort, levels=cohorts);

with(s, table(mutation, cohort))
with(s, fisher.test(table(mutation, cohort)))

fit <- glm(mutation ~ cohort, data = s, family=binomial("logit"));
summary(fit);

fit <- glm(mutation ~ cohort, data = s, family=quasibinomial);
summary(fit);

# calculate significance
fit <- glm(mutation ~ cohort, weights = weight_raw, data = s, family=quasibinomial);
summary(fit);

fit <- glm(mutation ~ cohort, weights = weight, data = s, family=quasibinomial);
summary(fit);

fit <- glm(mutation ~ cohort, weights = weight, data = s, family=binomial("logit"));
summary(fit);

# get t confidence intervals
#fit <- glm(mutation ~ cohort - 1, weights = weight, data = s, family=binomial("logit"));
fit <- glm(mutation ~ cohort - 1, weights = weight_raw, data = s, family=quasibinomial);
lest <- summary(fit)$coefficients[, 1];
lse <- summary(fit)$coefficients[, 2];
exp(lest - qnorm(1 - 0.20/2) * lse)
exp(lest)
exp(lest + qnorm(1 - 0.20/2) * lse)

sum(filter(s, cohort == "TCGA", mutation)$weight)
sum(filter(s, cohort == "Present", mutation)$weight)

####

