library(io);
library(dplyr);

# Analyze matched BM and primary with candidate CNA drivers
# using the mid p binomial test (e.g. a fixed-effect model)

clin <- qread("~/exigens/brain-mets/annot/patient-info_stage2.rds");
pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv");
patients <- as.character(unique(pheno$clinical_id));
clin <- clin[match(patients, clin$clinical_id), ];

cn <- qread("../comut/brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);

mut.prims <- qread("../comut/pkb-luad_prims_black-f.pset.mmaf", type="tsv", quote="\"");
mut.mets <- qread("../comut/pkb-luad_brain-mets_black-f.pset.mmaf", type="tsv", quote="\"");

out.fname <- filename("pkb-luad", tag="matched-cmp");

source("mid-p-binom.R");

mid_p_binomial_test(4, 4)
mid_p_binomial_test(4, 4, mid.p=FALSE)
binom.test(4, 4)

mid_p_binomial_test(4, 4, alternative = "greater")
mid_p_binomial_test(4, 4, mid.p=FALSE, alternative = "greater")
binom.test(4, 4,alternative = "greater")



amp.cut <- 8;

which(cn["MMP13", ] > amp.cut)
cn["MMP13", grep("PB0067", colnames(cn))]
cn["MMP13", grep("PB0086", colnames(cn))]
mid_p_binomial_test(4, 5, alternative = "greater")

which(cn["YAP1", ] > amp.cut)
cn["YAP1", grep("PB0067", colnames(cn))]
cn["YAP1", grep("PB0086", colnames(cn))]
mid_p_binomial_test(2, 2, alternative = "greater")

which(cn["MYC", ] > amp.cut)
cn["MYC", grep("LS-004", colnames(cn))]
cn["MYC", grep("LS-016", colnames(cn))]
cn["MYC", grep("PB0031", colnames(cn))]
cn["MYC", grep("PB0034", colnames(cn))]
cn["MYC", grep("PB0040", colnames(cn))]
cn["MYC", grep("PB0046", colnames(cn))]
cn["MYC", grep("PB0098", colnames(cn))]
cn["MYC", grep("PB0160", colnames(cn))]
cn["MYC", grep("PB0251", colnames(cn))]
cn["MYC", grep("PB0376", colnames(cn))]

mid_p_binomial_test(4, 4, alternative = "greater")


del.cut <- 0.5;
which(cn["CDKN2A", ] < del.cut)

samples <- names(which(cn["CDKN2B", ] < del.cut));
samples.e <- as.character(filter(pheno, clinical_id %in% filter(pheno, sample_id %in% samples)$clinical_id, sample_type %in% c("Brain metastasis", "Primary"))$sample_id);
data.frame(sample_id=samples.e, cdkn2b=cn["CDKN2B", samples.e])

mid_p_binomial_test(3, 3 + 1, alternative = "greater")

mid_p_binomial_test(3, 3 + 13, alternative = "two.sided")


# combined analysis

# BM-specific vs. Prim-specific
print(mid_p_binomial_test(13, 13 + 4));
# aberrations are more likely to be BM-specific rather than Prim-specific

# BM-late vs. BM early
print(mid_p_binomial_test(13, 13 + 24));
# aberrations trend toward begin early events
# (often shared betwen primary and metastasis)

