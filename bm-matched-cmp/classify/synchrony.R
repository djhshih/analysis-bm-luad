library(io)
library(dplyr)


calls <- qread("bm-comparts_calls.tsv");
pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
clin <- qread("~/exigens/brain-mets/annot/patient-info_stage2.rds");

patients <- as.character(unique(filter(pheno, pick)$clinical_id));
clin.f <- filter(clin, clinical_id %in% patients);
clin.called <- left_join(clin.f, calls);

with(clin.called, table(synchronous_bm, bm_cdkn2a_del))
fisher.test(with(clin.called, table(synchronous_bm, bm_cdkn2a_del)))

with(clin.called, table(synchronous_bm, bm_cdkn2b_del))
fisher.test(with(clin.called, table(synchronous_bm, bm_cdkn2b_del)))

with(clin.called, table(synchronous_bm, pr_cdkn2a_del))
fisher.test(with(clin.called, table(synchronous_bm, pr_cdkn2a_del)))

with(clin.called, table(synchronous_bm, pr_cdkn2b_del))
fisher.test(with(clin.called, table(synchronous_bm, pr_cdkn2b_del)))

