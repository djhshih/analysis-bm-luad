library(io);
library(dplyr);

clin <- qread("patient-info.rds");
ances <- qread("../ancesinfer/summary/pkb-bm_exac-pop.tsv");
emu <- qread("../emu/pkb-luad_smoking.tsv");
sex <- qread("../sex/brain-mets_normals_sex.tsv");
pheno <- qread("sample-info_wes.tsv");

is_unique <- function(x) {
	length(unique(x)) == 1
}

ances.g <- group_by(ances, clinical_id) %>% select(-sample_id) %>% 
	summarize(exac_pop=first(exac_pop), agree=is_unique(exac_pop)) %>% ungroup();

filter(ances.g, !agree)
stopifnot(ances.g$agree)

sex <- mutate(sex,
	clinical_id = pheno$clinical_id[match(sample_id, pheno$sample_id)]
) %>% select(-sample_id);

clin.m <- left_join(clin, select(ances.g, -agree), by="clinical_id") %>%
	left_join(select(emu, -smoker), by="clinical_id") %>%
	left_join(sex, by="clinical_id");

out.fname <- filename("patient-info", tag="stage2");

qwrite(clin.m, insert(out.fname, ext="tsv"), quote=TRUE);
qwrite(clin.m, insert(out.fname, ext="rds"));

