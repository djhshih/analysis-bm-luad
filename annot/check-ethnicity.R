library(io);
library(dplyr);

s.pheno <- qread("sample-info_wes_stage2.tsv");
p.pheno <- qread("patient-info.tsv");

pops.df <- qread("../ancesinfer/summary/pkb-bm_exac-pop.tsv");

is_unique <- function(x) {
	length(unique(x)) == 1
}

s.pop <- group_by(pops.df, clinical_id) %>% select(clinical_id, exac_pop) %>%
	summarize(exac_pop=first(exac_pop), exac_pop_samples_agree=is_unique(exac_pop)) %>% ungroup();
summary(s.pop$exac_pop_samples_agree)

patients <- unique(s.pheno$clinical_id);

p.pheno.eth <- select(p.pheno, clinical_id, race)
with(p.pheno.eth, table(race), useNA="always")

sp.pheno <- left_join(s.pop, p.pheno.eth, by="clinical_id") %>%
	filter(clinical_id %in% patients)

with(sp.pheno, table(exac_pop, race, useNA="always"))
# there are no discrepancies

filter(sp.pheno, exac_pop == "SAS")
filter(p.pheno, clinical_id == "PB0135")

filter(sp.pheno, exac_pop == "OTH")
# ExAC population did not cluster with other major populations

with(sp.pheno, table(exac_pop))
with(sp.pheno, table(exac_pop) / sum(!is.na(exac_pop)))

