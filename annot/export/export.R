library(io);
library(dplyr);
library(lubridate);

pheno <- qread("../sample-info_wes_stage3_pass_luad.tsv");
clin <- qread("../patient-info_stage2.tsv");

pheno.out.fname <- filename("sample-info", tag=c("pass", "luad", "export"), ext="tsv");
clin.out.fname <- filename("patient-info", tag=c("pass", "luad", "export"), ext="tsv");

pheno$pick[pheno$sample_id %in% c("LS-016-M1", "LS-016-P-2")] <- TRUE;

pheno.e <- select(pheno,
	sample_id, wes_platform, center, clinical_id,
	primary_cancer_type, primary_histotype,
	sample_type, specimen_material, sex,
	purity, ploidy, genome_doublings,
	coverage_for_80pct_power, ffpe_q, picard_oxoq,
	somatic_mutation_covered_bases_capture,
	contamination_percentage_consensus_capture,
	pass_qc,
	mutation_burden, mutation_burden_indel, pick
) %>% filter(sample_type %in% c("Brain metastasis", "Primary", "Normal"));
#) %>% filter(pick, sample_type %in% c("Brain metastasis", "Primary", "Normal"));

print(nrow(clin.e))
print(nrow(filter(pheno.e, sample_type == "Brain metastasis")))
print(nrow(filter(pheno.e, sample_type == "Primary")))

clin.e <- transmute(clin,
	clinical_id,
	age_at_primary_dx,
	gender, sex,
	race, exac_pop,
	smoker, smoking_signature_high,
	primary_cancer_type, primary_histotype,
	stage, tumor_stage, node_status, met_status,
	year_of_primary_dx = year(date_of_primary_dx),
	year_of_bm_dx = year(date_of_bm_dx),
	os_years, bm_pfs_years, 
	synchronous_bm, 
) %>% filter(clinical_id %in% pheno.e$clinical_id);

# look for discrepancies

with(clin.e, table(gender, sex))
filter(clin.e, gender == "Female", sex == "XY")

with(clin.e, table(race, exac_pop))

with(clin.e, cbind(year_of_bm_dx - year_of_primary_dx, bm_pfs_years))
with(clin.e, max(abs((year_of_bm_dx - year_of_primary_dx) - bm_pfs_years)))

# remove gender for now
clin.e <- select(clin.e, -gender);

qwrite(pheno.e, pheno.out.fname);
qwrite(clin.e, clin.out.fname);

