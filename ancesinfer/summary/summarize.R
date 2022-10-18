library(io);

pheno <- qread("~/share/exigens/brain-mets/annot/sample-info_wes_stage2.tsv");

ances.infer <- qread("~/share/exigens/brain-mets/ancesinfer", type="directory", pattern="_ancestry.tsv");
names(ances.infer) <- sub("_ancestry.tsv", "", names(ances.infer));

out.fname <- filename("pkb-bm");

pops <- unlist(lapply(ances.infer,
	function(x) {
		as.character(x$population[which.max(x$log_posterior)])
	}
));

pops.df <- data.frame(
	sample_id = names(pops),
	exac_pop = pops
);
rownames(pops.df) <- NULL;
pops.df$clinical_id <- pheno$clinical_id[match(pops.df$sample_id, pheno$sample_id)];

qwrite(ances.infer, insert(out.fname, "ancesinfer", ext="rds"));
qwrite(pops.df, insert(out.fname, "exac-pop", ext="tsv"));

