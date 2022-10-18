# Compare sample level statistics between TCGA and present cohort

source("presembl.R")

####

variables <- c(
	ffpe_q = "FFPE quality score",
	picard_oxoq = "OxoG quality score",
	contamination_percentage_consensus_capture = "Cross-contamination estimate (%)",
	covered_mb = "Covered bases (Mbp)",
	mutation_burden = "SNVs / Mbp",
	mutation_burden_exome = "Exome SNVs / Mbp",
	mutation_burden_indel = "Indels / Mbp",
	purity = "Tumor purity",
	ploidy = "Tumor ploidy",
	subclonal_genome_fraction = "Subclonal genome fraction"
);

ecdf_plot_combined <- function(variable, xlab) {
	qdraw(
		{
			plot(ecdf(tcga.pheno[, variable]), col=cols.ch["control"],
				las = 1, main = "",
				xlab = xlab,
				ylab = "Cumulative probability"
			);
			lines(ecdf(pkb.pheno[, variable]), col=cols.ch["case"]);
			legend("bottomright", legend=cohorts, pch=19, col=cols.ch, bty="n", inset=0.05);
		},
		file = insert(pdf.fname, tag=c("ecdf", "combined", variable))
	);
}

ecdf_plot_pkb <- function(variable, xlab) {
	qdraw(
		{
			plot(ecdf(pkb.pheno[, variable]), col=cols.ch["case"],
				las = 1, main = "",
				xlab = xlab,
				ylab = "Cumulative probability"
			);
			lines(ecdf(pkb.pheno[pkb.pheno$sample_type == "Primary", variable]), col=cols.st["Primary"]);
			lines(ecdf(pkb.pheno[pkb.pheno$sample_type == "Brain metastasis", variable]), col=cols.st["Brain metastasis"]);
			legend("bottomright", legend=c("Overall", "Primary", "Brain metastasis"), pch=19,
				col=c(cols.ch["case"], cols.st["Primary"], cols.st["Brain metastasis"]),
				bty="n", inset=0.05);
		},
		file = insert(pdf.fname, tag=c("ecdf", "pkb", variable))
	);
}

ecdf_plot_split <- function(variable, xlab) {
	qdraw(
		{
			plot(ecdf(tcga.pheno[, variable]), col=cols.ch["control"],
				las = 1, main = "",
				xlab = xlab,
				ylab = "Cumulative probability"
			);
			lines(ecdf(pkb.pheno[pkb.pheno$sample_type == "Primary", variable]), col=cols.st["Primary"]);
			lines(ecdf(pkb.pheno[pkb.pheno$sample_type == "Brain metastasis", variable]), col=cols.st["Brain metastasis"]);
			legend("bottomright", legend=c("TCGA", "Present, Primary", "Present, Brain metastasis"), pch=19,
				col=c(cols.ch["control"], cols.st["Primary"], cols.st["Brain metastasis"]), bty="n", inset=0.05);
		},
		file = insert(pdf.fname, tag=c("ecdf", "split", variable))
	);
}

ecdf_plot_all <- function(variable, xlab) {
	qdraw(
		{
			plot(ecdf(pkb.pheno[, variable]), col=cols.ch["case"],
				las = 1, main = "",
				xlab = xlab,
				ylab = "Cumulative probability"
			);
			lines(ecdf(tcga.pheno[, variable]), col=cols.ch["control"]);
			lines(ecdf(pkb.pheno[pkb.pheno$sample_type == "Primary", variable]), col=cols.st["Primary"]);
			lines(ecdf(pkb.pheno[pkb.pheno$sample_type == "Brain metastasis", variable]), col=cols.st["Brain metastasis"]);
			legend("bottomright", legend=c("TCGA", "Present", "Present, Primary", "Present, Brain metastasis"), pch=19,
				col=c(cols.ch["control"], cols.ch["case"], cols.st["Primary"], cols.st["Brain metastasis"]), bty="n", inset=0.05);
		},
		file = insert(pdf.fname, tag=c("ecdf", "all", variable))
	);
}

options(plot.device=NULL);

for (i in 1:length(variables)) {
	variable <- names(variables)[i];
	xlab <- variables[i];

	message(variable);

	ecdf_plot_combined(variable, xlab);
	ecdf_plot_pkb(variable, xlab);
	ecdf_plot_split(variable, xlab);
	ecdf_plot_all(variable, xlab);
}

options(plot.device=NA);


variable <- "mutation_burden_exome";
xlab <- "Exome SNVs / Mbp";

qdraw(
	{
		plot(ecdf(tcga.pheno[, variable]), col=cols.ch["control"],
			las = 1, main = "",
			xlab = xlab,
			ylab = "Cumulative probability"
		);
		lines(ecdf(pkb.pheno[pkb.pheno$specimen_material == "FF", variable]), col="green4");
		lines(ecdf(pkb.pheno[pkb.pheno$specimen_material == "FFPE", variable]), col="sienna");
		legend("bottomright", legend=c("TCGA", "Present, FF", "Present, FFPE"), pch=19,
			col=c(cols.ch["control"], "green4", "sienna"), bty="n", inset=0.05);
	},
	file = insert(pdf.fname, tag=c("ecdf", variable, "ffpe-vs-frozen"))
);

qdraw(
	{
		plot(ecdf(tcga.pheno[, variable]), col=cols.ch["control"],
			las = 1, main = "",
			xlab = xlab,
			ylab = "Cumulative probability"
		);
		lines(ecdf(pkb.pheno[pkb.pheno$wes_platform == "Agilent", variable]), col="green4");
		lines(ecdf(pkb.pheno[pkb.pheno$wes_platform == "ICE", variable]), col="sienna");
		legend("bottomright", legend=c("TCGA, Agilent", "Present, Agilent", "Present, ICE"), pch=19,
			col=c(cols.ch["control"], "green4", "sienna"), bty="n", inset=0.05);
	},
	file = insert(pdf.fname, tag=c("ecdf", variable, "agilent-vs-ice"))
);

####

burden.df <- rbind(
	transmute(tcga.pheno, cohort="TCGA", sample_id, sample_type, mutation_burden, mutation_burden_exome, mutation_burden_indel, specimen_material, covered_mb, wes_platform, center),
	transmute(pkb.pheno, cohort="Present", sample_id, sample_type, mutation_burden, mutation_burden_exome, mutation_burden_indel, specimen_material, covered_mb, wes_platform, center)
);
burden.df$cohort <- factor(burden.df$cohort, levels=cohorts);


qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=mutation_burden, fill=cohort)) +
		geom_density(alpha=0.6) + theme_bw() +
		scale_fill_manual(values=unname(cols.ch))
	,
	width = 6, height = 3,
	file = insert(pdf.fname, tag=c("density", "burden"))
)

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=mutation_burden, fill=sample_type)) +
		geom_density(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .) +
		scale_fill_manual(values=cols.st)
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "burden", "sample-type"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=covered_mb, y=mutation_burden, colour=center)) +
		geom_point(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .)
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "burden", "center"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=covered_mb, y=mutation_burden, colour=center)) +
		geom_point(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .) + scale_y_log10()
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "burden", "center", "log"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=mutation_burden_exome, fill=cohort)) +
		geom_density(alpha=0.6) + theme_bw() +
		scale_fill_manual(values=unname(cols.ch))
	,
	width = 6, height = 3,
	file = insert(pdf.fname, tag=c("density", "exome-burden"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=mutation_burden_exome, fill=sample_type)) +
		geom_density(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .) +
		scale_fill_manual(values=cols.st)
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "exome-burden", "sample-type"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=mutation_burden_exome, fill=specimen_material)) +
		geom_density(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .)
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "exome-burden", "specimen-material"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=mutation_burden_exome, fill=wes_platform)) +
		geom_density(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .)
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "exome-burden", "wes-platform"))
);

qdraw(
	ggplot(filter(burden.df, sample_type %in% c("Primary", "Brain metastasis")),
		aes(x=covered_mb, y=mutation_burden_exome, colour=center)) +
		geom_point(alpha=0.6) + theme_bw() +
		facet_grid(cohort ~ .)
	,
	width = 6, height = 5,
	file = insert(pdf.fname, tag=c("density", "exome-burden", "center"))
);

