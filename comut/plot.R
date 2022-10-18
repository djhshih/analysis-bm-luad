library(io);
library(dplyr);
library(ComplexHeatmap);
library(RColorBrewer);
library(circlize);  # colorRamp2

options(stringsAsFactors=FALSE);

clin <- qread("~/exigens/brain-mets/annot/patient-info_stage2.rds");
pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv");
patients <- as.character(unique(pheno$clinical_id));
clin <- clin[match(patients, clin$clinical_id), ];

cn <- qread("brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);

mut.prims <- qread("pkb-luad_prims_black-f.pset.mmaf", type="tsv", quote="\"");
mut.mets <- qread("pkb-luad_brain-mets_black-f.pset.mmaf", type="tsv", quote="\"");

out.fname <- filename("pkb-luad", date=NA);

#gene.set <- "PKB";
#gene.set <- "TCGA";
gene.set <- "PKB2"

if (gene.set == "PKB") {

	mut.genes <- as.character(qread("pkb-luad_mutsigcv_sig-genes.vtr"));
	cn.genes <- as.character(qread("pkb-luad_gistic_sig-genes.vtr"));
	comut.fname <- insert(out.fname, c("comut", "pkb-luad-genes"), ext="pdf");

	genes <- c(mut.genes, setdiff(cn.genes, mut.genes));
	gene.types <- factor(character(length(genes)), levels=c("MutSig2CV", "GISTIC"));
	gene.types[genes %in% cn.genes] <- "GISTIC";
	gene.types[genes %in% mut.genes] <- "MutSig2CV";
	gene.types[genes %in% c("MMP13", "CDKN2A")] <- "GISTIC";

} else if (gene.set == "TCGA") {

	mut.genes <- as.character(qread("tcga-luad/tcga-luad_mutsigcv_sig-genes.vtr"));
	cn.genes <- as.character(qread("tcga-luad/tcga-luad_gistic_sig-genes.vtr"));
	comut.fname <- insert(out.fname, c("comut", "tcga-luad-genes"), ext="pdf");

	genes <- c(mut.genes, setdiff(cn.genes, mut.genes));
	gene.types <- factor(character(length(genes)), levels=c("TCGA MutSig2CV", "TCGA GISTIC"));
	gene.types[genes %in% cn.genes] <- "TCGA GISTIC";
	gene.types[genes %in% mut.genes] <- "TCGA MutSig2CV";
	gene.types[genes %in% c("CDKN2A")] <- "TCGA GISTIC";

} else if (gene.set == "PKB2") {

	mut.genes <- c("TP53", "KRAS", "STK11", "KEAP1", "EGFR");
	cn.genes <- as.character(qread("pkb-luad_gistic_sig-genes.vtr"));
	crit.genes <- c("YAP1", "MMP13", "MYC", "CDKN2A", "CDKN2B");
	cn.genes <- c(crit.genes, setdiff(cn.genes, crit.genes));

	comut.fname <- insert(out.fname, c("comut", "pkb-luad-genes"), ext="pdf");

	genes <- c(cn.genes, setdiff(mut.genes, cn.genes));
	# move critical genes to the front

	gene.types <- factor(character(length(genes)), levels=c("MutSig2CV", "GISTIC"));
	gene.types[genes %in% cn.genes] <- "GISTIC";
	gene.types[genes %in% mut.genes] <- "MutSig2CV";
	gene.types[genes %in% c("CDKN2A")] <- "GISTIC";

}


ith <- function(xs, i) {
	unlist(lapply(xs, function(x) x[i]))
}

mut <- rbind(mut.prims, mut.mets) %>% 
	mutate(
		fh_sample_id = ith(strsplit(as.character(Tumor_Sample_Barcode), ":"), 1),
		sample_index = match(fh_sample_id, pheno$fh_sample_id),
		sample_id = pheno$sample_id[sample_index],
		sample_type = pheno$sample_type[sample_index],
		clinical_id = pheno$clinical_id[sample_index]
	);


# check all genes are present exactly once in the CN matrix
genes[which(! genes %in% rownames(cn))]
stopifnot(genes %in% rownames(cn))
stopifnot(sum(rownames(cn) %in% genes) == length(genes))

cn.sel <- cn[genes, ];



# mid-p binomial test
midp.binom.test <- function(x, n, prob, alternative = c("two.sided", "less", "greater")) {
	r <- x / n;
	p <- pbinom(x, n, prob=prob) - 0.5 * dbinom(x, n, prob=prob);
	alternative <- match.arg(alternative);
	if (alternative == "two.sided") {
		p <- p * 2;
	} else if (alternative == "greater") {
		p <- 1 - p;
	}
	list(p.value = p, statistic = r)
}


mut.sel <- filter(mut, Hugo_Symbol %in% genes);
mut.sel$gene <- factor(mut.sel$Hugo_Symbol, genes);

# missense
# nonframeshift indel
# nonsense, splicing, frameshift indel
# amplification
# homozygous deletion


G <- length(genes);
N <- length(patients);

get_mutations <- function(x, gene, mutation_type, sample_type) {
	as.integer(patients %in% as.character(
		x$clinical_id[x$Hugo_Symbol == gene &
			x$Variant_Classification %in% mutation_type &
			x$sample_type %in% sample_type]
	))
}

get_mutations(mut.sel, "KRAS", "Missense_Mutation", "Primary")
get_mutations(mut.sel, "KRAS", "Missense_Mutation", "Brain metastasis")

get_mutation_matrix <- function(mut.sel, mutation_type, sample_type) {
	t(matrix(
		unlist(lapply(genes,
			function(gene) {
				get_mutations(mut.sel, gene, mutation_type, sample_type)
			}
		)),
		nrow=N, ncol=G, dimnames=list(patient=patients, gene=genes)));
}

table(mut.sel$Variant_Classification)
table(mut.sel$sample_type)
table(mut.sel$gene, mut.sel$sample_type)

#missense_prim <- get_mutation_matrix(mut.sel, "Missense_Mutation", "Primary");
missense_bmet <- get_mutation_matrix(mut.sel, "Missense_Mutation", "Brain metastasis");
#nfsindel_prim <- get_mutation_matrix(mut.sel, c("In_Frame_Del", "In_Frame_Ins"), "Primary");
nfsindel_bmet <- get_mutation_matrix(mut.sel, c("In_Frame_Del", "In_Frame_Ins"), "Brain metastasis");
#disrupt_prim <- get_mutation_matrix(mut.sel, c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"), "Primary");
disrupt_bmet <- get_mutation_matrix(mut.sel, c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"), "Brain metastasis");

rowSums(missense_bmet)
rowSums(nfsindel_bmet)
rowSums(disrupt_bmet)

get_cna <- function(x, gene, cna_type, sample_type) {
	samples <- pheno$sample_id[pheno$sample_type %in% sample_type];
	if (cna_type == "amp") {
		idx <- x[gene, match(samples, colnames(x))] > 8;
	} else if (cna_type == "homodel") {
		idx <- x[gene, match(samples, colnames(x))] < 0.5;
	} else {
		stop("CNA type not supported")
	}
	as.integer(patients %in% pheno$clinical_id[match(samples[idx], pheno$sample_id)])
}

get_cna_matrix <- function(cn, cna_type, sample_type) {
	t(matrix(
		unlist(lapply(genes,
			function(gene) {
				get_cna(cn, gene, cna_type, sample_type)
			}
		)),
		nrow=N, ncol=G, dimnames=list(patient=patients, gene=genes)));
}

#amp_prim <- get_cna_matrix(cn.sel, "amp", "Primary");
amp_bmet <- get_cna_matrix(cn.sel, "amp", "Brain metastasis");
#del_prim <- get_cna_matrix(cn.sel, "homodel", "Primary");
del_bmet <- get_cna_matrix(cn.sel, "homodel", "Brain metastasis");


mats <- list(
	missense = missense_bmet,
	nfsindel = nfsindel_bmet,
	disrupt = disrupt_bmet,
	amp = amp_bmet,
	del = del_bmet
);

col <- c(
	background = "#CCCCCC",
	missense = "#009000",
	nfsindel = "orange3",
	disrupt = "purple",
	amp = "red2",
	del = "royalblue2"
);

exac.pop.cols <- brewer.pal(7, "Set1");
names(exac.pop.cols) <- c("FIN", "EAS", "AFR", "AMR", "NFE", "SAS", "OTH");

ha <- HeatmapAnnotation(
	df = as.data.frame(select(clin, sex, age=age_at_primary_dx, smoking=smoking_signature_high, exac_pop)),
	col = list(
		sex = c(XX="plum", XY="plum4"),
		age = colorRamp2(c(30, 90), c("lightpink", "red3")), 
		smoking = c(Low="goldenrod", High="goldenrod4"),
		exac_pop = exac.pop.cols
	),
	gap = unit(0.5, "mm"),
	show_annotation_name = TRUE
);

ht <- oncoPrint(mats,
	alter_fun = list(
			background = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["background"], col=NA)),
			amp = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["amp"], col = NA)),
			del = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["del"], col = NA)),
			missense = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.4, gp = gpar(fill = col["missense"], col = NA)),
			nfsindel = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.4, gp = gpar(fill = col["nfsindel"], col = NA)),
			disrupt = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.4, gp = gpar(fill = col["disrupt"], col = NA))
	),
	row_order = NULL,
	heatmap_legend_param = list(
		title = "Alterations",
		at = c("missense", "nfsindel", "disrupt", "amp", "del"),
		labels = c("Missense mutation", "In-frame indel", "Nonsense/frameshift/splice", "High amplification", "Homozygous deletion"),
		nrow = 2
	),
	top_annotation = ha,
	show_row_barplot = FALSE,
	split = gene.types,
	col = col
);

#pdf(file = tag(comut.fname), width=6, height=10);
pdf(file = tag(comut.fname), width=15, height=10);
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom");
dev.off();

