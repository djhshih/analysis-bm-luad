library(io);
library(dplyr);
library(ComplexHeatmap);
library(RColorBrewer);
library(circlize);

options(stringsAsFactors=FALSE);

clin <- qread("~/exigens/tcga/tcga-luad/annot/patient-info_tcga-luad_stage2.rds");
clin <- rename(clin, age_at_primary_dx = age_at_initial_pathologic_diagnosis);

#pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv");

out.fname <- filename("tcga-luad", date=NA);

cn <- qread("tcga-luad_pass_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);

to_tcga_patient_id <- function(x) {
	unlist(lapply(strsplit(x, "-"), function(x) paste(x[2:4], collapse="-")))
}

if (colnames(cn) != clin$clinical_id) {
	samples <- colnames(cn);
	patients <- to_tcga_patient_id(samples);
	colnames(cn) <- patients;
	clin <- clin[match(patients, clin$clinical_id), ];
} else {
	patients <- colnames(cn);
}


mut.orig <- qread("tcga-luad_pass_black-f_20171113T105434.pset.mmaf", type="tsv", quote="\"");


mut <- mutate(mut.orig, sample_id = Tumor_Sample_Barcode, clinical_id = to_tcga_patient_id(sample_id));


#gene.set <- "PKB";
#gene.set <- "TCGA";
gene.set <- "PKB2"

if (gene.set == "PKB") {

	mut.genes <- as.character(qread("../pkb-luad_mutsigcv_sig-genes.vtr"));
	cn.genes <- as.character(qread("../pkb-luad_gistic_sig-genes.vtr"));
	comut.fname <- insert(out.fname, c("comut", "pkb-luad-genes"), ext="pdf");

	genes <- c(mut.genes, setdiff(cn.genes, mut.genes));
	gene.types <- factor(character(length(genes)), levels=c("MutSig2CV", "GISTIC"));
	gene.types[genes %in% cn.genes] <- "GISTIC";
	gene.types[genes %in% mut.genes] <- "MutSig2CV";
	gene.types[genes %in% c("MMP13", "CDKN2A")] <- "GISTIC";

} else if (gene.set == "TCGA") {

	mut.genes <- as.character(qread("tcga-luad_mutsigcv_sig-genes.vtr"));
	cn.genes <- as.character(qread("tcga-luad_gistic_sig-genes.vtr"));
	comut.fname <- insert(out.fname, c("comut", "tcga-luad-genes"), ext="pdf");

	genes <- c(mut.genes, setdiff(cn.genes, mut.genes));
	gene.types <- factor(character(length(genes)), levels=c("TCGA MutSig2CV", "TCGA GISTIC"));
	gene.types[genes %in% cn.genes] <- "TCGA GISTIC";
	gene.types[genes %in% mut.genes] <- "TCGA MutSig2CV";
	gene.types[genes %in% c("CDKN2A")] <- "TCGA GISTIC";

} else if (gene.set == "PKB2") {

	mut.genes <- c("TP53", "KRAS", "STK11", "KEAP1", "EGFR");
	cn.genes <- as.character(qread("../pkb-luad_gistic_sig-genes.vtr"));
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

# check all genes are present exactly once in the CN matrix
stopifnot(genes %in% rownames(cn))
stopifnot(sum(rownames(cn) %in% genes) == length(genes))


mut.sel <- filter(mut, Hugo_Symbol %in% genes);
mut.sel$gene <- factor(mut.sel$Hugo_Symbol, genes);

# missense
# nonframeshift indel
# nonsense, splicing, frameshift indel
# amplification
# homozygous deletion


cn.sel <- cn[genes, ];

G <- length(genes);
N <- length(patients);

get_mutations <- function(x, gene, mutation_type) {
	as.integer(patients %in% as.character(
		x$clinical_id[x$Hugo_Symbol == gene &
			x$Variant_Classification %in% mutation_type]
	))
}

get_mutations(mut.sel, "KRAS", "Missense_Mutation")

get_mutation_matrix <- function(mut.sel, mutation_type) {
	t(matrix(
		unlist(lapply(genes,
			function(gene) {
				get_mutations(mut.sel, gene, mutation_type)
			}
		)),
		nrow=N, ncol=G, dimnames=list(patient=patients, gene=genes)));
}

table(mut.sel$Variant_Classification)

missense <- get_mutation_matrix(mut.sel, "Missense_Mutation");
nfsindel <- get_mutation_matrix(mut.sel, c("In_Frame_Del", "In_Frame_Ins"));
disrupt <- get_mutation_matrix(mut.sel, c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"));

rowSums(missense)
rowSums(nfsindel)
rowSums(disrupt)

get_cna <- function(x, gene, cna_type) {
	if (cna_type == "amp") {
		idx <- x[gene, ] > 8;
	} else if (cna_type == "homodel") {
		idx <- x[gene, ] < 0.5;
	} else {
		stop("CNA type not supported")
	}
	as.integer(idx)
}

get_cna_matrix <- function(cn, cna_type) {
	t(matrix(
		unlist(lapply(genes,
			function(gene) {
				get_cna(cn, gene, cna_type)
			}
		)),
		nrow=N, ncol=G, dimnames=list(patient=patients, gene=genes)));
}

amp <- get_cna_matrix(cn.sel, "amp");
del <- get_cna_matrix(cn.sel, "homodel");

rowSums(amp)

mats <- list(
	missense = missense,
	nfsindel = nfsindel,
	disrupt = disrupt,
	amp = amp,
	del = del
);

col <- c(
	background = "#CCCCCC",
	missense="#009000",
	nfsindel="orange3",
	disrupt="purple",
	amp="red2",
	del="royalblue2"
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
	#row_order = NULL,
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

#pdf(file = tag(comut.fname, "wide"), width=32, height=10);
pdf(file = tag(comut.fname), width=12, height=10);
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom");
dev.off();

