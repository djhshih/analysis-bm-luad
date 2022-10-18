library(io);
library(dplyr);
library(ComplexHeatmap);

options(stringsAsFactors=FALSE);

pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2_pass_luad.tsv");
mut.genes <- as.character(qread("pkb-luad_mutsigcv_sig-genes.vtr"));
cn.genes <- as.character(qread("pkb-luad_gistic_sig-genes.vtr"));

out.fname <- filename("pkb-luad", date=NA);

cn <- qread("brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv");
cn <- as.matrix(cn);


mut.prims <- qread("pkb-luad_prims_black-f_20171031T113017.pset.mmaf", type="tsv", quote="\"");
mut.mets <- qread("pkb-luad_brain-mets_black-f_20171031T112959.pset.mmaf", type="tsv", quote="\"");

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

genes <- union(mut.genes, cn.genes);

# check all genes are present exactly once in the CN matrix
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

patients <- as.character(unique(pheno$clinical_id));

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

missense_prim <- get_mutation_matrix(mut.sel, "Missense_Mutation", "Primary");
missense_bmet <- get_mutation_matrix(mut.sel, "Missense_Mutation", "Brain metastasis");
nfsindel_prim <- get_mutation_matrix(mut.sel, c("In_Frame_Del", "In_Frame_Ins"), "Primary");
nfsindel_bmet <- get_mutation_matrix(mut.sel, c("In_Frame_Del", "In_Frame_Ins"), "Brain metastasis");
disrupt_prim <- get_mutation_matrix(mut.sel, c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"), "Primary");
disrupt_bmet <- get_mutation_matrix(mut.sel, c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"), "Brain metastasis");

rowSums(missense_prim)
rowSums(missense_bmet)
rowSums(nfsindel_prim)
rowSums(nfsindel_bmet)
rowSums(disrupt_prim)
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

amp_prim <- get_cna_matrix(cn.sel, "amp", "Primary");
amp_bmet <- get_cna_matrix(cn.sel, "amp", "Brain metastasis");
del_prim <- get_cna_matrix(cn.sel, "homodel", "Primary");
del_bmet <- get_cna_matrix(cn.sel, "homodel", "Brain metastasis");


mats <- list(
	missense_prim = missense_prim, missense_bmet = missense_bmet,
	nfsindel_prim = nfsindel_prim, nfsindel_bmet = nfsindel_bmet,
	disrupt_prim = disrupt_prim, disrupt_bmet = disrupt_bmet,
	amp_prim = amp_prim, amp_bmet = amp_bmet,
	del_prim = del_prim, del_bmet = del_bmet
);

col <- c(
	background = "#CCCCCC",
	missense_prim="#009000", missense_bmet="#006000",
	nfsindel_prim="orange3", nfsindel_bmet="orange4",
	disrupt_prim="purple", disrupt_bmet="purple4",
	amp_prim="red2", amp_bmet="red4",
	del_prim="royalblue2", del_bmet="royalblue4"
);

pdf(file = tag(out.fname, c("comut", "dual"), ext="pdf"), width=8, height=10);

oncoPrint(mats,
	alter_fun = list(
			background = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["background"], col=NA)),
			amp_prim = function(x, y, w, h) grid.rect(x, y+h*0.2, w*0.8, h*0.3, gp = gpar(fill = col["amp_prim"], col = NA)),
			amp_bmet = function(x, y, w, h) grid.rect(x, y-h*0.2, w*0.8, h*0.3, gp = gpar(fill = col["amp_bmet"], col = NA)),
			del_prim = function(x, y, w, h) grid.rect(x, y+h*0.2, w*0.8, h*0.3, gp = gpar(fill = col["del_prim"], col = NA)),
			del_bmet = function(x, y, w, h) grid.rect(x, y-h*0.2, w*0.8, h*0.3, gp = gpar(fill = col["del_bmet"], col = NA)),
			missense_prim = function(x, y, w, h) grid.rect(x, y+h*0.2, w*0.8, h*0.2, gp = gpar(fill = col["missense_prim"], col = NA)),
			missense_bmet = function(x, y, w, h) grid.rect(x, y-h*0.2, w*0.8, h*0.2, gp = gpar(fill = col["missense_bmet"], col = NA)),
			nfsindel_prim = function(x, y, w, h) grid.rect(x, y+h*0.2, w*0.8, h*0.2, gp = gpar(fill = col["nfsindel_prim"], col = NA)),
			nfsindel_bmet = function(x, y, w, h) grid.rect(x, y-h*0.2, w*0.8, h*0.2, gp = gpar(fill = col["nfsindel_bmet"], col = NA)),
			disrupt_prim = function(x, y, w, h) grid.rect(x, y+h*0.2, w*0.8, h*0.2, gp = gpar(fill = col["disrupt_prim"], col = NA)),
			disrupt_bmet = function(x, y, w, h) grid.rect(x, y-h*0.2, w*0.8, h*0.2, gp = gpar(fill = col["disrupt_bmet"], col = NA))
	),
	row_order = NULL,
	col = col
);

dev.off();



# test for associations

sum(get_cna(cn.sel, "MYC", "amp", "Primary"))
sum(get_cna(cn.sel, "MYC", "amp", "Brain metastasis"))

ct.myc <- table(data.frame(
	prim = get_cna(cn.sel, "MYC", "amp", "Primary"),
	bmet = get_cna(cn.sel, "MYC", "amp", "Brain metastasis")
));

sum(get_cna(cn.sel, "YAP1", "amp", "Primary"))
sum(get_cna(cn.sel, "YAP1", "amp", "Brain metastasis"))

sum(get_cna(cn.sel, "YAP1", "amp", c("Primary", "Brain metastasis")))

ct.yap1 <- table(data.frame(
	prim = get_cna(cn.sel, "YAP1", "amp", "Primary"),
	bmet = get_cna(cn.sel, "YAP1", "amp", "Brain metastasis")
));

ct.yap1
fisher.test(ct.yap1);
mcnemar.test(ct.yap1);

ct.mmp13 <- table(data.frame(
	prim = get_cna(cn.sel, "MMP13", "amp", "Primary"),
	bmet = get_cna(cn.sel, "MMP13", "amp", "Brain metastasis")
));

ct.mmp13
mcnemar.test(ct.mmp13);


variants <- unique(mut.sel$Variant_Classification);
ct.stk11 <- table(data.frame(
	prim = get_mutations(mut.sel, "STK11", variants, "Primary"),
	bmet = get_mutations(mut.sel, "STK11", variants, "Brain metastasis")
));

ct.stk11
mcnemar.test(ct.stk11);

midp.binom.test(ct.myc[2,1], ct.myc[2,1] + ct.myc[1,2], prob=0.5, alternative="less")
midp.binom.test(ct.yap1[2,1], ct.yap1[2,1] + ct.yap1[1,2], prob=0.5, alternative="less")
midp.binom.test(ct.mmp13[2,1], ct.mmp13[2,1] + ct.mmp13[1,2], prob=0.5, alternative="less")
midp.binom.test(ct.stk11[2,1], ct.stk11[2,1] + ct.stk11[1,2], prob=0.5, alternative="less")

