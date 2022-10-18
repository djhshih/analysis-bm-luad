library(io);
library(GSVA);
library(ComplexHeatmap);
library(RColorBrewer);
library(magrittr);
library(limma);
library(dplyr);
library(ggplot2);
library(mmalign);
library(reshape2);

x <- qread("luad-bm_cn-expr-snv-pheno.rds");


emt.gmts <- qread("${APIARIUS_DATA}/msigdb/contrib/emt/");


absolute <- qread("../exome-seq/SLC.ABSOLUTE.table.01.20.2016.txt", type="tsv");
absolute <- absolute[match(x$pheno$pair_id, absolute$sample), ];
# TOMOVE
purity <- absolute[["Cancer DNA fraction"]];
x$pheno$purity <- purity;
x$pheno$ploidy <- absolute$ploidy;
x$pheno$genome_doublings <- absolute[["Genome doublings"]];

hallmark <- qread("${APIARIUS_DATA}/msigdb/release-51/h.all.v5.1.symbols.gmt")$data;

reactome <- qread("${APIARIUS_DATA}/msigdb/release-51/c2.cp.reactome.v5.1.symbols.gmt")$data;

oncogenic <- qread("${APIARIUS_DATA}/msigdb/release-51/c6.all.v5.1.symbols.gmt ")$data;

markers <- qread("${APIARIUS_DATA}/msigdb/contrib/brain/LEIN_BRAIN_MARKERS.gmt")$data;

myc.gmts <- qread("${APIARIUS_DATA}/msigdb/contrib/myc/");

yap.gmts <- qread("${APIARIUS_DATA}/msigdb/contrib/yap1/");

cpg <- qread("${APIARIUS_DATA}/msigdb/release-51/c2.cgp.v5.1.symbols.gmt")$data;

senes <- cpg[grep("SENESCENCE", names(cpg))];

pluri <- cpg[grep("PLURIPOTENCY|PLURIPOTENT|STEMNESS", names(cpg))];

# TOMOVE
x$expr <- x$expr[!is.na(rownames(x$expr)), ]
x$pheno$met_type <- factor(x$pheno$met_type, levels=c("Primary", "BM"));
x$pheno$patient_id <- factor(x$pheno$patient_id);


# annotate amplifications
amp <- x$cn > 6.5;

amp.genes <- c("YAP1", "MYC", "EDNRA", "RUNX1");

mmp <- rownames(x$cn)[grep("MMP", rownames(x$cn))];
mmp.amp <- apply(amp[mmp, ], 2, any);

x$amp <- t(rbind(amp[amp.genes, ], MMP=mmp.amp));

cn.median <- apply(x$cn, 2, median);


# rename samples in exp
colnames(x$expr) <- gsub("-NT-SM-.+", "", colnames(x$expr))

x$cn["MYC", x$cn["MYC",] > 6]


jaccard <- function(x, y) {
	length(intersect(x, y)) / length(union(x, y))
}

# compuete jaccard index between each x of xs against each x of xs
jaccard_matrix <- function(xs) {
	J <- matrix(
		mapply(
			jaccard,
			rep(xs, times=length(xs)),
			rep(xs, each=length(xs))
		),
		nrow = length(xs)
	);
	colnames(J) <- rownames(J) <- names(xs)
	J
}

gmts_to_genesets <- function(xs) {
	ys <- lapply(
		xs,
		function(x) {
			x$data			
		}
	);
	names(ys) <- NULL;
	unlist(ys, recursive=FALSE)
}

emt.gsets <- gmts_to_genesets(emt.gmts);
myc.gsets <- gmts_to_genesets(myc.gmts);
yap.gsets <- gmts_to_genesets(yap.gmts);


names(emt.gsets) %<>%
	gsub("EPITHELIAL_(TO_)?MESENCHYMAL_TRANSITION", "EMT", .) %>%
	gsub("_?SIGNATURE", "", .) %>%
	gsub("_?RECEPTOR_SIGNALING", "", .) %>%
	gsub("MESENCHYMAL_TRANSITION", "EMT_UP", .) %>%
	gsub("EMT_EMT", "EMT", .) %>%
	gsub("_CANCER", "", .) %>%
	gsub("TGF_BETA", "TGFB", .);

gene.sets.jaccards <- jaccard_matrix(emt.gsets);

Heatmap(gene.sets.jaccards, cluster_rows=FALSE, cluster_columns=FALSE);


g.emt <- gsva(x$expr, method = "gsva", mx.diff=TRUE, emt.gsets)$es.obs;
#g.emt <- gsva(x$expr, method = "ssgsea", mx.diff=TRUE, emt.gsets);
#g.emt <- gsva(x$expr, method = "zscore", mx.diff=TRUE, emt.gsets);

# simplify gene set names

met_type_col <- brewer.pal(12, "Set3");
names(met_type_col) <- levels(x$pheno$met_type);
met_type_col <- met_type_col[1:length(levels(x$pheno$met_type))];

ha.top <- HeatmapAnnotation(
	df = 
		#x$pheno[, c("met_type", "histology_group", "patient_id")],
		x$pheno[, c("met_type", "histology_group")],
		#x$pheno[, c("met_type", "patient_id")],
	col = list(met_type = met_type_col)
);

ha.bot <- HeatmapAnnotation(
	#df = cbind(x$mut[, -1], x$amp)
	df = x$amp
);

my.heatmap <- function(g, ...) {
	Heatmap(
		g,
		top_annotation=ha.top,
		bottom_annotation=ha.bot,
		clustering_method_columns = "average",
		clustering_method_rows = "average",
		row_names_gp = gpar(fontsize = 6),
		column_names_gp = gpar(fontsize= 6)
	);
}

my.lm <- function(g, mod) {
	fit <- lmFit(g, mod, block=x$pheno$patient_id);
	eBayes(fit)
}


my.heatmap(g.emt);

mod.met <- model.matrix(~ met_type, data=x$pheno);
eb.emt <- my.lm(g.emt, mod.met);
topTable(eb.emt);


mod.mmp <- model.matrix(~ MMP, data=as.data.frame(x$amp));
eb.emt.mmp <- my.lm(g.emt, mod.mmp);
topTable(eb.emt.mmp);

mod.myc <- model.matrix(~ MYC, data=as.data.frame(x$amp));
eb.emt.myc <- my.lm(g.emt, mod.myc);
topTable(eb.emt.myc);

mod.runx1 <- model.matrix(~ RUNX1, data=as.data.frame(x$amp));
eb.emt.runx1 <- my.lm(g.emt, mod.runx1);
topTable(eb.emt.runx1);


g.myc <- gsva(x$expr, method = "gsva", mx.diff=TRUE, myc.gsets)$es.obs;
my.heatmap(g.myc);

x$cn["MYC", ]
x$expr["MYC", ]

g.yap <- gsva(x$expr, method = "gsva", mx.diff=TRUE, yap.gsets)$es.obs;
my.heatmap(g.yap);

x$cn["YAP1", ]
x$expr["YAP1", ]


library(GSVAdata);
data(brainTxDbSets);


g.brain <- gsva(x$expr, method = "gsva", mx.diff=TRUE, brainTxDbSets)$es.obs;

my.heatmap(g.brain);

eb.brain <- my.lm(g.brain, mod.met);
topTable(eb.brain);


g.hallmark <- gsva(x$expr, method = "gsva", mx.diff=TRUE, hallmark)$es.obs;

my.heatmap(g.hallmark);

eb.hallmark <- my.lm(g.hallmark, mod.met);
topTable(eb.hallmark);

eb.hallmark.mmp <- my.lm(g.hallmark, mod.mmp);
topTable(eb.hallmark.mmp);

eb.hallmark.myc <- my.lm(g.hallmark, mod.myc);
topTable(eb.hallmark.myc);

eb.hallmark.runx1 <- my.lm(g.hallmark, mod.runx1);
topTable(eb.hallmark.runx1);


g.hallmark[grep("TRANSITION", rownames(g.hallmark)), x$mut$MMP13 == "Missense_Mutation"];

summary(g.hallmark[grep("TRANSITION", rownames(g.hallmark)), ]);


g.oncogenic <- gsva(x$expr, method = "gsva", mx.diff=TRUE, oncogenic)$es.obs;

my.heatmap(g.oncogenic);


g.markers <- gsva(x$expr, method = "gsva", mx.diff=TRUE, markers)$es.obs;

my.heatmap(g.markers);
eb.markers <- my.lm(g.markers, mod.met);
topTable(eb.markers);


reactome.up.sel <- c(
	"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
	"REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
	"REACTOME_COLLAGEN_FORMATION",
	"REACTOME_CELL_CELL_JUNCTION_ORGANIZATION",
	"REACTOME_SIGNALING_BY_HIPPO",
	"REACTOME_GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS"
);

g.reactome.up <- gsva(x$expr, method = "gsva", mx.diff=TRUE, reactome[reactome.up.sel])$es.obs;
my.heatmap(g.reactome.up);

reactome.down.sel <- c(
	"REACTOME_CELL_CYCLE",
	"REACTOME_G1_S_TRANSITION",
	"REACTOME_G2_M_CHECKPOINTS",
	"REACTOME_DNA_REPLICATION",
	"REACTOME_METABOLISM_OF_RNA",
	"REACTOME_TRANSCRIPTION",
	"REACTOME_TRANSLATION"
	"REACTOME_PEPTIDE_CHAIN_ELONGATION",
	"REACTOME_GLYCOLYSIS",
	"REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"
);

g.reactome.down <- gsva(x$expr, method = "gsva", mx.diff=TRUE, reactome[reactome.down.sel])$es.obs;
my.heatmap(g.reactome.down);

markers.brain.cells <- markers[c("LEIN_ASTROCYTE_MARKERS", "LEIN_OLIGODENDROCYTE_MARKERS", "LEIN_NEURON_MARKERS")];

g.reactome.down.and.markers <- gsva(x$expr, method="gsv", mx.diff=TRUE, c(markers.brain.cells, reactome[reactome.down.sel]))$es.obs;
my.heatmap(g.reactome.down.and.markers);

g.senes <- gsva(x$expr, method = "gsva", mx.diff=TRUE, senes)$es.obs;
my.heatmap(g.senes);

g.pluri <- gsva(x$expr, method = "gsva", mx.diff=TRUE, pluri)$es.obs;
my.heatmap(g.pluri);

g.reactome.down.and.pluri <- gsva(x$expr, method = "gsva", mx.diff=TRUE, c(pluri, reactome[reactome.down.sel]))$es.obs;
my.heatmap(g.reactome.down.and.pluri);

plot(g.markers["LEIN_OLIGODENDROCYTE_MARKERS", ], g.reactome.down["REACTOME_CELL_CYCLE", ])
cor(g.markers["LEIN_OLIGODENDROCYTE_MARKERS", ], g.reactome.down["REACTOME_CELL_CYCLE", ])

# astrocyte/glial/neural stem cell marker
cor(x$expr["GFAP", ], g.reactome.down["REACTOME_CELL_CYCLE", ])

# neuronal progenitor marker
cor(x$expr["NES", ], g.reactome.down["REACTOME_CELL_CYCLE", ])

# oligodendrocyte markers
cor(x$expr["CNP", ], g.reactome.down["REACTOME_CELL_CYCLE", ])
cor(x$expr["PMP2", ], g.reactome.down["REACTOME_CELL_CYCLE", ])

plot(x$expr["CNP", ], x$expr["PMP2", ])

cor(absolute[["Cancer DNA fraction"]], g.reactome.down["REACTOME_CELL_CYCLE", ])
plot(absolute[["Cancer DNA fraction"]], g.reactome.down["REACTOME_CELL_CYCLE", ])


gene.sets.jaccards <- jaccard_matrix(reactome[reactome.down.sel]);
Heatmap(gene.sets.jaccards);

hist(g.reactome.down["REACTOME_CELL_CYCLE", ], breaks=6)

reactome.quiescence.sel <- c(
	"REACTOME_CELL_CYCLE",
	"REACTOME_GLYCOLYSIS",
	"REACTOME_DNA_REPLICATION",
	"REACTOME_TRANSCRIPTION",
	"REACTOME_TRANSLATION"
);

g.reactome.quiescence <- g.reactome.down[reactome.quiescence.sel, ];

x$gsva <- list();
x$gsva$quiescence <- g.reactome.quiescence;


quiescence.score <- -apply(g.reactome.quiescence, 2, mean);
x$pheno$quiescence <- quiescence.score;


hist(quiescence.score, breaks=10)

cor(absolute[["Cancer DNA fraction"]], quiescence.score)
plot(absolute[["Cancer DNA fraction"]], quiescence.score)
cor(x$expr["GFAP", ], g.reactome.down["REACTOME_CELL_CYCLE", ])

#t.cell.markers <- c("TRA", "TRB", "CD3G", "CD3D", "CD3E");
t.cell.markers <- c("CD3G", "CD3D", "CD3E");

match(t.cell.markers, rownames(x$expr))
plot(x$expr["CD3G", ], x$expr["CD3D", ]);
plot(x$expr["CD3G", ], x$expr["CD3E", ]);
plot(x$expr["CD3D", ], x$expr["CD3E", ]);

t.cells <- apply(x$expr[t.cell.markers, ], 2, max);
hist(t.cells, breaks=10);

cor(t.cells, quiescence.score);
plot(t.cells, quiescence.score);

cor(x$expr["GFAP", ], quiescence.score)
cor(x$expr["GFAP", ], quiescence.score, method="kendall")
plot(x$expr["GFAP", ], quiescence.score)

# neuronal progenitor marker
# bad marker, not restricted to CNS
cor(x$expr["NES", ], quiescence.score)
# MAP2 can be probably expressed by neuroendocrine cells
cor(x$expr["MAP2", ], quiescence.score)
plot(x$expr["MAP2", ], quiescence.score)

# lack of markers of differentiated neurons suggest that neurons are largely
# absent
# alternatively, these genes could be more easily degraded

# NeuN (RBFOX3)
#plot(x$expr["RBFOX3", ], quiescence.score)
# neuronal-specific genes do not appear to be expressed highly

# Beta III tubulin (TUBB3)
#cor(x$expr["TUBB3", ], quiescence.score)



# oligodendrocyte markers
cor(x$expr["CNP", ], quiescence.score)
cor(x$expr["CNP", ], quiescence.score, method="kendall")
plot(x$expr["CNP", ], quiescence.score)

cor(x$expr["PMP2", ], quiescence.score)
cor(x$expr["PMP2", ], quiescence.score, method="kendall")

wilcox.test(quiescence.score ~ met_type, data = x$pheno);

wilcox.test(quiescence.score ~ YAP1, data = cbind(x$pheno, x$amp));
wilcox.test(quiescence.score ~ MYC, data = cbind(x$pheno, x$amp));
wilcox.test(quiescence.score ~ RUNX1, data = cbind(x$pheno, x$amp));
plot(quiescence.score ~ RUNX1, data = cbind(x$pheno, x$amp))

plot(quiescence.score ~ RUNX1, data = cbind(x$pheno, t(x$expr)))
cor(x$expr["RUNX1", ], quiescence.score);
cor(x$expr["RUNX1", ], quiescence.score, method="kendall");

plot(x$expr["RUNX1", ], quiescence.score);
plot(x$expr["RUNX2", ], quiescence.score);
plot(x$expr["RUNX3", ], quiescence.score);

fit <- lm(quiescence.score ~ met_type + GFAP + NES + CNP + t.cells + YAP1 + MYC + RUNX1, data = cbind(x$pheno, t(x$expr[c("GFAP", "NES", "CNP"), ]), t.cells, x$amp));
summary(fit)

fit <- lm(quiescence.score ~ met_type + GFAP + RUNX1, data = cbind(x$pheno, t(x$expr[c("GFAP", "NES", "CNP", "RUNX1"), ]), t.cells));
summary(fit)

fit <- lm(quiescence.score ~ GFAP + RUNX1, data = cbind(x$pheno, t(x$expr[c("GFAP", "NES", "CNP", "RUNX1"), ]), t.cells));
summary(fit)

# CD45 (a glycosylated form of CD45, B220, is a marker of B cells)
plot(x$expr["PTPRC", ], t.cells);
plot(x$expr["PTPRC", ], quiescence.score);
plot(x$expr["PTPRC", ], x$expr["RUNX1", ]);

# leukocyte marker (CD18)
plot(x$expr["ITGB2", ], quiescence.score);
cor(x$expr["ITGB2", ], quiescence.score);
cor(x$expr["ITGB2", ], quiescence.score, method="kendall");

mod.q <- model.matrix(~ quiescence.score);
fit <- lmFit(x$expr, mod.q);
eb.quiescence <- eBayes(fit);
tt.q <- topTable(eb.quiescence, number=Inf);
tt.q[tt.q$logFC > 0, ][1:10, ]
tt.q[tt.q$logFC < 0, ][1:10, ]

plot(x$expr["MYC", ], quiescence.score);

# top genes correlated with quiescence.score in current data
plot(x$expr["KPNA2", ], quiescence.score);
unlist(lapply(reactome[reactome.quiescence.sel], function(x) any(x == "KPNA2")))

# Downregulation of KPNA2 in non-small-cell lung cancer is associated with
# Oct4 (POU5F1) expression
# Li 2013

#plot(x$expr["KPNA2", ], x$expr["POU5F1", ]);


plot(x$expr["C1orf116", ], quiescence.score);
plot(x$expr["RHOU", ], quiescence.score);
plot(x$expr["TMEM63A", ], quiescence.score);
plot(x$expr["LAMB2", ], quiescence.score);

# quiescence associated transcription factors
# Foxp1, Foxo1, and Klf2
# Wong. Mol Immunol. 2015;68:223-33
# http://www.ncbi.nlm.nih.gov/pubmed/26350416
tt.q["FOXP1", ]
tt.q["FOXO1", ]
tt.q["KLF2", ]

tt.q["KLF4", ]

plot(x$expr["RUNX1", ], x$expr["FOXP1", ])
plot(x$expr["RUNX1", ], x$expr["FOXO1", ])
plot(x$expr["RUNX1", ], x$expr["KLF2", ])
plot(x$expr["KLF2", ], quiescence.score);
plot(x$expr["RUNX1", ], quiescence.score);


# leukocyte marker
# CD45 (PTPRC)

# leukocyte marker (CD18)
# ITGB2

# activated B cells and monocytes
# CD80

# macrophage and microglia marker
# ITGAM (CD11B)

plot(x$expr["PTPRC", ], x$expr["ITGB2", ]);
abline(a=0, b=1);
# PTPRC appears more highly expressed

plot(x$expr["PTPRC", ], x$expr["ITGAM", ]);
plot(x$expr["GFAP", ], x$expr["NES", ]);
plot(x$expr["GFAP", ], x$expr["CNP", ]);

plot(x$expr["PMP2", ], x$expr["CNP", ]);
abline(a=0, b=1);

markers <- c("GFAP", "NES", "CNP", "CD3G", "PTPRC", "ITGAM");

normal.markers <- x$expr[markers, ];


fit <- lm(quiescence.score ~ met_type + met_type : purity, data=x$pheno);
summary(fit);


qwrite(x, filename("luad-bm", tag="cn-expr-snv-pheno-gsva", ext="rds"));



qc <- qread("Brastianos_BMET_RNAQC.txt", type="tsv");
qc$rnaseq_id <- qc$"Collaborator Sample ID";
qc <- qc[match(x$pheno$rnaseq_id, qc$rnaseq_id), ];
x$pheno$dv200 <- qc$DV200;


p <- pca(x$expr);
pca_plot(p, x$pheno, aes(colour=quiescence)) + 
	scale_colour_gradient2(low="green", mid="black", high="red");

ks <- 1:nrow(p$Z);

dv200.cor <- unlist(lapply(ks, function(k) cor(p$Z[k, ], x$pheno$dv200)));
dv200.cor.p <- unlist(lapply(ks, function(k) cor.test(p$Z[k, ], x$pheno$dv200)$p.value));
g <- ggplot(, aes(x=ks, y=dv200.cor, fill=-log10(dv200.cor.p))) + geom_bar(stat="identity");
g

quies.cor <- unlist(lapply(ks, function(k) cor(p$Z[k, ], x$pheno$quiescence)));
quies.cor.p <- unlist(lapply(ks, function(k) cor.test(p$Z[k, ], x$pheno$quiescence)$p.value));
g <- ggplot(, aes(x=ks, y=quies.cor, fill=-log10(quies.cor.p))) + geom_bar(stat="identity");
g

pca_plot(p, x$pheno, aes(colour=purity))
purity.cor <- unlist(lapply(ks, function(k) cor(p$Z[k, ], x$pheno$purity)));
purity.cor.p <- unlist(lapply(ks, function(k) cor.test(p$Z[k, ], x$pheno$purity)$p.value));
g <- ggplot(, aes(x=ks, y=purity.cor, fill=-log10(purity.cor.p))) + geom_bar(stat="identity");
g


fit1 <- lm(p$Z[1,] ~ purity + dv200, data=x$pheno);
summary(fit1);

fit2 <- lm(p$Z[1,] ~  purity + dv200 + quiescence, data=x$pheno);
summary(fit2);

anova(fit1, fit2, test="Chisq");

expr <- x$expr;
colnames(expr) <- x$pheno$pair_id;
expr.m <- melt(expr, varnames=c("gene", "pair_id"));
expr.m$quiescence <- x$pheno$quiescence[match(expr.m$pair_id, x$pheno$pair_id)];
expr.m$quiescent <- expr.m$quiescence > 0;
g <- ggplot(expr.m, aes(x=value, fill=quiescent)) + geom_density(alpha=0.2);
g


fit <- lm(p$Z[2,] ~ quiescence + purity + dv200, data=x$pheno);
summary(fit);

fit <- lm(p$Z[3,] ~ quiescence + purity + dv200, data=x$pheno);
summary(fit);

plot(quiescence ~ dv200, data=x$pheno);
with(x$phen, cor(quiescence, dv200));
with(x$phen, cor(quiescence, dv200, method="kendall"));

plot(quiescence ~ x$expr["MYC", ], data=x$pheno);
plot(x$expr["MYC", ] ~ dv200, data=x$pheno);
with(x$phen, cor(quiescence, x$expr["MYC", ]));
with(x$phen, cor(quiescence, x$expr["MYC", ], method="kendall"));

summary(lm(quiescence ~ dv200 + x$expr["MYC", ], data=x$pheno));
fit <- lm(quiescence ~ dv200, data=x$pheno);
summary(fit);
quiescence.c <-  x$pheno$quiescence - coef(fit)[2] * x$pheno$dv200;
summary(lm(quiescence.c ~ x$expr["MYC", ]))

plot(quiescence ~ x$expr["MYC", ], data=x$pheno);
plot(quiescence.c ~ x$expr["MYC", ]);
cor(quiescence.c, x$expr["MYC", ]);

plot(quiescence ~ purity, data=x$pheno);
plot(purity ~ dv200, data=x$pheno);


