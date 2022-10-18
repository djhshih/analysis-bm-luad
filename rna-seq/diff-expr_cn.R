library(io);
library(GSVA);
library(ComplexHeatmap);
library(RColorBrewer);
library(magrittr);
library(limma);


x <- qread("luad-bm_cn-expr-snv-pheno.rds");

emt.gmts <- qread("${APIARIUS_DATA}/msigdb/contrib/emt/");

hallmark <- qread("${APIARIUS_DATA}/msigdb/release-51/h.all.v5.1.symbols.gmt")$data;

markers <- qread("${APIARIUS_DATA}/msigdb/contrib/brain/LEIN_BRAIN_MARKERS.gmt")$data;

myc.gmts <- qread("${APIARIUS_DATA}/msigdb/contrib/myc/");

absolute <- qread("../exome-seq/SLC.ABSOLUTE.table.01.20.2016.txt", type="tsv");
absolute <- absolute[match(x$pheno$pair_id, absolute$sample), ];


x$pheno$met_type <- factor(x$pheno$met_type, levels=c("Primary", "BM"));
x$pheno$patient_id <- factor(x$pheno$patient_id);

x$pheno$purity <- absolute[["Cancer DNA fraction"]];


x$expr <- x$expr[!is.na(rownames(x$expr)), ];

# annotate amplifications
amp <- x$cn > 6.5;

amp.tfs <- c("YAP1", "MYC", "RUNX1");

mmp <- rownames(x$cn)[grep("MMP", rownames(x$cn))];
mmp.amp <- apply(amp[mmp, ], 2, any);

x$amp <- as.data.frame(t(rbind(amp[amp.tfs, ], MMP=mmp.amp)));

#any.tf.amp <- apply(x$amp[, amp.tfs], 1, sum) > 0
#x$amp.excl <- as.data.frame(lapply(
#	x$amp.all,
#	function(z) {
#		z[!z & any.tf.amp] <- NA;
#		z
#	}
#));
#rownames(x$amp.excl) <- rownames(x$amp);


# rename samples in exp
colnames(x$expr) <- gsub("-NT-SM-.+", "", colnames(x$expr))

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

names(emt.gsets) %<>%
	gsub("EPITHELIAL_(TO_)?MESENCHYMAL_TRANSITION", "EMT", .) %>%
	gsub("_?SIGNATURE", "", .) %>%
	gsub("_?RECEPTOR_SIGNALING", "", .) %>%
	gsub("MESENCHYMAL_TRANSITION", "EMT_UP", .) %>%
	gsub("EMT_EMT", "EMT", .) %>%
	gsub("_CANCER", "", .) %>%
	gsub("TGF_BETA", "TGFB", .);


mod.myc <- model.matrix(~ MYC, data=x$amp);
fit.myc <- eBayes(lmFit(x$expr, mod.myc));
tt.myc <- topTable(fit.myc, number=Inf);

tt.myc["MYC", ]
head(tt.myc);

yamanaka.factors <- c("SOX2", "POU5F1", "MYC", "KLF4");

hippo <- c("YAP1", "AMOTL2", "AMOT", "LATS2", "YWHAB", "DVL2", "WWC1", "AMOTL1", "LATS1");

tt.myc[yamanaka.factors, ]


expr.myc <- x$expr[rownames(x$expr) %in% hallmark$HALLMARK_MYC_TARGETS_V1, ];
fit.myc.sp <- eBayes(lmFit(expr.myc, mod.myc));
topTable(fit.myc.sp)

#sig.myc.sp <- decideTests(fit.myc.sp, p.value=0.1);
#expr.myc <- x$expr[sig.myc[,2] > 0, ];
#expr.myc

qwrite(tt.myc, "myc-amp.rnk", value.var="t");


mod.yap1 <- model.matrix(~ YAP1, data=x$amp);
fit.yap1 <- eBayes(lmFit(x$expr, mod.yap1));
tt.yap1 <- topTable(fit.yap1, number=Inf);

tt.yap1["YAP1", ]
head(tt.yap1);

sig.yap1 <- decideTests(fit.yap1, p.value=0.1);
expr.sig.yap1 <- x$expr[sig.yap1, ];

mod.yap1.p <- model.matrix(~ YAP1 + met_type:purity, data=cbind(x$amp, x$pheno));
fit.yap1.p <- eBayes(lmFit(x$expr, mod.yap1.p));
tt.yap1.p <- topTable(fit.yap1.p, coef="YAP1TRUE", number=Inf);

tt.yap1.p["YAP1", ]
tt.yap1.p[1:100, ]


#Heatmap(expr.sig.yap1, bottom_annotation=HeatmapAnnotation(df=x$amp), cluster_rows=FALSE, cluster_columns=FALSE);

tt.yap1[yamanaka.factors, ]

# 2008 Zhao, YAP1-inducible genes in MCF10A (mammary epithelial)
yap.zhao.up.targets <- c("CTGF", "ITGB2");
tt.yap1[toupper(yap.zhao.up.targets), ]
tt.yap1.p[toupper(yap.zhao.up.targets), ]

# 2014 Kapoor, upregulated in YAP1-bypassed mouse tumours (pancreatic ductal
# adenocarcinoma, iKras-) 
yap.kapoor.up.targets <- c("Ube2c", "Ccnb1", "Ccnb2", "Cenpi", "Cdc20", "Mcm3", "Mcm5", "Bub1", "Cdc6", "Rpa1", "Aurkb", "Mcm10", "Cdk1", "Cdc7", "Cdc25c", "Rad17", "Mcm6", "Pola1", "Ccnf", "Rad51", "Rfc3", "Rpa2", "Aurka", "Mcm2", "Rad1", "Rfc5", "Cdc25b", "Atr", "Chek1", "Smc2", "Bub3", "Cdkn1b", "Cenpl", "Uhrf1", "PsmB5", "Rfc2", "Mdm2", "Cdkn1a", "Birc5")
tt.yap1[toupper(yap.kapoor.up.targets), ]
tt.yap1.p[toupper(yap.kapoor.up.targets), ]

# 2014 Kapoor, overexpressed in YAP1 stably transfected and validated by ChIP
yap.kapoor.up.chip.targets <- c("Aurka", "Aurkb", "Ccnd1", "Ccnf", "Cdk1", "Mcm6", "PolA1", "Rad51", "Uhrf1");
tt.yap1[toupper(yap.kapoor.up.chip.targets), ]
tt.yap1.p[toupper(yap.kapoor.up.chip.targets), ]

# 2012 von Gise, Yap1-induced in cardiomyocytes
yap.vongise.up.targets <- c("AURKA", "AURKB", "CCNA2", "CCNB1", "CDK1", "CDC20", "CDC25B");
tt.yap1[toupper(yap.vongise.up.targets), ]
tt.yap1.p[toupper(yap.vongise.up.targets), ]

# 2014 Lorenzetto, downregulated by siYAP1 in YAP1-amplified cell lines (3)
# Cell lines: 
# CaSki, Cervical squamous cell carcinoma; metastasized to small intestine
# RO82, Follicular thyroid carcinoma
# EKVX, lung adenocarcinoma
yap.lorenzetto.up.targets <- c("CTGF", "SKP2", "GADD45A", "CCNA2");
tt.yap1[toupper(yap.lorenzetto.up.targets), ]
tt.yap1.p[toupper(yap.lorenzetto.up.targets), ]

# 2014 Lorenzetto, upregulated by siYAP1 in YAP1-amplified cell lines (3)
yap.lorenzetto.down.targets <- c("TEAD2", "VEGFA");
tt.yap1[toupper(yap.lorenzetto.down.targets), ]
tt.yap1.p[toupper(yap.lorenzetto.down.targets), ]

# 2015 Stein, SF268 (anaplastic astrocytoma)
# RNAseq and ChIP
yap.stein.up.targets <- c("CTGF", "ANKRD1", "CYR61", "NPPB", "CCND1", "AXL", "DKK1", "ITGB2", "WWC1", "KISS1", "NEXN", "PAWR", "S1PR1", "SNAPC1", "SKP2");
tt.yap1[toupper(yap.stein.up.targets), ]
tt.yap1.p[toupper(yap.stein.up.targets), ]

# Shunsuke Kitajima, David Barbie work in progress on A546
yap.kitajima.up.targets <- c("SETAD4?", "ANKRD1", "CTGF", "CDKN2C", "TNS1");
tt.yap1[toupper(yap.kitajima.up.targets), ]


tt.yap1[yamanaka.factors, ]
tt.yap1[hippo, ];


expr.no.mycamp <- x$expr[, !x$amp$MYC];
amp.no.mycamp <- x$amp[!x$amp$MYC, ];

mod.yap.no.mycamp <- model.matrix(~ YAP1, data=amp.no.mycamp);
fit.yap.no.mycamp <- eBayes(lmFit(expr.no.mycamp, mod.yap.no.mycamp));
tt.yap1.no.mycamp <- topTable(fit.yap.no.mycamp, number=Inf);

tt.yap1.no.mycamp[toupper(yap.zhao.up.targets), ]
tt.yap1.no.mycamp[toupper(yap.kapoor.up.targets), ]
tt.yap1.no.mycamp[toupper(yap.kapoor.up.chip.targets), ]
tt.yap1.no.mycamp[toupper(yap.vongise.up.targets), ]
tt.yap1.no.mycamp[toupper(yap.lorenzetto.up.targets), ]
tt.yap1.no.mycamp[toupper(yap.lorenzetto.down.targets), ]
tt.yap1.no.mycamp[toupper(yap.stein.up.targets), ]
tt.yap1.no.mycamp[yamanaka.factors, ];
tt.yap1.no.mycamp[hippo, ];


tt.yap1.no.mycamp["PPARGC1A", ]
tt.yap1.no.mycamp["SALL1", ]
tt.yap1.no.mycamp["FKBP5", ]


qwrite(tt.yap1, "yap1-amp.rnk", value.var="t");
qwrite(tt.yap1.no.mycamp, "yap1-amp_no-myc-amp.rnk", value.var="t");


mod.runx1 <- model.matrix(~ RUNX1, data=x$amp);
fit.runx1 <- eBayes(lmFit(x$expr, mod.runx1));
tt.runx1 <- topTable(fit.runx1, number=Inf);

tt.runx1["RUNX1", ]
head(tt.runx1);

tt.runx1[yamanaka.factors, ]


mod.runx1.no.mycamp <- model.matrix(~ RUNX1, data=amp.no.mycamp);
fit.runx1.no.mycamp <- eBayes(lmFit(expr.no.mycamp, mod.runx1.no.mycamp));
tt.runx1.no.mycamp <- topTable(fit.runx1.no.mycamp, number=Inf);

tt.runx1.no.mycamp["RUNX1", ]
tt.runx1.no.mycamp[yamanaka.factors, ]
tt.runx1.no.mycamp[hippo, ]


tt.runx1.no.mycamp[tt.runx1.no.mycamp$adj.P.Val < 0.1, ]

tt.runx1.no.mycamp[hippo, ]

qwrite(tt.runx1.no.mycamp, "runx1-amp_no-myc-amp.rnk", value.var="t");

