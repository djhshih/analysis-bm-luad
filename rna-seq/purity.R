library(io);
library(RColorBrewer);
library(magrittr);
library(limma);
library(dplyr);
library(ggplot2);
library(MASS);
library(reshape2);
library(tntrna);

x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");

pdf.fname <- filename("luad-bm_tntrna", ext="pdf");

colnames(x$expr) <- x$pheno$tumor_id;

markers <- c("GFAP", "NES", "CNP", "CD3G", "PTPRC", "ITGAM");
normal.markers <- x$expr[markers, ];
normal.markers.df <- data.frame(t(normal.markers));

purity <- x$pheno$purity;
quiescence.score <- x$pheno$quiescence;


fits <- apply(normal.markers, 1, function(y) lm(y ~ met_type + met_type:purity, data=x$pheno))
lapply(fits, summary)

fits <- apply(normal.markers, 1, function(y) lm(purity ~ met_type * y, data=x$pheno))
lapply(fits, summary)

wilcox.test(purity ~ met_type, data=x$pheno);
t.test(purity ~ met_type, data=x$pheno);
plot(purity ~ met_type, data=x$pheno);


d <- cbind(x$pheno, normal.markers.df);
fit0 <- lm(purity ~ 1, data = d);
fit1 <- lm(purity ~ GFAP + CNP + CD3G + PTPRC + ITGAM, data=d);
summary(fit1);
anova(fit0, fit1, test="Chisq");


met.idx <- x$pheno$met_type == "BM";

plot(x$expr["GFAP", met.idx], purity[met.idx])
cor(x$expr["GFAP", met.idx], purity[met.idx])^2
cor(x$expr["GFAP", met.idx], purity[met.idx], method="kendall")

# very low expression in primary: positive correlation is probably spurious
plot(x$expr["GFAP", !met.idx], purity[!met.idx])
cor(x$expr["GFAP", !met.idx], purity[!met.idx], method="kendall")

plot(x$expr["NES", !met.idx], purity[!met.idx])
cor(x$expr["NES", !met.idx], purity[!met.idx], method="kendall")

plot(x$expr["CNP", met.idx], purity[met.idx])
cor(x$expr["CNP", met.idx], purity[met.idx], method="kendall")

plot(x$expr["CNP", !met.idx], purity[!met.idx])
cor(x$expr["CNP", !met.idx], purity[!met.idx], method="kendall")

plot(x$expr["CD3G", met.idx], purity[met.idx])
cor(x$expr["CD3G", met.idx], purity[met.idx], method="kendall")

plot(x$expr["PTPRC", met.idx], purity[met.idx])
cor(x$expr["PTPRC", met.idx], purity[met.idx])^2
cor(x$expr["PTPRC", met.idx], purity[met.idx], method="kendall")

plot(x$expr["PTPRC", !met.idx], purity[!met.idx])
cor(x$expr["PTPRC", !met.idx], purity[!met.idx], method="kendall")

plot(x$expr["ITGAM", met.idx], purity[met.idx])
cor(x$expr["ITGAM", met.idx], purity[met.idx])^2
cor(x$expr["ITGAM", met.idx], purity[met.idx], method="kendall")


plot(x$expr["NKX2-1", met.idx], purity[met.idx])
cor(x$expr["NKX2-1", met.idx], purity[met.idx], method="kendall")

fit <- lm(x$expr["NKX2-1", !met.idx] ~ purity[!met.idx]);
plot(purity[!met.idx], x$expr["NKX2-1", !met.idx], xlim=c(0, 1), ylim=c(0, 7));
abline(a=coef(fit)[1], b=coef(fit)[2], col="red");

cor(x$expr["NKX2-1", !met.idx], purity[!met.idx])^2
cor(x$expr["NKX2-1", !met.idx], purity[!met.idx], method="kendall")
# NKX2 is a marker of primary lung adenocarcinoma and not metastatic lung
# adenocarcinoma

apply(normal.markers, 1, function(y) cor(y, purity));

# met
apply(normal.markers, 1, function(y) cor(y[met.idx], purity[met.idx], method="kendall"));

plot(x$expr["PTPRC", ], x$expr["ITGAM", ])

fit0 <- lm(purity[met.idx] ~ 1, data = cbind(x$pheno, normal.markers.df)[met.idx, ]);
fit1 <- lm(purity[met.idx] ~ GFAP + CNP + PTPRC, data=cbind(x$pheno, normal.markers.df)[met.idx, ]);
summary(fit1)
anova(fit0, fit1, test="Chisq");
0.33951 / 0.62489

# primary
apply(normal.markers, 1, function(y) cor(y[!met.idx], purity[!met.idx], method="kendall"));

plot(x$expr["GFAP", !met.idx], purity[!met.idx])

plot(x$expr["NES", !met.idx], purity[!met.idx])
# Eh?


# highly specific to met samples
plot(GFAP ~ met_type, data = cbind(x$pheno, normal.markers.df));

plot(NES ~ met_type, data = cbind(x$pheno, normal.markers.df));
plot(CNP ~ met_type, data = cbind(x$pheno, normal.markers.df));

# higher in primary
plot(CD3G ~ met_type, data = cbind(x$pheno, normal.markers.df));
plot(PTPRC ~ met_type, data = cbind(x$pheno, normal.markers.df));

# higher in met?
plot(ITGAM ~ met_type, data = cbind(x$pheno, normal.markers.df));

plot(NES ~ met_type, data = cbind(x$pheno, normal.markers.df));
plot(MAP2 ~ met_type, data = cbind(x$pheno, MAP2=x$expr["MAP2", ]));


plot(x$expr["ITGB2", ], quiescence.score);

normal.markers.mean <- apply(normal.markers, 2, mean);
plot(normal.markers.mean, purity);
cor(normal.markers.mean, purity, method="kendall");


# lung adenocarcinoma marker
# NKX2-1 (TTF-1)
# NKX2-1 is lineage specific marker for lung adenocarcinoma
# http://www.ncbi.nlm.nih.gov/pubmed/12023581
# http://www.ncbi.nlm.nih.gov/pubmed/11156325
# NKX2-1 (TTF-1) distinguish between primary and met lung adenocarcinoma
# http://link.springer.com/article/10.1007%2FBF02893461
plot(NKX2.1 ~ met_type, data = cbind(x$pheno, NKX2.1=x$expr["NKX2-1", ]));
t.test(NKX2.1 ~ met_type, data = cbind(x$pheno, NKX2.1=x$expr["NKX2-1", ]));
wilcox.test(NKX2.1 ~ met_type, data = cbind(x$pheno, NKX2.1=x$expr["NKX2-1", ]));


mod.mp <- model.matrix(~ purity[met.idx]);
fit.mp <- lmFit(x$expr[, met.idx], mod.mp);
eb.mp <- eBayes(fit.mp);
tt.mp <- topTable(eb.mp, number=Inf);
tt.mp[tt.mp$logFC < 0, ][1:100, ]
tt.mp[tt.mp$logFC > 0, ][1:100, ]


tt.mp["AGO2", ]
# particularly high slope
tt.mp["H19", ]

tt.mpi <- topTable(eb.mp, coef="(Intercept)", number=Inf);
nonstromal.genes <- rownames(tt.mpi)[abs(tt.mpi$logFC) < 0.1];

tt.mp[tt.mp$logFC > 0 & rownames(tt.mp) %in% nonstromal.genes, ]

mod.mpz <- model.matrix(~ purity[met.idx] - 1);
fit.mpz <- lmFit(x$expr[, met.idx], mod.mpz);
eb.mpz <- eBayes(fit.mpz);
tt.mpz <- topTable(eb.mpz, number=Inf);
tt.mpz[tt.mpz$logFC < 0, ][1:100, ]
tt.mpz[tt.mpz$logFC > 0, ][1:100, ]

tt.mpz["AGO2", ]


mod.emp <- model.matrix(~ purity[met.idx]);
fit.emp <- lmFit(expm1(x$expr[, met.idx]), mod.emp);
eb.emp <- eBayes(fit.emp);
tt.emp <- topTable(eb.emp, number=Inf);
tt.emp[tt.emp$logFC < 0, ][1:100, ]
tt.emp[tt.emp$logFC > 0, ][1:100, ]

tt.emp["AGO2", ]

logit <- function(x) {
	log( x / (1 - x) )
}
mod.lmp <- model.matrix(~ logit(purity[met.idx]));
fit.lmp <- lmFit(x$expr[, met.idx], mod.lmp);
eb.lmp <- eBayes(fit.lmp);
tt.lmp <- topTable(eb.lmp, number=Inf);
tt.lmp[tt.lmp$logFC < 0, ][1:100, ]
tt.lmp[tt.lmp$logFC > 0, ][1:100, ]

tt.lmp["AGO2", ]





# high slope of 11.6
plot(purity[met.idx], x$expr["H19", met.idx]);
# roboust method would get rid of spurious high effect size

fit <- lm(SEPHS1 ~ purity, data=data.frame(purity=purity[met.idx], SEPHS1=x$expr["SEPHS1", met.idx]));
plot(purity[met.idx], x$expr["SEPHS1", met.idx], xlim=c(0, 1))
abline(a=coef(fit)[1], b=coef(fit)[2], col="red")
tt.mp["SEPHS1", ]

tt.mp["GFAP", ]
tt.mp["CNP", ]
tt.mp["ITGAM", ]
tt.mp["PTPRC", ]
plot(x$expr["ITGAM", met.idx], purity[met.idx])

# FCGR2A and PLAUR and many other top genes negatively associated with purity
# are immune specific genes

plot(x$expr["PTPRC", met.idx], purity[met.idx])
cor(x$expr["PTPRC", met.idx], purity[met.idx])^2
cor(x$expr["PTPRC", met.idx], purity[met.idx], method="kendall")

plot(x$expr["FCGR2A", met.idx], purity[met.idx])
cor(x$expr["FCGR2A", met.idx], purity[met.idx])^2
cor(x$expr["FCGR2A", met.idx], purity[met.idx], method="kendall")

plot(x$expr["PLAUR", met.idx], purity[met.idx])
cor(x$expr["PLAUR", met.idx], purity[met.idx])^2
cor(x$expr["PLAUR", met.idx], purity[met.idx], method="kendall")


plot(x$expr["SEPHS1", met.idx], purity[met.idx])
cor(x$expr["SEPHS1", met.idx], purity[met.idx])^2

g <- ggplot(cbind(x$pheno, SEPHS1=x$expr["SEPHS1", ]), aes(y=SEPHS1, x=met_type, alpha=purity)) + geom_jitter(width=0.1)
g



plot(x$expr["AGO2", met.idx], purity[met.idx])
cor(x$expr["AGO2", met.idx], purity[met.idx])^2

plot(AGO2 ~ met_type, data = cbind(x$pheno, AGO2=x$expr["AGO2", ]));
t.test(AGO2 ~ met_type, data = cbind(x$pheno, AGO2=x$expr["AGO2", ]));
wilcox.test(AGO2 ~ met_type, data = cbind(x$pheno, AGO2=x$expr["AGO2", ]));

g <- ggplot(cbind(x$pheno, AGO2=x$expr["AGO2", ]), aes(y=AGO2, x=met_type, alpha=purity)) +
	geom_jitter(width=0.1)
g

d <- cbind(x$pheno, AGO2=x$expr["AGO2", ]);
fit <- lm(AGO2 ~ met_type * purity, data = d);
summary(fit)
# paper shows that AGO2 overexpression in lung adenocarcinoma cell line reduces growth

# same thing
fit <- lm(AGO2 ~ met_type + met_type:purity, data = d);
summary(fit)


fit <- lm(AGO2 ~ met_type:purity, data = d);
summary(fit)


# but AGO2 expression does not appear to be positively associated with quiescence
fit <- lm(AGO2 ~ quiescence.score + met_type + met_type:purity, data = d);
summary(fit)

plot(x$expr["AGO2", ], quiescence.score);
plot(x$expr["AGO2", met.idx], quiescence.score[met.idx]);
cor(x$expr["AGO2", met.idx], quiescence.score[met.idx], method="kendall");
plot(x$expr["AGO2", !met.idx], quiescence.score[!met.idx]);
cor(x$expr["AGO2", !met.idx], quiescence.score[!met.idx], method="kendall");

# another gene highly correlated with purity in met
plot(x$expr["TOX3", met.idx], purity[met.idx])
cor(x$expr["TOX3", met.idx], purity[met.idx])^2


mod.mpp <- model.matrix(~ met_type + met_type:purity, data=x$pheno);
fit.mpp <- lmFit(x$expr, mod.mpp);
eb.mpp <- eBayes(fit.mpp);
tt.mpp <- topTable(eb.mpp, coef="met_typeBM", number=Inf);
tt.mpp[tt.mpp$logFC < 0, ][1:10, ]
tt.mpp[tt.mpp$logFC > 0, ][1:10, ]

tt.mpp["AGO2", ]
tt.mpp["CD163", ]

g <- ggplot(cbind(x$pheno, CD163=x$expr["CD163", ]), aes(y=CD163, x=met_type, alpha=purity)) +
	geom_jitter(width=0.1)
g

# it might make more sense to focus on genes that are preferentially
# associated with purity in the met compartment
tt.mppm <- topTable(eb.mpp, coef="met_typeBM:purity", number=Inf);
tt.mppm[tt.mppm$logFC < 0, ][1:10, ]
tt.mppm[tt.mppm$logFC > 0, ][1:10, ]

tt.mppm["AGO2", ]

# top gene is likely driven by two outliers
g <- ggplot(cbind(x$pheno, SEPHS1=x$expr["SEPHS1", ]), aes(y=SEPHS1, x=met_type, alpha=purity)) +
	geom_jitter(width=0.1)
g

# not sure what this does
g <- ggplot(cbind(x$pheno, STRBP=x$expr["STRBP", ]), aes(y=STRBP, x=met_type, alpha=purity)) +
	geom_jitter(width=0.1)
g

# low expression: could be spurious
g <- ggplot(cbind(x$pheno, KCNH2=x$expr["KCNH2", ]), aes(y=KCNH2, x=met_type, alpha=purity)) +
	geom_jitter(width=0.1)
g

# co-receptor of WNT signaling
g <- ggplot(cbind(x$pheno, LRP6=x$expr["LRP6", ]), aes(y=LRP6, x=met_type, alpha=purity)) +
	geom_jitter(width=0.1)
g



mod.pp <- model.matrix(~ purity[!met.idx]);
fit.pp <- lmFit(x$expr[, !met.idx], mod.pp);
eb.pp <- eBayes(fit.pp);
tt.pp <- topTable(eb.pp, number=Inf);
tt.pp[tt.pp$logFC < 0, ][1:10, ]
tt.pp[tt.pp$logFC > 0, ][1:10, ]

tt.pp["NKX2-1", ]
tt.pp["NAPSA", ]
tt.pp["AGO2", ]

# CLDN3 is specifically expressed in lung adenocarcinoma
# http://www.ncbi.nlm.nih.gov/pubmed/25935651
plot(x$expr["CLDN3", met.idx], purity[met.idx])
cor(x$expr["CLDN3", met.idx], purity[met.idx], method="kendall")
plot(x$expr["CLDN3", !met.idx], purity[!met.idx])
cor(x$expr["CLDN3", !met.idx], purity[!met.idx])^2
cor(x$expr["CLDN3", !met.idx], purity[!met.idx], method="kendall")

plot(x$expr["GJB3", met.idx], purity[met.idx])
plot(x$expr["GJB3", !met.idx], purity[!met.idx])


plot(purity, quiescence.score)
plot(purity[met.idx], quiescence.score[met.idx])
plot(purity[!met.idx], quiescence.score[!met.idx])

plot(x$expr["AGO2", met.idx], quiescence.score[met.idx]);


# Which genes are more preferentially degraded?
# Longer genes?
# Exome-capture: so expect no 3' to 5' bias; typically based by polyA
# selection or oligodT priming

fit.met <- rlm(AGO2 ~ purity, data=data.frame(purity=purity[met.idx], AGO2=x$expr["AGO2", met.idx]));
plot(purity[met.idx], x$expr["AGO2", met.idx], xlim=c(0, 1), ylim=c(0, 4.5))
abline(a=coef(fit.met)[1], b=coef(fit.met)[2], col="red")

pred.met <- predict(fit.met, data.frame(purity=1.0), se.fit=TRUE);


fit.pr <- rlm(AGO2 ~ purity, data=data.frame(purity=purity[!met.idx], AGO2=x$expr["AGO2", !met.idx]));
plot(purity[!met.idx], x$expr["AGO2", !met.idx], xlim=c(0, 1), ylim=c(0, 4.5))
abline(a=coef(fit.pr)[1], b=coef(fit.pr)[2], col="red")

pred.pr <- predict(fit.pr, data.frame(purity=1.0), se.fit=TRUE);

se <- sqrt(pred.met$se.fit^2 + pred.pr$se.fit^2);
z <- (pred.met$fit - pred.pr$fit) / se;

# right-tail test
pval <- 1 - pnorm(z)


purified_t_test(x$expr["AGO2", ], x$pheno$met_type, x$pheno$purity);
purified_t_test(x$expr["AGO2", ], x$pheno$met_type, x$pheno$purity, alternative="less");

met.genes <- rownames(tt.mp)[tt.mp$logFC > 0 & tt.mp$P.Value < 0.01];

hs.met <- apply(x$expr[met.genes, ], 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="less");

ps.met <- unlist(lapply(hs.met, function(x) x$p.value));
#names(ps) <- sub(".z", "", names(ps), fixed=TRUE);
d.met <- data.frame(
	logFC = unlist(lapply(hs.met, function(x) x$effect)),
	y2 = unlist(lapply(hs.met, function(x) x$estimate[2])),
	t = unlist(lapply(hs.met, function(x) x$statistic)),
	p = ps.met,
	q = p.adjust(ps.met, method="BH"),
	row.names = names(ps.met)
);
d.met <- d.met[order(d.met$p), ]

y1 <- unlist(lapply(hs.met, function(x) x$estimate[1]));
y2 <- unlist(lapply(hs.met, function(x) x$estimate[2]));
hist(y1);
hist(y2);


met.strom.genes <- rownames(tt.mp)[tt.mp$logFC < 0 & tt.mp$P.Value < 0.01];

hs.met.strom <- apply(x$expr[met.strom.genes, ], 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="less",
	target.purity=0);

ps.met.strom <- unlist(lapply(hs.met.strom, function(x) x$p.value));
#names(ps) <- sub(".z", "", names(ps), fixed=TRUE);
d.met.strom <- data.frame(
	logFC = unlist(lapply(hs.met.strom, function(x) x$effect)),
	y2 = unlist(lapply(hs.met.strom, function(x) x$estimate[2])),
	t = unlist(lapply(hs.met.strom, function(x) x$statistic)),
	p = ps.met.strom,
	q = p.adjust(ps.met.strom, method="BH"),
	row.names = names(ps.met.strom)
);
d.met.strom <- d.met.strom[order(d.met.strom$p), ]



hs.met.z <- apply(x$expr[met.genes, ], 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, zero.intercept=TRUE, alternative="less");

# plausible
purified_group_plot(x$expr["LRP6", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["LRP6", ], x$pheno$met_type, x$pheno$purity);

purified_group_plot(x$expr["SEPHS1", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["SEPHS1", ], x$pheno$met_type, x$pheno$purity);

# plausible, but
# 3 samples in primary causing signifiance...
purified_group_plot(x$expr["AGO2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["AGO2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(expm1(x$expr["AGO2", ]), x$pheno$met_type, x$pheno$purity);
#purity_plot(x$expr["AGO2", ], x$pheno$met_type, logit(x$pheno$purity))

# very high variability
purity_plot(x$expr["H19", ], x$pheno$met_type, x$pheno$purity);

# primary trend is down: not plausible
purity_plot(x$expr["TTC30B", ], x$pheno$met_type, x$pheno$purity);

purified_group_plot(x$expr["SH3BP4", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["SH3BP4", ], x$pheno$met_type, x$pheno$purity);
purity_plot(expm1(x$expr["SH3BP4", ]), x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["SH3BP4", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);

# low expression... spurious?
# negative y-intercept
purified_group_plot(x$expr["KCNH2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["KCNH2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["KCNH2", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);

# negative y-intercept
purified_group_plot(x$expr["DGCR6", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["DGCR6", ], x$pheno$met_type, x$pheno$purity);


# non-stromal genes that are correlated with purity in BM
purity_plot(x$expr["ACSS1", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);

purity_plot(x$expr["SIX1", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);

purity_plot(x$expr["PRKAR2B", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);

purity_plot(x$expr["STK33", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);


# genes suspected to be expressed in the stroma
purity_plot(x$expr["RUNX1", ], x$pheno$met_type, x$pheno$purity);


# astrocyte
g <- purity_plot(x$expr["GFAP", ], x$pheno$met_type, x$pheno$purity);
qdraw(g, filename("purity-plot", tag="GFAP", ext="pdf"));

purity_plot(x$expr["NES", ], x$pheno$met_type, x$pheno$purity);

purified_t_test(x$expr["GFAP", ], x$pheno$met_type, x$pheno$purity, alternative="less", target.purity=0.0);

# pericyte
purity_plot(x$expr["PDGFRB", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["ANPEP", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["DES", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["ACTA2", ], x$pheno$met_type, x$pheno$purity);

plot(factor(x$mut$XRCC5), x$expr["XRCC5", ]);

# mutated gene (DNA repair)
purity_plot(x$expr["XRCC1", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["XRCC2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["XRCC3", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["XRCC4", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["XRCC5", ], x$pheno$met_type, x$pheno$purity);


# amp gnees

d <- data.frame(
	y = x$expr["MYC", ],
	group = x$pheno$met_type, 
	purity = x$pheno$purity,
	amp = factor(x$amp[, "MYC"]),
	cn = log2(x$cn["MYC", ] + 1)
);
purity_plot(x$expr["MYC", ], x$pheno$met_type, x$pheno$purity, points=FALSE) +
	geom_point(aes(x = purity, y = y, colour = group, shape=amp, size=cn), d);

d <- data.frame(
	y = x$expr["YAP1", ],
	group = x$pheno$met_type, 
	purity = x$pheno$purity,
	amp = factor(x$amp[, "YAP1"]),
	cn = log2(x$cn["YAP1", ] + 1)
);
purity_plot(x$expr["YAP1", ], x$pheno$met_type, x$pheno$purity, points=FALSE) +
	geom_point(aes(x = purity, y = y, colour = group, shape=amp, size=cn), d);


# other genes
purity_plot(x$expr["HIF1A", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["HIF1AN", ], x$pheno$met_type, x$pheno$purity);


pri.genes <- rownames(tt.pp)[tt.pp$logFC > 0 & tt.pp$P.Value < 0.01];

hs.pri <- apply(x$expr[pri.genes, ], 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="greater");

ps.pri <- unlist(lapply(hs.pri, function(x) x$p.value));
#names(ps) <- sub(".z", "", names(ps), fixed=TRUE);
d.pri <- data.frame(
	logFC = unlist(lapply(hs.pri, function(x) x$effect)),
	y = unlist(lapply(hs.pri, function(x) x$estimate[1])),
	t = unlist(lapply(hs.pri, function(x) x$statistic)),
	p = ps.pri,
	q = p.adjust(ps.pri, method="BH"),
	row.names = names(ps.pri)
);
d.pri <- d.pri[order(d.pri$p), ]

hist(d.pri$y);


pri.strom.genes <- rownames(tt.pp)[tt.pp$logFC < 0 & tt.pp$P.Value < 0.01];

hs.pri.strom <- apply(x$expr[pri.strom.genes, ], 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="greater",
	target.purity=0);

ps.pri.strom <- unlist(lapply(hs.pri.strom, function(x) x$p.value));
#names(ps) <- sub(".z", "", names(ps), fixed=TRUE);
d.pri.strom <- data.frame(
	logFC = unlist(lapply(hs.pri.strom, function(x) x$effect)),
	y = unlist(lapply(hs.pri.strom, function(x) x$estimate[1])),
	t = unlist(lapply(hs.pri.strom, function(x) x$statistic)),
	p = ps.pri.strom,
	q = p.adjust(ps.pri.strom, method="BH"),
	row.names = names(ps.pri.strom)
);
d.pri.strom <- d.pri.strom[order(d.pri.strom$p), ]

hist(d.pri.strom$y);
hist(unlist(lapply(hs.pri.strom, function(x) x$estimate[2])));



d.pri["NKX2-1", ]

purified_t_test(x$expr["NKX2-1", ], x$pheno$met_type, x$pheno$purity, alternative="greater");
purified_t_test(x$expr["CLDN3", ], x$pheno$met_type, x$pheno$purity, alternative="greater");

purified_group_plot(x$expr["NKX2-1", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["NKX2-1", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["NKX2-1", ], x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);
purity_plot(expm1(x$expr["NKX2-1", ]), x$pheno$met_type, x$pheno$purity, zero.intercept=TRUE);

g <- purity_plot(x$expr["NKX2-1", ], x$pheno$met_type, x$pheno$purity, points=FALSE) +
	geom_point(aes(x = purity, y = y, colour = group, shape = x$pheno$histology_group)) +
	scale_shape_discrete(name="histotype") +
	geom_text(aes(x = purity, y = y, label = x$pheno$patient_id), nudge_y = 0.2, size=1.5) +
	ylab("log2(TPM + 1) expression") + xlab("purity");

qdraw(g, file=tag(pdf.fname, c("NKX2-1")), width=7);



luad.idx <- x$pheno$histology_group == "Lung adenocarcinoma";
x.luad <- list(
	expr = x$expr[, luad.idx],
	pheno = x$pheno[luad.idx, ],
	cn = x$cn[, luad.idx]
);

g <- purity_plot(x.luad$expr["NKX2-1", ], factor(rep("Lung adenocarcinoma", nrow(x.luad$pheno))),
	x.luad$pheno$purity, points=FALSE) +
	guides(fill=FALSE) +
	geom_point(aes(x = x$pheno$purity, y = x$expr["NKX2-1", ], colour = x$pheno$met_type,
		shape = x$pheno$histology_group)) +
	scale_shape_discrete(name="histotype") +
	geom_text(aes(x = x$pheno$purity,  y = x$expr["NKX2-1", ], label = x$pheno$patient_id),
		nudge_y = 0.2, size=1.5) +
	ylab("log2(TPM + 1) expression") + xlab("purity") + ggtitle("NKX2-1 (TTF-1)");
qdraw(g, file=tag(pdf.fname, c("NKX2-1", "fit-luad")), width=7);

g <- purity_plot(x.luad$expr["NKX2-1", ], factor(rep("Lung adenocarcinoma", nrow(x.luad$pheno))),
	x.luad$pheno$purity, points=FALSE) +
	guides(fill=FALSE) +
	geom_point(aes(x = x$pheno$purity, y = x$expr["NKX2-1", ], colour = x$pheno$met_type,
		shape = x$pheno$histology_group)) +
	scale_shape_discrete(name="histotype") +
	geom_text(aes(x = x$pheno$purity,  y = x$expr["NKX2-1", ], label = x$pheno$patient_id),
		nudge_y = 0.2, size=1.5) +
	ylab("log2(TPM + 1) expression") + xlab("purity") + ggtitle("NKX2-1 (TTF-1)");
qdraw(g, file=tag(pdf.fname, c("NKX2-1", "fit-luad")), width=7);


g <- purity_plot(x.luad$expr["CLDN3", ], factor(rep("Lung adenocarcinoma", nrow(x.luad$pheno))),
	x.luad$pheno$purity, points=FALSE) +
	guides(fill=FALSE) +
	geom_point(aes(x = x$pheno$purity, y = x$expr["CLDN3", ], colour = x$pheno$met_type,
		shape = x$pheno$histology_group)) +
	scale_shape_discrete(name="histotype") +
	geom_text(aes(x = x$pheno$purity,  y = x$expr["CLDN3", ], label = x$pheno$patient_id),
		nudge_y = 0.2, size=1.5) +
	ylab("log2(TPM + 1) expression") + xlab("purity") + ggtitle("CLDN3");
qdraw(g, file=tag(pdf.fname, c("CLDN3", "fit-luad")), width=7);

g <- purity_plot(x.luad$expr["NAPSA", ], factor(rep("Lung adenocarcinoma", nrow(x.luad$pheno))),
	x.luad$pheno$purity, points=FALSE) +
	guides(fill=FALSE) +
	geom_point(aes(x = x$pheno$purity, y = x$expr["NAPSA", ], size = x$cn["NAPSA", ], colour = x$pheno$met_type,
		shape = x$pheno$histology_group)) +
	scale_shape_discrete(name="histotype") +
	geom_text(aes(x = x$pheno$purity,  y = x$expr["NAPSA", ], label = x$pheno$patient_id),
		nudge_y = 0.2, size=1.5) +
	ylab("log2(TPM + 1) expression") + xlab("purity") + ggtitle("NAPSA");
qdraw(g, file=tag(pdf.fname, c("NAPSA", "fit-luad")), width=7);


napsa <-  x$expr["NAPSA", ] / x$cn["NAPSA", ];
g <- purity_plot(x.luad$expr["NAPSA", ] / x.luad$cn["NAPSA", ], factor(rep("Lung adenocarcinoma", nrow(x.luad$pheno))),
	x.luad$pheno$purity, points=FALSE) +
	guides(fill=FALSE) +
	geom_point(aes(x = x$pheno$purity, y = napsa, colour = x$pheno$met_type,
		shape = x$pheno$histology_group)) +
	scale_shape_discrete(name="histotype") +
	geom_text(aes(x = x$pheno$purity,  y = napsa, label = x$pheno$patient_id),
		nudge_y = 0.2, size=1.5) +
	ylab("log2(TPM + 1) expression") + xlab("purity") + ggtitle("NAPSA");
qdraw(g, file=tag(pdf.fname, c("NAPSA", "fit-luad")), width=7);



d <- data.frame(
	y = x$expr["NKX2-1", ],
	group = x$pheno$met_type, 
	purity = x$pheno$purity,
	cn = log2(x$cn["NKX2-1", ] + 1)
);
purity_plot(x$expr["NKX2-1", ], x$pheno$met_type, x$pheno$purity, points=FALSE) +
	geom_point(aes(x = purity, y = y, colour = group, size=cn), d);


purified_group_plot(x$expr["CLDN3", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["CLDN3", ], x$pheno$met_type, x$pheno$purity);


purity_plot(x$expr["FAM127B", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["ZNF84", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["HOOK3", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["RHOU", ], x$pheno$met_type, x$pheno$purity);

purified_group_plot(x$expr["DYNC1LI2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["DYNC1LI2", ], x$pheno$met_type, x$pheno$purity);

purified_group_plot(x$expr["GRIP1", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["GRIP1", ], x$pheno$met_type, x$pheno$purity);

purified_group_plot(x$expr["PHC2", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["PHC2", ], x$pheno$met_type, x$pheno$purity);

purity_plot(x$expr["ICA1", ], x$pheno$met_type, x$pheno$purity);
purified_t_test(x$expr["ICA1", ], x$pheno$met_type, x$pheno$purity, alternative="greater");

purity_plot(x$expr["SNTB1", ], x$pheno$met_type, x$pheno$purity);


# negatively correlated with primary purity
# value should not become negative at purity of 1.0
# GJB2 and GJB3 are gap junction proteins: expressed in both normal lung and
# cerebral tissue at low level
# tumour cells downregulate this gene presumably during EMT
purity_plot(x$expr["GJB3", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["GJB2", ], x$pheno$met_type, x$pheno$purity);
# low expression
purity_plot(x$expr["TNIP3", ], x$pheno$met_type, x$pheno$purity);
# prion protein...
purity_plot(x$expr["PRNP", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["CD44", ], x$pheno$met_type, x$pheno$purity);

# negatively correlated with met purity
purity_plot(x$expr["FCGR2A", ], x$pheno$met_type, x$pheno$purity);
purity_plot(expm1(x$expr["FCGR2A", ]), x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["FCER1G", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["FCGR3A", ], x$pheno$met_type, x$pheno$purity);
# plasminogen activator
purity_plot(x$expr["PLAUR", ], x$pheno$met_type, x$pheno$purity);
# membrane protein 1
purity_plot(x$expr["MPP1", ], x$pheno$met_type, x$pheno$purity);

purity_plot(x$expr["OLR1", ], x$pheno$met_type, x$pheno$purity);
# complement proteins
purity_plot(x$expr["C1QB", ], x$pheno$met_type, x$pheno$purity);
purity_plot(x$expr["C1QC", ], x$pheno$met_type, x$pheno$purity);

# macrophage marker
purity_plot(x$expr["CD68", ], x$pheno$met_type, x$pheno$purity);

# cathepsin (probably expressed by neutrophils)
purity_plot(x$expr["CTSB", ], x$pheno$met_type, x$pheno$purity);


# Plot coefficients

coefs.mp <- data.frame(coef(eb.mp));
colnames(coefs.mp) <- c("beta0", "beta1");
coefs.mp$p <- tt.mp$P.Value[match(rownames(coefs.mp), rownames(tt.mp))];

coefs.pp <- data.frame(coef(eb.pp));
colnames(coefs.pp) <- c("beta0", "beta1");
coefs.pp$p <- tt.pp$P.Value[match(rownames(coefs.pp), rownames(tt.pp))];

coefs.pp.mp <- rbind(
	data.frame(gene = rownames(coefs.pp), coefs.pp, met_type="Primary"),
	data.frame(gene = rownames(coefs.mp), coefs.mp, met_type="BM")
);
coefs.pp.mp$met_type <- factor(coefs.pp.mp$met_type, levels=levels(x$pheno$met_type));

g <- ggplot(coefs.pp.mp, aes(x=beta0, y=beta0+beta1*1, alpha=-log10(p), colour=met_type)) +
	geom_point() + theme_minimal() + 
	geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey60") + 
	xlab("non-tumour expression") + ylab("tumour expression");
g

g <- ggplot(filter(coefs.pp.mp, p < 0.05), aes(x=beta0, y=beta0+beta1*1, alpha=-log10(p), colour=met_type)) +
	geom_point() + theme_minimal() + 
	geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey60") + 
	xlab("non-tumour expression") + ylab("tumour expression");
g

ggplot(data.frame(x=as.numeric(x$expr)), aes(x=x)) + geom_histogram(bins=100);
ggplot(coefs.pp.mp, aes(x=beta01)) + geom_histogram(bins=100);

d <- with(coefs.pp.mp, data.frame(
	nontumour = beta0,
	tumour = beta0 + beta1
));
d <- melt(d);
ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=0.5, linetype=0);


ggplot(x$pheno, aes(x=met_type, y=purity)) + geom_jitter(width=0.1) + theme_minimal();


# take the negative of t since t is done by Primary - BM
coefs.mp.t <- cbind(
	coefs.mp[match(rownames(d.met), rownames(coefs.mp)), ],
	t = - d.met$t, y2=d.met$y2,
	gene = ifelse(d.met$q < 0.25, rownames(d.met), "")
);

g <- ggplot(coefs.mp.t, aes(x=beta0, y=beta0+beta1*1, colour=t, label=gene)) +
	geom_point() + 
	scale_colour_gradient2(low="green4", mid="white", high="firebrick") +
	geom_text(hjust=0, nudge_x = 0.1, size=3, colour="grey30") +
	theme_bw() + 
	geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey60") + 
	xlab("non-tumour expression") + ylab("tumour expression");
g

coefs.mp.s.t <- cbind(
	coefs.mp[match(rownames(d.met.strom), rownames(coefs.mp)), ],
	t = - d.met.strom$t, y2=d.met.strom$y2,
	gene = ifelse(d.met.strom$q < 0.1, rownames(d.met.strom), "")
);

g <- ggplot(coefs.mp.s.t, aes(x=beta0, y=beta0+beta1*1, colour=t, label=gene)) +
	geom_point() + 
	scale_colour_gradient2(low="green4", mid="white", high="firebrick") +
	geom_text(hjust=0, nudge_x = 0.1, size=3, colour="grey30") +
	theme_bw() + 
	geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey60") + 
	xlab("non-tumour expression") + ylab("tumour expression");
g


coefs.pp.t <- cbind(
	coefs.pp[match(rownames(d.pri), rownames(coefs.pp)), ],
	t=d.pri$t, y=d.pri$y,
	gene = ifelse(d.pri$q < 0.25, rownames(d.pri), "")
);

g <- ggplot(coefs.pp.t, aes(x=beta0, y=beta0+beta1*1, colour=t, label=gene)) +
	geom_point() + 
	scale_colour_gradient2(low="green", mid="white", high="red") +
	geom_text(hjust=0, nudge_x = 0.1, size=3, colour="grey30") +
	theme_minimal() + 
	geom_abline(intercept=0, slope=1) + xlab("non-tumour expression") + ylab("tumour expression");
g



coefs.pp.s.t <- cbind(
	coefs.pp[match(rownames(d.pri.strom), rownames(coefs.pp)), ],
	t=d.pri.strom$t, y=d.pri.strom$y,
	gene = ifelse(d.pri.strom$q < 0.25, rownames(d.pri.strom), "")
);

g <- ggplot(coefs.pp.s.t, aes(x=beta0, y=beta0+beta1*1, colour=t, label=gene)) +
	geom_point() + 
	scale_colour_gradient2(low="green", mid="white", high="red") +
	geom_text(hjust=0, nudge_x = 0.1, size=3, colour="grey30") +
	theme_minimal() + 
	geom_abline(intercept=0, slope=1) + xlab("non-tumour expression") + ylab("tumour expression");
g


# ST6GALNAC5 is upregulated in brain metastasis from breast cancer
# Bos, Nature 2009
st6 <- grep("^ST6GALNAC", rownames(x$expr), value=TRUE);
tt.pp[st6, ]
tt.mp[st6, ]

hs.met.st6 <- apply(x$expr[st6, ], 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="less");


####

hs.all <- apply(x$expr, 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="less");

ps.all <- unlist(lapply(hs.all, function(x) x$p.value));
d.all <- data.frame(
	logFC = unlist(lapply(hs.all, function(x) x$effect)),
	y1 = unlist(lapply(hs.all, function(x) x$estimate[1])),
	y2 = unlist(lapply(hs.all, function(x) x$estimate[2])),
	t = unlist(lapply(hs.all, function(x) x$statistic)),
	p = ps.all,
	q = p.adjust(ps.all, method="BH"),
	row.names = names(ps.all)
);

hist(d.all$y1, breaks=1000);
hist(d.all$y2, breaks=1000);
smoothScatter(d.all$y1, d.all$y2);


hs.all.strom <- apply(x$expr, 1, purified_t_test,
	group = x$pheno$met_type, purity = x$pheno$purity, alternative="less", target.purity=0);

ps.all.strom <- unlist(lapply(hs.all.strom, function(x) x$p.value));
d.all.strom <- data.frame(
	logFC = unlist(lapply(hs.all.strom, function(x) x$effect)),
	y1 = unlist(lapply(hs.all.strom, function(x) x$estimate[1])),
	y2 = unlist(lapply(hs.all.strom, function(x) x$estimate[2])),
	t = unlist(lapply(hs.all.strom, function(x) x$statistic)),
	p = ps.all.strom,
	q = p.adjust(ps.all.strom, method="BH"),
	row.names = names(ps.all.strom)
);

hist(d.all.strom$y1, breaks=1000);
hist(d.all.strom$y2, breaks=1000);
smoothScatter(d.all.strom$y1, d.all.strom$y2);
cor(d.all.strom$y1, d.all.strom$y2);

beta11 <- d.all$y1 - d.all.strom$y1
beta12 <- d.all$y2 - d.all.strom$y2
hist(beta11, breaks=1000);
hist(beta12, breaks=1000);

smoothScatter(beta11, beta12);
abline(h=0)
abline(v=0)

