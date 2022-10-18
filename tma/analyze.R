library(io);
library(dplyr);
library(ggplot2);

prim <- qread("primary-tma.tsv");
bm <- qread("bm-tma.tsv");

d <- rbind(
	data.frame(prim, group = "Primary"),
	data.frame(bm, group = "Brain metastasis")
);

qwrite(d, "tma_primary-bm.tsv");

with(d, table(group, YAP1))
with(d, table(group, MYC))
with(d, table(group, PIK3CB))

with(d, fisher.test(table(group, YAP1)))
with(d, fisher.test(table(group, MYC)))

with(d, table(group, YAP1 == 1))
with(d, table(group, MYC == 1))
with(d, table(group, PIK3CB == 1))

fit <- glm(MYC == 1 ~ group, family=binomial(link = "logit"), data=d);
summary(fit)

fit <- glm(MYC == 1 ~ group, family=quasibinomial, data=d);
summary(fit)

with(d, fisher.test(table(group, YAP1 == 1), alternative="greater"))
with(d, fisher.test(table(group, MYC == 1), alternative="greater"))

contrast <- function(x, group0, group1) {
	x[! x %in% c(group0, group1)] <- NA;
	x == group1
}

with(d, table(group, contrast(YAP1, 0, 1)))
with(d, table(group, contrast(MYC, 0, 1)))
with(d, table(group, contrast(PIK3CB, 0, 1)))

with(d, fisher.test(table(group, contrast(YAP1, 0, 1))))
with(d, fisher.test(table(group, contrast(MYC, 0, 1))))
with(d, fisher.test(table(group, contrast(PIK3CB, 0, 1))))

with(d, fisher.test(table(group, contrast(YAP1, 0, 1)), alternative="greater"))
with(d, fisher.test(table(group, contrast(MYC, 0, 1)), alternative="greater"))

d <- mutate(d, yap1_amp = contrast(YAP1, 0, 1));
d <- mutate(d, myc_amp = contrast(MYC, 0, 1));

with(d, table(group, yap1_amp))
with(d, table(group, myc_amp))

fit <- glm(yap1_amp ~ group, family=binomial(link = "logit"), data=d);
summary(fit)

fit <- glm(myc_amp ~ group, family=binomial(link = "logit"), data=d);
summary(fit)

fit <- glm(myc_amp ~ group, family=quasibinomial, data=d);
summary(fit)


qdraw(
	ggplot(d, aes(x = group, fill = factor(YAP1))) +
		geom_bar(position="fill") + theme_bw() +
		ggtitle("YAP1") + ylab("Proportion") + xlab("")
	,
	file = "tma-bar_YAP1.pdf"
);

qdraw(
	ggplot(d, aes(x = group, fill = factor(MYC))) +
		geom_bar(position="fill") + theme_bw() +
		ggtitle("MYC") + ylab("Proportion") + xlab("")
	,
	file = "tma-bar_MYC.pdf"
);

