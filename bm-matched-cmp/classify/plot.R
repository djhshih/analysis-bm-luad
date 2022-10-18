library(io);
library(ggplot2);
library(reshape2)
library(dplyr)
library(binom)

count.amp.drivers <- qread("bm-comparts_amp-drivers.rds");
count.del.drivers <- qread("bm-comparts_del-drivers.rds");
count.amp <- qread("bm-comparts_amp-genes.rds");
count.del<- qread("bm-comparts_del-genes.rds");

out.fname <- filename("pkb-luad", tag="bm-matched-cmp");
pdf.fname <- insert(out.fname, ext="pdf");
tsv.fname <- insert(out.fname, ext="tsv");

compart.levels <- c("Primary", "Shared", "Brain metastasis");

# rearrange columns
count.amp <- count.amp[, compart.levels];
count.del <- count.del[, compart.levels];
count.amp.drivers <- count.amp.drivers[, compart.levels];
count.del.drivers <- count.del.drivers[, compart.levels];

alpha.level <- 0.20;

source("~/exigens/brain-mets/common/theme.R");

####

count_matrix_to_d <- function(x) {
	d <- melt(x, varnames = c("gene", "compartment"), value.name = "count");
	d$compartment <- factor(d$compartment, levels = compart.levels);
	d
}

count_matrix_to_d_marginal <- function(x) {
	z <- apply(x, 2, sum);
	d <- data.frame(
		gene = NA,
		compartment = names(z),
		count = z
	)
	rownames(d) <- NULL;
	d
}

tests_to_string <- function(hs, show.statistic=FALSE) {
	paste(
		unlist(lapply(hs, test_to_string, show.statistic=show.statistic)),
		collapse = "\n"
	);
}

test_to_string <- function(h, show.statistic=FALSE) {
	if (show.statistic) {
		s0 <- sprintf("%s = %.2f", names(h$statistic), h$statistic);
	} else {
		s0 <- NULL;
	}
	s1 <- p_to_string(h$p.value);
	if (is.null(s0)) {
		s <- s1;
	} else {
		s <- paste(s0, s1, sep=",  ");
	}
	s
}

p_to_string <- function(p, abbrev=TRUE) {
	if (abbrev) {
		ifelse(p < 0.001,
			"p < 0.001",
			ifelse(p < 0.01,
				sprintf("p = %.4f", p),
				ifelse(p < 0.05,
					ifelse(round(p, 2) - 0.05 == 0,
						"p < 0.05",
						sprintf("p = %.3f", p)
					),
					ifelse(p < 0.10,
						sprintf("p = %.3f", p),
						sprintf("p = %.2f", p)
					)
				)
			)
		)
	} else {
		sprintf("p = %f", p)
	}
}

####


compart_to_branch <- function(x) {
	factor(x, levels = compart.levels, labels = sub("Shared", "Trunk", compart.levels))
}


####

compart.amp.drivers <- count_matrix_to_d(count.amp.drivers);
compart.amp.genes <- count_matrix_to_d_marginal(count.amp);
compart.amp <- rbind(
	mutate(compart.amp.drivers, group = "Driver genes"),
	mutate(compart.amp.genes, group = "Other genes")
);

mean.amp <- apply(count.amp, 2, mean);

compart.amp2 <- rbind(
	compart.amp.drivers,
	data.frame(
		gene = "Expected",
		compartment = names(mean.amp),
		count = mean.amp
	)
);
rownames(compart.amp2) <- NULL;
compart.amp2$branch <- compart_to_branch(compart.amp2$compartment);


####

compart.del.drivers <- count_matrix_to_d(count.del.drivers);
compart.del.genes <- count_matrix_to_d_marginal(count.del);
compart.del <- rbind(
	mutate(compart.del.drivers, group = "Driver genes"),
	mutate(compart.del.genes, group = "Other genes")
);

mean.del <- apply(count.del, 2, mean);

compart.del2 <- rbind(
	compart.del.drivers,
	data.frame(
		gene = "Expected",
		compartment = names(mean.del),
		count = mean.del
	)
);
rownames(compart.del2) <- NULL;
compart.del2$branch <- compart_to_branch(compart.del2$compartment);

####

gmarginal.amp <- apply(count.amp, 1, sum);
gmarginal.amp.drivers <- apply(count.amp.drivers, 1, sum);


mean(gmarginal.amp)
fit <- glm(gmarginal.amp ~ 1, family="quasipoisson");
exp(coef(fit))

gmarginal.amp.d <- data.frame(
	group = c(rep("other", length(gmarginal.amp)), rep("driver", length(gmarginal.amp.drivers))),
	count = c(gmarginal.amp, gmarginal.amp.drivers)
);

fit <- glm(count ~ group, data=gmarginal.amp.d, family="quasipoisson");
summary(fit)

lapply(gmarginal.amp.drivers,
	function(x) {
		poisson.test(x, r = mean(gmarginal.amp), alternative="greater")
	}
);

h.amp.pois <- poisson.test(mean(gmarginal.amp.drivers), r = mean(gmarginal.amp), alternative="greater");

count.amp.d <- rbind(
	data.frame(count.amp, group = "other"),
	data.frame(count.amp.drivers, group ="driver")
);

lapply(rownames(count.amp.drivers),
	function(x) {
		x <- count.amp.drivers[x, "Brain metastasis"];
		poisson.test(x, r = mean(count.amp[, "Brain metastasis"]), alternative="greater")
	}
);
# YAP1 is not significant (only 2 events...)

fit <- lm(cbind(Brain.metastasis, Primary, Shared) ~ group, data=count.amp.d);
summary(fit)
# only brain metastasis counts are significantly different between groups

hs.amp.glm <- lapply(1:ncol(count.amp),
	function(v) {
		fit <- glm(count.amp.d[[v]] ~ count.amp.d$group, family="quasipoisson");
		z <- coef(summary(fit))[2, ];
		list(statistic = unname(z[3]), p.value = unname(z[4]), coefficients=coef(summary(fit)))
	}
);
names(hs.amp.glm) <- colnames(count.amp);


####

gmarginal.del <- apply(count.del, 1, sum);
gmarginal.del.drivers <- apply(count.del.drivers, 1, sum);

mean(gmarginal.del)
fit <- glm(gmarginal.del ~ 1, family="quasipoisson");
exp(coef(fit))

lapply(gmarginal.del.drivers,
	function(x) {
		poisson.test(x, r = mean(gmarginal.del), alternative="greater")
	}
);

h.del.pois <- poisson.test(mean(gmarginal.del.drivers), r = mean(gmarginal.del), alternative="greater");

count.del.d <- rbind(
	data.frame(count.del, group = "other"),
	data.frame(count.del.drivers["CDKN2B", , drop=FALSE], group ="driver")
);

fit <- lm(cbind(Brain.metastasis, Primary, Shared) ~ group, data=count.del.d);
summary(fit)
# brain metastasis and shared counts are significantly different between groups


lapply(rownames(count.del.drivers),
	function(x) {
		x <- count.del.drivers[x, "Brain metastasis"];
		poisson.test(x, r = mean(count.del[, "Brain metastasis"]), alternative="greater")
	}
);

fit <- glm(Brain.metastasis ~ group, data=count.del.d, family="quasipoisson");
z <- coef(summary(fit))[2, ];
h.del.pois <- list(statistic = z[3], p.value = z[4]);

hs.del.glm <- lapply(1:ncol(count.del),
	function(v) {
		fit <- glm(count.del.d[[v]] ~ count.del.d$group, family="quasipoisson");
		z <- coef(summary(fit))[2, ];
		list(statistic = unname(z[3]), p.value = unname(z[4]), coefficients=coef(summary(fit)))
	}
);
names(hs.del.glm) <- colnames(count.del);


###

qdraw(
	ggplot(compart.amp.drivers, aes(x=compartment, y=count, fill=gene)) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		theme(
			legend.text = element_text(face = "bold.italic"),
			axis.text.x = element_text(angle=30, hjust=1)
		) +
		xlab("") + ylab("# patients")
	,
	width = 2.4, height=3,
	file = insert(pdf.fname, c("comparts", "bar"))
);

####

qdraw(
	ggplot(compart.amp2, aes(x=gene, y=count, fill=branch)) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		theme(
			axis.text.x = element_text(angle=30, hjust=1, face = "bold.italic")
		) +
		ylim(0, 11) +
		annotate(x = 0.5, y = 11, geom="text", label = test_to_string(h.amp.pois), hjust=0) +
		annotate(x = 0.5, y = 11*0.75, geom="text", label = tests_to_string(hs.amp.glm), hjust=0) +
		xlab("") + ylab("# patients")
	,
	width = 3.4, height = 3,
	file = insert(pdf.fname, c("amp-genes", "bar"))
);

compart.del2.f <- filter(compart.del2, gene != "CDKN2A");

qdraw(
	ggplot(compart.del2.f, aes(x=gene, y=count, fill=branch)) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		theme(
			axis.text.x = element_text(angle=30, hjust=1, face = "bold.italic")
		) +
		ylim(0, 20) +
		annotate(x = 0.5, y = 20, geom="text", label = test_to_string(h.del.pois), hjust=0) +
		annotate(x = 0.5, y = 20*0.75, geom="text", label = tests_to_string(hs.del.glm), hjust=0) +
		xlab("") + ylab("# patients")
	,
	width = 2.8, height = 3,
	file = insert(pdf.fname, c("del-genes", "bar"))
);

####

marginal.amp.drivers <- apply(count.amp.drivers, 2, sum);
marginal.amp <- apply(count.amp, 2, sum);
ct.amp <- matrix(c(marginal.amp.drivers, marginal.amp), ncol=2,
	dimnames = list(compartment = names(marginal.amp), group = c("Driver genes", "Other genes"))
);

print(ct.amp)
fisher.test(ct.amp)

compart.levels.sh.pr <- c("Shared", "Primary");
compart.levels.bm.sh <- c("Brain metastasis", "Shared");

# event in the primary tumor is inherited by brain metastasis
h.amp.sh.pr <- fisher.test(ct.amp[compart.levels.sh.pr, ]);

# late vs. early event
h.amp.bm.sh <- fisher.test(ct.amp[compart.levels.bm.sh, ]);


marginal.del.drivers <- apply(count.del.drivers, 2, sum);
marginal.del <- apply(count.del, 2, sum);
ct.del <- matrix(c(marginal.del.drivers, marginal.del), ncol=2,
	dimnames = list(compartment = names(marginal.del), group = c("Driver genes", "Other genes"))
);

print(ct.del)
fisher.test(ct.del)

# event in the primary tumor is inherited by brain metastasis
h.del.sh.pr <- fisher.test(ct.del[compart.levels.sh.pr, ]);

# late vs. early event
h.del.bm.sh <- fisher.test(ct.del[compart.levels.bm.sh, ]);

####

collapsed.amp <- group_by(compart.amp, compartment, group) %>% summarize(count = sum(count));

collapsed.amp.sh.pr <- filter(collapsed.amp, compartment %in% compart.levels.sh.pr);
total.amp.sh.pr <- collapsed.amp.sh.pr %>% group_by(group) %>% summarize(total = sum(count));
prop.amp.sh.pr <- left_join(collapsed.amp.sh.pr, total.amp.sh.pr);
bi <- with(prop.amp.sh.pr, binom.agresti.coull(count, total, conf.level = 1 - alpha.level));
prop.amp.sh.pr <- cbind(as.data.frame(prop.amp.sh.pr), as.data.frame(bi[, -(1:3)]));

collapsed.amp.bm.sh <- filter(collapsed.amp, compartment %in% compart.levels.bm.sh);
total.amp.bm.sh <- collapsed.amp.bm.sh %>% group_by(group) %>% summarize(total = sum(count));
prop.amp.bm.sh <- left_join(collapsed.amp.bm.sh, total.amp.bm.sh);
bi <- with(prop.amp.bm.sh, binom.agresti.coull(count, total, conf.level = 1 - alpha.level));
prop.amp.bm.sh <- cbind(as.data.frame(prop.amp.bm.sh), as.data.frame(bi[, -(1:3)]));


qdraw(
	ggplot(filter(prop.amp.sh.pr, compartment == "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper, fill=group)) +
		geom_bar(stat="identity") + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		scale_fill_manual(values=c("firebrick", "grey60")) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.amp.sh.pr), hjust=0) +
		ylab("P. of CNAs in primary shared by matched BM") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("amp-genes", "prop-sh-pr"))
);

qdraw(
	ggplot(filter(prop.amp.bm.sh, compartment != "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper, fill=group)) +
		geom_bar(stat="identity") + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		scale_fill_manual(values=c("firebrick", "grey60")) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.amp.bm.sh), hjust=0) +
		ylab("P. of CNAs in BM depleted from matched primary") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("amp-genes", "prop-bm-sh"))
);

qdraw(
	ggplot(filter(prop.amp.sh.pr, compartment == "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.amp.sh.pr), hjust=0) +
		ylab("P. of CNAs in primary shared by matched BM") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("amp-genes", "prop-sh-pr", "vdot"))
);

qdraw(
	ggplot(filter(prop.amp.bm.sh, compartment != "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.amp.bm.sh), hjust=0) +
		ylab("P. of CNAs in BM depleted from matched primary") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("amp-genes", "prop-bm-sh", "vdot"))
);

qdraw(
	ggplot(filter(prop.amp.sh.pr, compartment == "Shared"), aes(x = factor(group, rev(unique(group))), y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_dot() + coord_flip() +
		theme(
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 2, y = 0, geom="text", label = test_to_string(h.amp.sh.pr), hjust=0) +
		ylab("P. of CNAs in primary shared by matched BM") + xlab("")
	,
	width = 3, height = 1.2,
	file = insert(pdf.fname, c("amp-genes", "prop-sh-pr", "dot"))
);

qdraw(
	ggplot(filter(prop.amp.bm.sh, compartment != "Shared"), aes(x = factor(group, rev(unique(group))), y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_dot() + coord_flip() +
		theme(
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 1, y = 0.5, geom="text", label = test_to_string(h.amp.bm.sh), hjust=0) +
		ylab("P. of CNAs in BM depleted from matched primary") + xlab("")
	,
	width = 3, height = 1.2,
	file = insert(pdf.fname, c("amp-genes", "prop-bm-sh", "dot"))
);

####

collapsed.del <- group_by(filter(compart.del, is.na(gene) | gene != "CDKN2A"), compartment, group) %>% summarize(count = sum(count));

collapsed.del.sh.pr <- filter(collapsed.del, compartment %in% compart.levels.sh.pr);
total.del.sh.pr <- collapsed.del.sh.pr %>% group_by(group) %>% summarize(total = sum(count));
prop.del.sh.pr <- left_join(collapsed.del.sh.pr, total.del.sh.pr);
bi <- with(prop.del.sh.pr, binom.agresti.coull(count, total, conf.level = 1 - alpha.level));
prop.del.sh.pr <- cbind(as.data.frame(prop.del.sh.pr), as.data.frame(bi[, -(1:3)]));
prop.del.sh.pr$lower <- pmax(prop.del.sh.pr$lower, 0);
prop.del.sh.pr$upper <- pmin(prop.del.sh.pr$upper, 1);

collapsed.del.bm.sh <- filter(collapsed.del, compartment %in% compart.levels.bm.sh);
total.del.bm.sh <- collapsed.del.bm.sh %>% group_by(group) %>% summarize(total = sum(count));
prop.del.bm.sh <- left_join(collapsed.del.bm.sh, total.del.bm.sh);
bi <- with(prop.del.bm.sh, binom.agresti.coull(count, total, conf.level = 1 - alpha.level));
prop.del.bm.sh <- cbind(as.data.frame(prop.del.bm.sh), as.data.frame(bi[, -(1:3)]));

qdraw(
	ggplot(filter(prop.del.sh.pr, compartment == "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper, fill=group)) +
		geom_bar(stat="identity") + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		scale_fill_manual(values=c("royalblue3", "grey60")) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.del.sh.pr), hjust=0) +
		ylab("P. of CNAs in primary shared by matched BM") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("del-genes", "prop-sh-pr"))
);

qdraw(
	ggplot(filter(prop.del.bm.sh, compartment != "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper, fill=group)) +
		geom_bar(stat="identity") + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		scale_fill_manual(values=c("royalblue3", "grey60")) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.del.bm.sh), hjust=0) +
		ylab("P. of CNAs in BM depleted from matched primary") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("del-genes", "prop-bm-sh"))
);

qdraw(
	ggplot(filter(prop.del.sh.pr, compartment == "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.del.sh.pr), hjust=0) +
		ylab("P. of CNAs in primary shared by matched BM") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("del-genes", "prop-sh-pr", "vdot"))
);

qdraw(
	ggplot(filter(prop.del.bm.sh, compartment != "Shared"), aes(x = group, y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_bw() +
		theme(
			axis.text.x = element_text(angle=30, hjust=1),
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 0.5, y = 1, geom="text", label = test_to_string(h.del.bm.sh), hjust=0) +
		ylab("P. of CNAs in BM depleted from matched primary") + xlab("")
	,
	width = 1.4, height = 3,
	file = insert(pdf.fname, c("del-genes", "prop-bm-sh", "vdot"))
);

qdraw(
	ggplot(filter(prop.del.sh.pr, compartment == "Shared"), aes(x = factor(group, rev(unique(group))), y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_dot() + coord_flip() +
		theme(
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 2, y = 0, geom="text", label = test_to_string(h.del.sh.pr), hjust=0) +
		ylab("P. of CNAs in primary shared by matched BM") + xlab("")
	,
	width = 3, height = 1.2,
	file = insert(pdf.fname, c("del-genes", "prop-sh-pr", "dot"))
);

qdraw(
	ggplot(filter(prop.del.bm.sh, compartment != "Shared"), aes(x = factor(group, rev(unique(group))), y = mean, ymin = lower, ymax = upper)) +
		geom_point() + geom_errorbar(width=0.1) + theme_dot() + coord_flip() +
		theme(
			panel.grid = element_blank(),
			legend.position = "none"
		) +
		ylim(0, 1) +
		annotate(x = 2, y = 0.5, geom="text", label = test_to_string(h.del.bm.sh), hjust=0) +
		ylab("P. of CNAs in BM depleted from matched primary") + xlab("")
	,
	width = 3, height = 1.2,
	file = insert(pdf.fname, c("del-genes", "prop-bm-sh", "dot"))
);

qwrite(prop.amp.bm.sh, insert(tsv.fname, tag=c("amp", "prop-bm-sh")));
qwrite(prop.amp.sh.pr, insert(tsv.fname, tag=c("amp", "prop-sh-pr")));
qwrite(prop.del.bm.sh, insert(tsv.fname, tag=c("del", "prop-bm-sh")));
qwrite(prop.del.sh.pr, insert(tsv.fname, tag=c("del", "prop-sh-pr")));


# Statistics important for power analysis

hs.amp.glm[["Brain metastasis"]]$coefficients
exp(hs.amp.glm[["Brain metastasis"]]$coefficients[1,1])
exp(hs.amp.glm[["Brain metastasis"]]$coefficients[2,1])
exp(sum(hs.amp.glm[["Brain metastasis"]]$coefficients[,1]))
table(count.amp.d$group)
prop.amp.bm.sh
count.amp.drivers
count.amp.drivers / gmarginal.amp.drivers
gmarginal.amp.drivers
gmarginal.amp.drivers / 57
summary(glm(count.amp.d[["Brain.metastasis"]] ~ count.amp.d$group, family="quasipoisson"))

hs.del.glm[["Brain metastasis"]]$coefficients
exp(hs.del.glm[["Brain metastasis"]]$coefficients[1,1])
exp(hs.del.glm[["Brain metastasis"]]$coefficients[2,1])
exp(sum(hs.del.glm[["Brain metastasis"]]$coefficients[,1]))
table(count.del.d$group)
marginal.del <- apply(count.del, 2, sum);
prop.del.bm.sh
count.del.drivers
count.del.drivers / gmarginal.del.drivers
gmarginal.del.drivers
gmarginal.del.drivers / 57
summary(glm(count.del.d[["Brain.metastasis"]] ~ count.del.d$group, family="quasipoisson"))

