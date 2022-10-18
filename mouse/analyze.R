library(io);
library(dplyr);
library(reshape2);
library(ggplot2);
library(binom);

x <- qread("xenograft_2018-06-02.tsv");

out.fname <- filename("brain-met-incidence", tag=c("mouse-intracardiac-inject"));

source("~/exigens/brain-mets/common/params.R")


# analyze experiment C by itself

y.c <- filter(x, experiment == "C") %>% select(-experiment, -time_days);
ct.c <- as.matrix(y.c[, c("brain_met_positive", "brain_met_negative")]);
rownames(ct.c) <- y.c$treatment;

ct.c

fisher.test(ct.c)

fisher.test(ct.c[c("LacZ", "MMP13"), ]);
fisher.test(ct.c[c("LacZ", "MYC"), ]);
fisher.test(ct.c[c("LacZ", "YAP1"), ])


# combine experiments

y <- select(x, -experiment) %>% 
	group_by(treatment, time_days) %>%
	summarize(
		brain_met_positive=sum(brain_met_positive),
		brain_met_negative=sum(brain_met_negative)
	);

# create contingency table
y.24 <- filter(y, time_days == 24) %>% select(-time_days);
ct.24 <- as.matrix(y.24[, c("brain_met_negative", "brain_met_positive")]);
rownames(ct.24) <- y.24$treatment;

fisher.test(ct.24)

fisher.test(ct.24[c("LacZ", "MMP13"), ]);
fisher.test(ct.24[c("LacZ", "MYC"), ]);
fisher.test(ct.24[c("LacZ", "YAP1"), ])


# create contingency table
y.12 <- filter(y, time_days == 12) %>% select(-time_days);
ct.12 <- as.matrix(y.12[, c("brain_met_negative", "brain_met_positive")]);
rownames(ct.12) <- y.12$treatment;

ct.12

fisher.test(ct.12)

fisher.test(ct.12, alternative="greater")

ct.12[c("LacZ", "MMP13"), ]
fisher.test(ct.12[c("LacZ", "MMP13"), ]);
fisher.test(ct.12[c("LacZ", "MYC"), ]);
fisher.test(ct.12[c("LacZ", "YAP1"), ])

fisher.test(ct.12[c("LacZ", "MMP13"), ], alternative="greater")
fisher.test(ct.12[c("LacZ", "MYC"), ], alternative="greater")
fisher.test(ct.12[c("LacZ", "YAP1"), ], alternative="greater")


y.e2.e3 <- filter(x, experiment %in% c("B", "C")) %>% 
	select(-experiment) %>% 
	group_by(treatment, time_days) %>%
	summarize(
		brain_met_positive=sum(brain_met_positive),
		brain_met_negative=sum(brain_met_negative)
	);
y.d12.e2.e3 <- filter(y.e2.e3, time_days == 12) %>% select(-time_days);
ct.d12.e2.e3 <- as.matrix(y.d12.e2.e3[, c("brain_met_negative", "brain_met_positive")]);
rownames(ct.d12.e2.e3) <- y.d12.e2.e3$treatment;

fisher.test(ct.d12.e2.e3)

ct.d12.e2.e3[c("LacZ", "MMP13"), ]
fisher.test(ct.d12.e2.e3[c("LacZ", "MMP13"), ]);
ct.d12.e2.e3[c("LacZ", "MYC"), ]
fisher.test(ct.d12.e2.e3[c("LacZ", "MYC"), ]);
fisher.test(ct.d12.e2.e3[c("LacZ", "YAP1"), ])


# prepare data for glm
d <- melt(lapply(
	as.list(as.data.frame(t(ct.12))),
	function(counts) {
		c(rep(1, counts[1]), rep(0, counts[2]))
	}
));
colnames(d) <- c("brain_met", "treatment");
d$treatment <- factor(d$treatment, levels=rownames(ct.12));

table(d$treatment, d$brain_met)

# assess treatment effect with log likelihood ratio test

fit0 <- glm(brain_met ~ 1, data=d, family=binomial);
fit <- glm(brain_met ~ treatment, data=d, family=binomial);
anova(fit0, fit, test="LRT")

# standard error estimate is wide in control group,
# because no event is observed
summary(fit)


# prepare data for glm
d2 <- mapply(
	function(treatment, time, pos, neg) {
		data.frame(
			treatment,
			time,	
			event = c(rep(1, pos), 	rep(0, neg))
		)
	},
	y$treatment,
	paste0("d", y$time_days),
	y$brain_met_positive,
	y$brain_met_negative,
	SIMPLIFY=FALSE
);
d2 <- do.call(rbind, d2);

fit.null <- glm(event ~ 1, data=d2, family=binomial);
fit.treatment <- glm(event ~ treatment, data=d2, family=binomial);
fit.time <- glm(event ~ time, data=d2, family=binomial);
fit <- glm(event ~ treatment + time, data=d2, family=binomial);
fit.full <- glm(event ~ treatment + time + treatment:time, data=d2, family=binomial);
anova(fit.null, fit.time, fit, fit.full, test="LRT")
anova(fit.null, fit.treatment, fit, fit.full, test="LRT")

summary(fit)

# difficult to estimate coefficient in control group because
# no events occurred in the control group
summary(fit.full)

# use d24, LacZ group as the control group
d3 <- d2;
d3$time <- relevel(d2$time, "d24");
fit.full <- glm(event ~ treatment + time + treatment:time, data=d3, family=binomial);
summary(fit.full)


# create contingency table but combine day 12 and day 24, since treatment
# effect persist over time, and time effect is not significant
y.c <- select(x, -experiment) %>% 
	group_by(treatment) %>%
	summarize(
		brain_met_positive=sum(brain_met_positive),
		brain_met_negative=sum(brain_met_negative)
	);
ct.c <- as.matrix(y.c[, c("brain_met_negative", "brain_met_positive")]);
rownames(ct.c) <- y.c$treatment;

ct.c

fisher.test(ct.c)

ct.c[c("LacZ", "MMP13"), ]
fisher.test(ct.c[c("LacZ", "MMP13"), ]);
fisher.test(ct.c[c("LacZ", "MYC"), ]);
fisher.test(ct.c[c("LacZ", "YAP1"), ])

fisher.test(ct.c[c("LacZ", "MMP13"), ], alternative="greater")
fisher.test(ct.c[c("LacZ", "MYC"), ], alternative="greater")
fisher.test(ct.c[c("LacZ", "YAP1"), ], alternative="greater")


binom_confint <- function(x, n, conf.level) {
	y <- binom.confint(x, n, conf.level=conf.level, method="exact");
	y[, c("mean", "lower", "upper")]
}


#y.use <- y.d12.e2.e3;
#ct.use <- ct.d12.e2.e3;

y.use <- y.12;
ct.use <- ct.12;

risks <- binom_confint(y.use$brain_met_positive, y.use$brain_met_negative, conf.level=0.80)
risks <- data.frame(treatment = y.use$treatment, risks);

treatments <- c("LacZ", "MYC", "MMP13", "YAP1");

ps <- vapply(treatments[-1],
	function(treatment) {
		#fisher.test(ct.use[c("LacZ", treatment), ])$p.value
		fisher.test(ct.use[c("LacZ", treatment), ], alternative="greater")$p.value
	},
	0
);

p_to_str <- function(x) {
	sprintf("p = %s", format(x, digits=1))
}

ps.str <- p_to_str(ps);

ct.use <- ct.use[treatments, ];
ns <- apply(ct.use, 1, sum);

risks$n <- ns;
risks$treatment_label <- factor(risks$treatment, levels=treatments, labels=paste0(treatments, "\n(", ns, ")"));

tickh <- 2;
ynudge <- 1.5;
ybase <- 45;
yspace <- 8;

yhbar <- ybase + yspace * 0:2;

qdraw(
	ggplot(risks, aes(x=treatment_label, y=mean*100, ymin=lower*100, ymax=upper*100)) +
		geom_bar(stat="identity", fill="orange") + theme_bw() +
		geom_errorbar(width=0.2, colour="grey30") +
		xlab("") + ylab("Brain metastasis incidence (%)") +
		ylim(0, max(yhbar) + yspace) +
		theme(legend.position="none", panel.grid=element_blank()) +
		geom_path(x=c(1, 1, 2, 2), y=c(yhbar[1]-tickh, yhbar[1], yhbar[1], yhbar[1]-tickh)) +
		geom_path(x=c(1, 1, 3, 3), y=c(yhbar[2]-tickh, yhbar[2], yhbar[2],  yhbar[2]-tickh)) +
		geom_path(x=c(1, 1, 4, 4), y=c(yhbar[3]-tickh, yhbar[3], yhbar[3], yhbar[3]-tickh)) +
		geom_text(x=mean(c(1, 2)), y=yhbar[1] + ynudge, label=ps.str[[1]], vjust=0) +
		geom_text(x=mean(c(1, 3)), y=yhbar[2] + ynudge, label=ps.str[[2]], vjust=0) +
		geom_text(x=mean(c(1, 4)), y=yhbar[3] + ynudge, label=ps.str[[3]], vjust=0)
	,
	width = 3,
	height = 4,
	file = insert(out.fname, tag=c("d12", "all"), ext="pdf")
);

risks$group <- "A";
risks$group[risks$treatment != "LacZ"] <- "B";

qdraw(
	ggplot(risks, aes(x=treatment_label, y=mean*100, ymin=lower*100, ymax=upper*100)) +
		geom_errorbar(width=0.2, colour="grey30") +
		#geom_point(aes(colour=group)) + theme_bw() +
		geom_point(colour="darkorange") + theme_bw() +
		scale_colour_manual(values=col.cohorts) + 
		theme(legend.position="none", panel.grid=element_blank()) +
		xlab("") + ylab("Brain metastasis incidence (%)") +
		ylim(0, max(yhbar) + yspace) +
		geom_path(x=c(1, 1, 2, 2), y=c(yhbar[1]-tickh, yhbar[1], yhbar[1], yhbar[1]-tickh)) +
		geom_path(x=c(1, 1, 3, 3), y=c(yhbar[2]-tickh, yhbar[2], yhbar[2],  yhbar[2]-tickh)) +
		geom_path(x=c(1, 1, 4, 4), y=c(yhbar[3]-tickh, yhbar[3], yhbar[3], yhbar[3]-tickh)) +
		geom_text(x=mean(c(1, 2)), y=yhbar[1] + ynudge, label=ps.str[[1]], vjust=0) +
		geom_text(x=mean(c(1, 3)), y=yhbar[2] + ynudge, label=ps.str[[2]], vjust=0) +
		geom_text(x=mean(c(1, 4)), y=yhbar[3] + ynudge, label=ps.str[[3]], vjust=0)
	,
	width = 3,
	height = 4,
	file = insert(out.fname, tag=c("d12", "all", "dot"), ext="pdf")
);


qwrite(select(risks, -treatment_label), insert(out.fname, c("d12", "all"), ext="tsv"));

