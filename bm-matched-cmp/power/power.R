library(io)
library(ggplot2)
library(powerMediation)
library(pwr)

paired_analysis_main_effect <- function(prevalence, p.exclusive, beta0, n) {
	f <- prevalence * p.exclusive * n;

	# main effect
	# log(2.5) -- log(3.5)
	# f == exp(beta1 - beta0)
	# log(f) == beta1 - beta0
	# beta1 == log(f) + beta0
	log(f) + beta0
}

# Power calculation for paired analysis.
#
# Test alternative hypothesis of beta1 > 0
#
# @param prevalence   frequency of event of interest across paired samples
# @param p.exclusive       probability of event occuring exclusively in the target compartment
#                     (e.g. metastasis, primary, or shared)
# @param beta0        rate of background events
# @param n            number of pairs; sample size
# @param phi          overdispersion parameter
# @param alpha        significance level
# @return power (numeric)
paired_analysis_power <- function(prevalence, p.exclusive, beta0, n, phi=1, alpha=0.05) {
	beta1 <- paired_analysis_main_effect(prevalence, p.exclusive, beta0, n);

	# number of recurrently aberrant gene
	# recurrently amplified: 2485 + 3
	# recurrently deleted: 299 + 1
	# since the power calculated below is insensitity to this parameter (100 - 10000)
	# we fix this parameter to 1000
	ngenes <- 1000;

	# binary predictor: driver gene vs. non-driver genes
	# Bernoulii distributed variable
	mu.x <- 1 / ngenes;
	sigma2.x <- mu.x * (1 - mu.x);

	ifelse(beta1 <= 0,
		NA,
		# double the significance level, since
		# the function considers a two-tailed test
		powerPoisson(beta0, beta1, mu.x, sigma2.x,
			alpha = alpha*2, phi = phi, N = ngenes)
	)
}

####

# significance level
alpha <- 0.05;

# nubmer of pairs
n <- 57;

# overdispersion parameter
phi <- 1;

# observed values for amplified drivers (YAP1, MMP13, and MYC)
# prevalence: 0.15
# p.exclusive (probability of event occuring only in the brain metastasis): 0.4
# beta0 (rate of event occuring only in the BM for non-driver genes): 1.0
message(
	"realized power for amplified drivers: ",
	paired_analysis_power(prevalence=0.15, p.exclusive=0.4, beta0=log(1.0), n=n, phi=phi, alpha=alpha)
)

# observed values for amplified drivers (CDKN2A, CDKN2B)
# prevalence: 0.30
# p.exclusive (probability of event occuring only in the brain metastasis): 0.2
# beta0 (rate of event occuring only in the BM for non-driver genes): 0.9
message(
	"realized power for deleted drivers: ",
	paired_analysis_power(prevalence=0.30, p.exclusive=0.2, beta0=log(0.9), n=n, phi=phi, alpha=alpha)
)

####

# significance level
alpha <- 0.05;

# overdispersion parameter
phi <- 1;

prevalences <- c(0.05, 0.1, 0.2, 0.4);
ps.exclusive <- c(0.1, 0.2, 0.4, 0.8);
ns <- seq(10, 500, by=10);
beta0 <- log(1.0);

pa.d <- expand.grid(prevalences, ps.exclusive, beta0, ns);
names(pa.d) <- c("prevalence", "p.exclusive", "beta0", "n");

pa.d$power <- with(pa.d,
	paired_analysis_power(
		prevalence=prevalence, p.exclusive=p.exclusive, beta0=beta0, n=n, phi=phi, alpha=alpha
	)
);

p.exclusive.labeller <- function(value) {
	sprintf("P(late) = %s", value)
}

obs.n.pairs <- 57;

qdraw(
	ggplot(pa.d, aes(x=n, y=power, colour=factor(prevalence))) +
		geom_hline(yintercept=0.8, colour="black") +
		geom_vline(xintercept=obs.n.pairs, colour="grey60") +
		geom_line(size=1) + theme_bw() +
		facet_wrap(~ p.exclusive,
			strip.position="right",
			labeller = labeller(p.exclusive = p.exclusive.labeller)
		) +
		scale_colour_brewer(palette="Set2") +
		# requires a minimum of 3 samples per patient: normal, primary, and metastasis
		scale_x_continuous(sec.axis = sec_axis(~ . * 3, name = "Number of samples")) +
		scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 1)) +
		theme(
			strip.background = element_blank(),
			panel.grid = element_blank()
		) +
		labs(colour = "Driver frequency") +
		ylab("Power to detect driver") +
		xlab("Number of patients")
	,
	width = 6.5, height = 4.8,
	file = "power_paired-analysis.pdf"
);


case_control_analysis_power <- function(p1, p2, n1, n2, g, alpha=0.05) {
	pwr.2p2n.test(
		h = ES.h(p1 = p1, p2 = p1*g + p2*(1-g)),
		n1 = n1,
		n2 = n2,
		sig.level = alpha,
		alternative = "greater"
	)$power
}

# observed values uing TCGA as the control cohort
g <- 0.3;
n2 <- 464;
p2 <- 0.01;

cca.d <- expand.grid(prevalences, p2, ns, n2, g);
names(cca.d) <- c("p1", "p2", "n1", "n2", "g");
cca.d$power <- with(cca.d, case_control_analysis_power(p1, p2, n1, n2, g, alpha));

obs.n.cases <- 73;

qdraw(
	ggplot(cca.d, aes(x=n1, y=power, colour=factor(p1))) +
		geom_hline(yintercept=0.8, colour="black") +
		geom_vline(xintercept=obs.n.cases, colour="grey60") +
		geom_line(size=1) + theme_bw() +
		scale_colour_brewer(palette="Set2") +
		# requires a minimum of 2 samples per patient: normal and tumour
		scale_x_continuous(sec.axis = sec_axis(~ . * 2, name = "Number of samples")) +
		scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 1)) +
		theme(
			strip.background = element_blank(),
			panel.grid = element_blank()
		) +
		labs(colour = "Driver frequency") +
		ylab("Power to detect driver") +
		xlab("Number of patients")
	,
	width = 4, height = 3,
	file = "power_case-control-analysis.pdf"
);


case_control_analysis_power(p1 = 0.2, p2 = 0.01, n1 = 73, n2 = 464, g = 0.3, alpha = 0.025)
case_control_analysis_power(p1 = 0.2, p2 = 0.01, n1 = 73, n2 = 464, g = 0.3, alpha = 0.05)

case_control_analysis_power(p1 = 0.1, p2 = 0.01, n1 = 73, n2 = 464, g = 0.3, alpha = 0.025)
case_control_analysis_power(p1 = 0.1, p2 = 0.01, n1 = 73, n2 = 464, g = 0.3, alpha = 0.05)

# MYC, MMP13, YAP1 for BM-LUAD vs. TCGA-LUAD
case_control_analysis_power(p1 = 0.14, p2 = 0.01, n1 = 73, n2 = 464, g = 0.3, alpha = 0.025)

# MYC, MMP13, YAP1 for BM-LUAD-V vs. Tracerx
case_control_analysis_power(p1 = 0.14, p2 = 0.01, n1 = 205, n2 = 61, g = 0.3, alpha = 0.025)

# MYC, MMP13, YAP1 for BM-LUAD-V vs. TCGA-LUAD
case_control_analysis_power(p1 = 0.14, p2 = 0.01, n1 = 205, n2 = 464, g = 0.3, alpha = 0.025)

