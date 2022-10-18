library(io)
library(ggplot2)
library(pwr)

# case contamination of control cohort
gs <- seq(0, 1, by=0.01);

p1 <- 0.2;
p2 <- 0.01;
alpha <- 0.05;
n1 <- 73;
n2 <- 464;

p1s <- c(0.05, 0.1, 0.2, 0.4);

powers <- lapply(p1s,
	function(p1) {
		unlist(lapply(gs,
			function(g) {
				pwr.2p2n.test(
					h = ES.h(p1 = p1, p2 = p1*g + p2*(1-g)),
					n1 = n1,
					n2 = n2,
					sig.level = alpha,
					alternative = "greater"
				)$power
			}
		))
	}
);

d <- do.call(rbind,
	mapply(
		function(p1, p) {
			data.frame(p1 = p1, contamination = gs, power = p)
		},
		p1s,
		powers,
		SIMPLIFY=FALSE
	)
);

qdraw(
	ggplot(d, aes(x=contamination, y=power, colour=factor(p1))) +
		annotate("rect", xmin=0.1, xmax=0.61, ymin=-Inf, ymax=Inf, fill="black", alpha=0.1) +
		geom_line(size=1) + theme_bw() + 
		geom_vline(xintercept=0.3, size=1, colour="grey60") +
		scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 1)) +
		scale_x_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 1)) +
		scale_colour_brewer(palette="Set2") +
		theme(panel.grid = element_blank()) +
		xlab("Fraction of brain metastasis patients in TCGA-LUAD") +
		ylab("Power to detect increase in driver frequency") +
		labs(colour = "Driver frequency")
	,
	width = 6,
	file = "power-vs-contamination.pdf"
);

sim_ct <- function(m, n1, n2, p1, p2) {
	x1 <- rbinom(m, n1, p1);
	x2 <- rbinom(m, n2, p2);

	data.frame(
		x1,
		y1 = n1 - x1,
		x2,
		y2 = n2 - x2
	)
}

test <- function(z) {
	#fisher.test(matrix(z, nrow=2), alternative="greater")$estimate
	fisher.test(matrix(z, nrow=2), alternative="greater")$p.value
}

m <- 100;

# cache simulation results
if (!file.exists("fprs.rds")) {
	fprs <- lapply(p1s,
		function(p1) {
			unlist(lapply(gs,
				function(g) {

					# since the null hypothesis is true, p1 == p2
					# p1 and p2 cannot be too low;
					# otherwise, low n1 causes significant p-value to more difficult to achieve
					p2 <- p1;
					p2c <- p1*g + p2*(1-g);

					cts <- sim_ct(m, n1, n2, p1, p2c);
					results <- apply(cts, 1, test);
					mean(results < alpha)
				}
			))
		}
	);
	qwrite(fprs, "fprs.rds");
} else {
	fprs <- qread("fprs.rds");
}

d1 <- do.call(rbind,
	mapply(
		function(p1, f) {
			data.frame(p1 = p1, contamination = gs, fpr = f)
		},
		p1s,
		fprs,
		SIMPLIFY=FALSE
	)
);

qdraw(
	ggplot(d1, aes(x=contamination, y=fpr, colour=factor(p1))) +
		annotate("rect", xmin=0.1, xmax=0.61, ymin=-Inf, ymax=Inf, fill="black", alpha=0.1) +
		geom_line(size=1) + theme_bw() + 
		theme(panel.grid = element_blank()) +
		geom_vline(xintercept=0.3, size=1, colour="grey60") +
		scale_x_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 1)) +
		scale_y_continuous(breaks=c(0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0), limits=c(0, 1)) +
		geom_hline(yintercept=alpha, colour="black") +
		scale_colour_brewer(palette="Set2") +
		xlab("Fraction of brain metastasis patients in TCGA-LUAD") +
		ylab("False positive rate") +
		labs(colour = "Driver frequency")
	,
	width = 6,
	file = "fpr-vs-contamination.pdf"
);

