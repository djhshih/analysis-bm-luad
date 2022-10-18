library(io);
library(danitizer);
library(dplyr);
library(ggplot2);

control.fname <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_absolute/scores.gistic";
case.fname <- "../luad_absolute_bmet-only/scores.gistic";

s.control <- qread(control.fname, type = "tsv") %>% normalize_df_field_names();
s.case <- qread(case.fname, type="tsv") %>% normalize_df_field_names();

d <- data.frame(
	g_score = c(s.control$g_score, s.case$g_score),
	cna_type = factor(c(s.control$type, s.case$type), labels=c("Amp", "Del")),
	group = c(rep("control", nrow(s.control)), rep("case", nrow(s.case)))
);

# compare distribution of G scores
qdraw(
	ggplot(d, aes(x = g_score)) + geom_histogram(bins=100) + facet_grid(group ~ cna_type) + theme_bw(),
	file = "gscore-hist_bm-luad_tcga-luad.pdf"
);


chrom <- 8;
cna.type <- "Amp";
pos <- 48100000;

s1 <- filter(s.control, chromosome == chrom, end >= pos, type == cna.type);
s2 <- filter(s.case, chromosome == chrom, end >= pos, type == cna.type);

# a very small proportion of overlapping positions between the two groups
length(intersect(s1$start, s2$start)) / length(union(s1$start, s2$start))
length(intersect(s1$end, s2$end)) / length(union(s1$end, s2$end))

add_noise <- function(x, ...) {
	x + rnorm(length(x), ...)
}

# add noise to position to make them unique
#s1$position <- add_noise(s1$position);
#s2$position <- add_noise(s2$position);

#stopifnot(length(intersect(s1$position, s2$position)) == 0)

data <- list(
	J = 2*nrow(s1) + 2*nrow(s2),
	x = c(s1$start, s1$end, s2$start, s2$end) / 10e6,
	g = c(rep(0, 2*nrow(s1)), rep(1, 2*nrow(s2))),
	y = c(s1$g_score, s1$g_score, s2$g_score, s2$g_score)
);

# order by x position
idx <- order(data$x);
data$x <- data$x[idx];
data$g <- data$g[idx];
data$y <- data$y[idx];

qdraw(
	with(data,
		{
			plot(NA, xlim=range(x), ylim=range(y));
			lines(x[g == 0], y[g == 0], col="grey", type="b", pch=20);
			lines(x[g == 1], y[g == 1], col="orange", type="b", pch=20);
		}
	),
	width = 10,
	file = "gscore_chr8q-amp_bm-luad_tcga-luad_ends.pdf"
);


qwrite(data, "gistic-chr8q-amp_bm-luad_tcga-luad_ends.rds");

