library(io);
library(danitizer);
library(dplyr);
library(ggplot2);

control.fname <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_absolute/scores.gistic";
case.fname <- "../luad_absolute_bmet-only/scores.gistic";

region_center <- function(start, end) {
	start + floor((end - start)/2)
}

s.control <- qread(control.fname, type = "tsv") %>% normalize_df_field_names() %>%
	mutate(position = region_center(start, end))

s.case <- qread(case.fname, type="tsv") %>% normalize_df_field_names() %>%
	mutate(position = region_center(start, end));

d <- data.frame(
	g_score = c(s.control$g_score, s.case$g_score),
	group = c(rep("control", nrow(s.control)), rep("case", nrow(s.case)))
);

# compare distribution of G scores
ggplot(d, aes(x = g_score)) + geom_histogram(bins=100) + facet_grid(group ~ .) + theme_bw()


#chrom <- 8;
#cna.type <- "Amp";
#pos <- 48100000;

#chrom <- 11;
#cna.type <- "Amp";
#pos <- 55700000;

#chrom <- 4;
#cna.type <- "Amp";
#pos <- 52700000;

#chrom <- 14;
#cna.type <- "Amp";
#pos <- 19100000;

#chrom <- 5;
#cna.type <- "Amp";
#pos <- 46100000;

chrom <- 19;
cna.type <- "Amp";
pos <- 28600000;

s1 <- filter(s.control, chromosome == chrom, position >= pos, type == cna.type);
s2 <- filter(s.case, chromosome == chrom, position >= pos, type == cna.type);

#s1 <- filter(s.control, chromosome == chrom, position <= pos, type == cna.type);
#s2 <- filter(s.case, chromosome == chrom, position <= pos, type == cna.type);

# a very small proportion of overlapping positions between the two groups
length(intersect(s1$position, s2$position)) / length(union(s1$position, s2$position))

add_noise <- function(x, ...) {
	x + rnorm(length(x), ...)
}

# add noise to position to make them unique
#s1$position <- add_noise(s1$position);
#s2$position <- add_noise(s2$position);

#stopifnot(length(intersect(s1$position, s2$position)) == 0)

data <- list(
	J = nrow(s1) + nrow(s2),
	x = c(s1$position, s2$position) / 1e6,
	g = c(rep(-0.5, nrow(s1)), rep(0.5, nrow(s2))),
	y = c(s1$g_score, s2$g_score)
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
			lines(x[g <= 0], y[g <= 0], col="grey", type="b", pch=20);
			lines(x[g > 0], y[g > 0], col="orange", type="b", pch=20);
		}
	),
	width = 10,
	#file = "gscore_chr8q-amp_bm-luad_tcga-luad.pdf"
	#file = "gscore_chr11q-amp_bm-luad_tcga-luad.pdf"
	#file = "gscore_chr4q-amp_bm-luad_tcga-luad.pdf"
	#file = "gscore_chr14q-amp_bm-luad_tcga-luad.pdf"
	#file = "gscore_chr5p-amp_bm-luad_tcga-luad.pdf"
	file = "gscore_chr19q-amp_bm-luad_tcga-luad.pdf"
);

#qwrite(data, "gistic-chr8-amp_bm-luad_tcga-luad.rds");
#qwrite(data, "gistic-chr8q-amp_bm-luad_tcga-luad.rds");
#qwrite(data, "gistic-chr11q-amp_bm-luad_tcga-luad.rds");
#qwrite(data, "gistic-chr4q-amp_bm-luad_tcga-luad.rds");
#qwrite(data, "gistic-chr14q-amp_bm-luad_tcga-luad.rds");
#qwrite(data, "gistic-chr5p-amp_bm-luad_tcga-luad.rds");
qwrite(data, "gistic-chr19q-amp_bm-luad_tcga-luad.rds");

