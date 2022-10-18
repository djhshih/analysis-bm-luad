library(io);
library(dplyr);

x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");
pheno <- qread("../annot/sample-info_wes_stage2.tsv");

out.fname <- filename("luad-bm_cor");

###

mapping <- match(x$pheno$tumor_id, pheno$fh_sample_id);

# update sample names with standardized annotation
colnames(x$expr) <- pheno$sample_id[mapping];

x$pheno <- mutate(x$pheno,
	sample_id = pheno$sample_id[mapping],
	sample_type = pheno$sample_type[mapping],
	primary_histotype = pheno$primary_histotype[mapping]
);

###

sample.types <- c("Primary", "Brain metastasis");

# select only LUAD
luad.idx <- x$pheno$primary_histotype == "Lung adenocarcinoma" & x$pheno$sample_type %in% sample.types;
x$pheno <- x$pheno[luad.idx, ];
x$cn <- x$cn[, luad.idx];
fit <- lm(z[j, ] ~ z[i, ] + purity + group + purity:group + z[i,]:group);
x$expr <- x$expr[, luad.idx];
x$mut <- x$mut[luad.idx, ];
x$amp <- x$amp[luad.idx, ];
x$gsva <- lapply(x$gsva, function(z) z[, luad.idx]);

str(x)

iqrs <- apply(x$expr, 1, IQR);

table(x$pheno$sample_type)
x$pheno$sample_type <- factor(x$pheno$sample_type, levels = sample.types);


gene.idx <- iqrs > 0.5;
mean(gene.idx)

#z <- x$expr[gene.idx, ];
z <- x$expr[order(iqrs, decreasing=TRUE), ];
z <- z[1:100, ];
purity <- x$pheno$purity;

ref.group <- "Primary";
group <- x$pheno$sample_type != ref.group;

z1 <- x$expr[gene.idx, group];
z2 <- x$expr[gene.idx,!group];
purity1 <- x$pheno$purity[group];
purity2 <- x$pheno$purity[!group];

#lm(z1[j, ] ~ z1[i, ] + purity1)
#lm(z2[j, ] ~ z2[i, ] + purity2)

# number of parameters
d <- 6;

i <- 1;
j <- 2;

fit <- lm(z[j, ] ~ z[i, ] + purity + group + z[i,]:group);
coef(summary(fit))[p, 3]

# t test on coefficient is done on n - p degrees of freedom
coef

G <- nrow(z);
s <- matrix(0, nrow=G, ncol=G, dimnames=list(x=rownames(z), y=rownames(z)));
p <- matrix(1, nrow=G, ncol=G, dimnames=list(x=rownames(z), y=rownames(z)));

for (i in 1:G) {
	if (i %% 10 == 0) {
		message("i = ", i);
	}
	for (j in 1:G) {
		if (i != j) {
			fit <- lm(z[j, ] ~ z[i, ] + purity + group + purity:group + z[i,]:group);
			# record the t value of the interaction term (difference in slope of
			# group 2 vs. group 1
			coefs <- coef(summary(fit));
			s[i, j] <- coefs[d, 3];
			p[i, j] <- coefs[d, 4];
		}
	}
}

i <- 48;
j <- 82;
fit <- lm(z[j, ] ~ z[i, ] + purity + purity:group + group + z[i,]:group);
summary(fit);

fit <- lm(z[i, ] ~ z[j, ] + purity + purity:group + group + z[j,]:group);
summary(fit);

rownames(z)[c(i, j)]
# PLP1 CNS proteolipid
# SFTPA2 lung surfactant

plot(group, z[j, ])
dev.off()

plot(purity[!group], z[j, !group])
dev.off()

plot(z[i, !group], z[j, !group])
dev.off()

plot(z[j, !group], z[i, !group])
dev.off()



# @param m  number of rows
# @param n  number of columns
# lower triangular matrix not including the diagonal
ind2sub_ltri <- function(k, m, n) {
	# number of rows in each column
	ms <- seq.int(m-1, by=-1, length.out=n);
	# cumulative number of elements including all elements of column i 
	# in column-major order
	cms <- cumsum(ms);
	j <- which(k < cms)[1];
	i <- k - cms[j-1] + j;
	c(i, j)
}

idx <- lower.tri(s);
s.min <- pmin(abs(s[idx]), abs(t(s)[idx]));
max(s.min)
min.idx <- ind2sub_ltri(which.max(s.min), nrow(s), ncol(s));


df.res <- ncol(z) - p;
alpha <- 0.05;
t.crit <- qt(1 - alpha/2, df.res);
summary(s.min)
mean(s.min > t.crit)


which(s.min > 3)
min.idx <- ind2sub_ltri(2296, nrow(s), ncol(s));


i <- min.idx[1];
j <- min.idx[2];
s[i, j]
fit <- lm(z[j, ] ~ z[i, ] + purity + group + purity:group + z[i,]:group);
summary(fit)

i <- min.idx[2];
j <- min.idx[1];
s[i, j]
fit <- lm(z[j, ] ~ z[i, ] + purity + group + purity:group + z[i,]:group);
summary(fit)

i <- min.idx[1];
j <- min.idx[2];
fit <- lm(z[j, ] ~ z[i, ] + purity + group + purity:group + z[i,]:group);
summary(fit);
rownames(z)[c(i, j)]

library(io);
library(ggplot2);
d <- data.frame(
	x = z[i, ],
	y = z[j, ],
	group = x$pheno$sample_type,
	purity = x$pheno$purity
);
qdraw(
	ggplot(d, aes(x=x, y=y, colour=purity, shape=group)) + geom_point() + theme_bw(),
	file = "tmp.pdf",
	width = 6
);

