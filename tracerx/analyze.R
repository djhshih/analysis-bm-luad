library(io)
library(dplyr)

pheno <- qread("tracerx-sample-info.tsv");
cn <- qread("tracerx-cn.tsv");

amp.genes <- c(
	MYC = "chr8:128748315-128753680",
	YAP1 = "chr11:101981192-102104154",
	MMP13 = "chr11:102813721-102826463"
);

del.genes <- c(
	CDKN2A = "chr9:21967751-21994490",
	CDKN2B = "chr9:22002902-22009312"
);

string_to_coordinate <- function(x) {
	xs <- strsplit(sub("chr", "", x, fixed=TRUE), ":|-");
	y <- data.frame(do.call(rbind, xs), stringsAsFactors=FALSE);
	colnames(y) <- c("chromosome", "start", "end");
	y$chromosome <- as.integer(y$chromosome);
	y$start <- as.integer(y$start);
	y$end <- as.integer(y$end);
	y
}

#extreme <- function(x) {
#	x[which.max(abs(x))]
#}

collapse_sample_group <- function(x) {
	sub("-.+$", "", x)
}

groups.luad <- filter(pheno, grepl("adenocarcinoma", Histology))$TRACERxID;
print(length(groups.luad))

amp.coords <- string_to_coordinate(amp.genes);
del.coords <- string_to_coordinate(del.genes);

pad <- 1e6;

cn.f <- cn %>%
	# exclude lymph nodes
	filter(!grepl("-LN", sample)) %>%
	mutate(group = collapse_sample_group(sample)) %>%
	# use only selected sample groups
	filter(group %in% groups.luad);


print(collapse_sample_group(cn.f$sample) %>% unique() %>% length())

extreme <- max;
coords <- amp.coords;

ys <- lapply(1:nrow(coords),
	function(i) {
		x <- coords[i, ];

		f <- filter(cn.f,
				chr == x$chromosome,
				pmax(startpos, x$start - pad) <= pmin(endpos, x$end + pad)
			) %>%
			# collapse multiple regions
			select(group, cn = cnTotal) %>%
			#transmute(group, cn = nAraw + nBraw) %>%
			group_by(group) %>%
			summarize(cn = extreme(cn));

		f$cn[match(groups.luad, f$group)];
	}
);
names(ys) <- rownames(coords);

y <- data.frame(do.call(cbind, ys));
rownames(y) <- groups.luad;

y.amp <- y;


extreme <- min;
coords <- del.coords;

ys <- lapply(1:nrow(coords),
	function(i) {
		x <- coords[i, ];

		f <- filter(cn.f,
				chr == x$chromosome,
				pmax(startpos, x$start - pad) <= pmin(endpos, x$end + pad)
			) %>%
			# collapse multiple regions
			select(group, cn = cnTotal) %>%
			#transmute(group, cn = nAraw + nBraw) %>%
			group_by(group) %>%
			summarize(cn = extreme(cn));

		f$cn[match(groups.luad, f$group)];
	}
);
names(ys) <- rownames(coords);

y <- data.frame(do.call(cbind, ys));
rownames(y) <- groups.luad;

y.del <- y;

####

# thresholds chosen to reproduce counts called
# in the TracerX paper
amp.cut <- 10;
del.cut <- 0.5;

n <- length(groups.luad);

m.amp <- apply(y.amp, 2, function(z) sum(z > amp.cut));
m.amp
m.amp / n

h.amp <- lapply(m.amp, binom.test, n = n, conf.level=0.8);
print(h.amp)

m.del <- apply(y.del, 2, function(z) sum(z < del.cut));
m.del
m.del / n

h.del <- lapply(m.del, binom.test, n = n, conf.level=0.8);
print(h.del)


