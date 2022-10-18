library(io);
library(dplyr);

amp.cut <- 8;
amp.cut2 <- 6;
del.cut <- 0.5;
del.cut2 <- 0.6;

filter.chroms <- FALSE;
cluster.genes <- TRUE;

amp.drivers <- c("YAP1", "MMP13", "MYC");
del.drivers <- c("CDKN2A", "CDKN2B");

drivers <- c(amp.drivers, del.drivers); 
names(drivers) <- drivers;

compart.levels <- c("Brain metastasis", "Primary");
names(compart.levels) <- compart.levels;

cn <- as.matrix(qread("~/exigens/brain-mets/comut/brain-mets_pass_luad_absolute-1-4_gene_corrected_CN.txt", type="tsv"));

pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage3_pass_luad.tsv");
pheno <- pheno[pheno$sample_type %in% compart.levels, ];
#pheno <- pheno[pheno$pick, ];

out.fname <- filename("bm-comparts");
rds.fname <- insert(out.fname, ext="rds");

####

# construct compartments
patients <- as.character(unique(pheno$clinical_id));
names(patients) <- patients;

comparts <- lapply(patients,
	function(patient) {
		idx <- which(pheno$clinical_id == patient);
		if (length(idx) > 0) {
			samples <- as.character(pheno$sample_id[idx]);
			sample_types <- pheno$sample_type[idx];
			names(sample_types) <- samples;
			sample_types
		} else {
			NULL
		}
	}
);


idx.valid <- unlist(lapply(comparts,
	function(compart) {
		all(compart.levels %in% compart)
	}
));
patients.valid <- patients[idx.valid];
length(patients.valid)

samples.valid <- as.character(pheno$sample_id[pheno$clinical_id %in% patients.valid]);
cn.f <- cn[, colnames(cn) %in% samples.valid];


## Identify recurrently amplified or deleted genes

samples <- filter(pheno, clinical_id %in% patients.valid, sample_type %in% compart.levels)$sample_id %>%
	as.character();
samples <- intersect(samples, colnames(cn));

# do not consider CN outlier samples
cn.burden <- apply(cn.f, 2, function(z) sum(z > amp.cut | z < del.cut));
cn.burden.cut <- quantile(cn.burden, 0.95);
idx <- cn.burden >= cn.burden.cut;
outliers <- names(idx)[idx];
samples <- setdiff(samples, outliers);
cn.sel <- cn.f[, samples];

# select genes that are recurrently amplified: amplified in > 1 patient
patients.sel <- pheno$clinical_id[match(samples, pheno$sample_id)];
#amp.idx <- apply(cn.sel, 1, function(z) sum(z > amp.cut) > 1);
amp.npats<- apply(cn.sel, 1, function(z) {
	sum(tapply(z > amp.cut, patients.sel, any), na.rm=TRUE)
});
amp.idx <- amp.npats > 1;
amp.genes <- rownames(cn.sel)[amp.idx];

# select genes that are recurrently deleted: deleted in > 1 patient
#del.idx <- apply(cn.sel, 1, function(z) sum(z < del.cut) > 1);
del.npats<- apply(cn.sel, 1, function(z) {
	sum(tapply(z < del.cut, patients.sel, any), na.rm=TRUE)
});
del.idx <- del.npats > 1;
del.genes <- rownames(cn.sel)[del.idx];

fragile.site.genes <- c("WWOX", "FHIT", "CSMD1", "MACROD2", "PTPRD");

print(intersect(fragile.site.genes, del.genes))


if (filter.chroms) {

	# remove genes on same chromosome as driver genes, to avoid calling genes near
	# driver genes, which would be co-amplified or co-deleted

	library(org.Hs.eg.db);

	amp.genes.annot <- AnnotationDbi::select(org.Hs.eg.db, keys=amp.genes, keytype="SYMBOL", columns=c("CHR"));
	# YAP1 and MMP13 are on chr11
	# MYC is on chr8
	amp.genes.remove <- amp.genes.annot$SYMBOL[is.na(amp.genes.annot$CHR) | amp.genes.annot$CHR %in% c("8", "11")];
	amp.genes <- setdiff(amp.genes, amp.genes.remove);

	del.genes.annot <- AnnotationDbi::select(org.Hs.eg.db, keys=del.genes, keytype="SYMBOL", columns=c("CHR"));
	# CDKN2A is on chr9
	del.genes.remove <- del.genes.annot$SYMBOL[is.na(del.genes.annot$CHR) | del.genes.annot$CHR %in% c("9")];
	del.genes <- setdiff(del.genes, del.genes.remove);

	detach("package:org.Hs.eg.db", unload=TRUE);

}

if (cluster.genes) {

	collapse_clusters <- function(genes, x, blacklist) {
		d <- dist(x[genes, ], method="euclidean");
		hc <- hclust(d, method="single");
		clusters <- cutree(hc, h = 0.001);
		K <- max(clusters);
		clusters.l <- lapply(1:K, function(k) which(clusters == k));
		remove.idx <- which(unlist(lapply(clusters.l, function(z) any(blacklist %in% names(z)))));

		names(unlist(lapply(clusters.l[-remove.idx], function(x) x[1])));
	}

	amp.genes <- collapse_clusters(amp.genes, cn.f, amp.drivers);
	del.genes <- collapse_clusters(del.genes, cn.f, del.drivers);

} else {

	amp.genes <- setdiff(amp.genes, amp.drivers);
	del.genes <- setdiff(del.genes, del.drivers);

}

##

classes <- c(compart.levels, "Shared");
names(classes) <- classes;

ph <- pheno[match(colnames(cn.f), pheno$sample_id), ];

call_gene <- function(gene, direction) {
	if (direction > 0) {
		idx <- which(cn.f[gene, ] > amp.cut);
	} else {
		idx <- which(cn.f[gene, ] < del.cut);
	}

	pats <- unique(as.character(ph$clinical_id[idx]));
	names(pats) <- pats;

	# initial events with 0
	events <- lapply(pats, function(pat) {
		y <- rep(0, length(compart.levels));
		names(y) <- compart.levels;
		y
	});

	# call events
	for (i in idx) {
		events[[ as.character(ph$clinical_id[i]) ]][[ as.character(ph$sample_type[i]) ]] <- 1;
	}

	# rescue events based for patients with at least one called sample
	samples.sel <- ph$sample_id[ph$clinical_id %in% pats];
	ph2 <- ph[match(samples.sel, ph$sample_id), , drop=FALSE];
	cn.f2 <- cn.f[, colnames(cn.f) %in% samples.sel, drop=FALSE];
	# use secondary threshold to rescue events
	if (direction > 0) {
		idx2 <- which(cn.f2[gene, ] > amp.cut2);
	} else {
		idx2 <- which(cn.f2[gene, ] < del.cut2);
	}
	for (i in idx2) {
		pa <- as.character(ph2$clinical_id[i]);
		cp <- as.character(ph2$sample_type[i]);
		if (events[[pa]][[cp]] == 0) {
			events[[pa]][[cp]] <- 0.5;
		}
	}

	events
}

call_samples <- function(gene, direction) {
	if (direction > 0) {
		idx <- which(cn.f[gene, ] > amp.cut);
	} else {
		idx <- which(cn.f[gene, ] < del.cut);
	}

	calls <- numeric(ncol(cn.f));
	names(calls) <- colnames(cn.f);
	calls[idx] <- 1;


	pats <- unique(as.character(ph$clinical_id[idx]));
	names(pats) <- pats;

	# initial events with 0
	events <- lapply(pats, function(pat) {
		y <- rep(0, length(compart.levels));
		names(y) <- compart.levels;
		y
	});

	# call events
	for (i in idx) {
		events[[ as.character(ph$clinical_id[i]) ]][[ as.character(ph$sample_type[i]) ]] <- 1;
	}

	# rescue events based for patients with at least one called sample
	samples.sel <- ph$sample_id[ph$clinical_id %in% pats];
	ph2 <- ph[match(samples.sel, ph$sample_id), , drop=FALSE];
	cn.f2 <- cn.f[, colnames(cn.f) %in% samples.sel, drop=FALSE];
	# use secondary threshold to rescue events
	if (direction > 0) {
		idx2 <- which(cn.f2[gene, ] > amp.cut2);
	} else {
		idx2 <- which(cn.f2[gene, ] < del.cut2);
	}
	for (i in idx2) {
		pa <- as.character(ph2$clinical_id[i]);
		cp <- as.character(ph2$sample_type[i]);
		sa <- as.character(ph2$sample_id[i]);
		if (events[[pa]][[cp]] == 0) {
			events[[pa]][[cp]] <- 0.5;
			calls[sa] <- 0.5;
		}
	}

	calls
}

classify_gene <- function(gene, direction) {
	events <- call_gene(gene, direction);
	calls <- unlist(lapply(events, function(z) {
		if (sum(unlist(z) > 0) >= 2) {
			"Shared"
		} else {
			names(z)[which(z > 0)]
		}
	}));
	calls
}

count_classes <- function(genes, direction) {
	callss <- lapply(genes, classify_gene, direction = direction)
	counts <- lapply(
		callss,
		function(calls) {
			unlist(lapply(classes, function(z) sum(calls == z)));
		}
	);

	t(matrix(unlist(counts), nrow=length(classes), dimnames=list(classes, genes)));
}

####

call_gene("YAP1", direction = 1)
call_gene("MMP13", direction = 1)
call_gene("MYC", direction = 1)
call_gene("CDKN2B", direction = -1)
call_gene("CDKN2A", direction = -1)

####

call_patients <- function(gene, direction) {
	calls <- call_gene(gene, direction);
	if (direction > 0) {
		dir.type <- "amp";
	} else {
		dir.type <- "del";
	}
	d <- do.call(rbind, calls);
	colnames(d) <- paste(abbreviate(tolower(compart.levels), 2), tolower(gene), dir.type, sep="_");
	d <- data.frame(
		clinical_id = rownames(d),
		d
	);
	# append other samples
	pats <- as.character(unique(ph$clinical_id));
	other.pats <- setdiff(pats, d$clinical_id);
	other.d <- data.frame(clinical_id = other.pats);
	for (j in colnames(d)[-1]) {
		other.d[[j]] <- 0;
	}
	y <- rbind(d, other.d);
	rownames(y) <- NULL;

	# standardize patient order
	y <- y[match(pats, y$clinical_id), ];

	y
}

pat.amp.calls <- lapply(amp.drivers, call_patients, direction = 1);
pat.del.calls <- lapply(del.drivers, call_patients, direction = -1);

pat.calls <- c(pat.amp.calls, pat.del.calls);

pat.calls.merged <- do.call(cbind, pat.calls);
pat.calls.merged <- data.frame(
	clinical_id = pat.calls.merged[, 1],
	pat.calls.merged[, colnames(pat.calls.merged) != "clinical_id"]
);
summary(pat.calls.merged)

qwrite(pat.calls.merged, insert(out.fname, tag=c("calls"), ext="tsv"))

####

smp.amp.calls <- lapply(amp.drivers, call_samples, direction = 1);
names(smp.amp.calls) <- amp.drivers;

smp.del.calls <- lapply(del.drivers, call_samples, direction = -1);
names(smp.del.calls) <- del.drivers;

smp.calls.merged <- do.call(cbind, c(smp.amp.calls, smp.del.calls));
smp.calls.merged <- data.frame(sample_id = rownames(smp.calls.merged), smp.calls.merged);
rownames(smp.calls.merged) <- NULL;

# filter out samples that are neither primary nor brain metastasis
samples.sel <- as.character(ph$sample_id[ph$sample_type %in% c("Primary", "Brain metastasis")]);
length(unique(ph$clinical_id[ph$sample_id %in% samples.sel]));

qwrite(smp.calls.merged, insert(out.fname, tag=c("calls", "samples"), ext="tsv"));

####

counts.amp.drivers.m <- count_classes(amp.drivers, direction = 1);
print(counts.amp.drivers.m)

counts.del.drivers.m <- count_classes(del.drivers, direction = -1);
print(counts.del.drivers.m)

counts.amp.m <- count_classes(amp.genes, direction = 1);
marginals.amp <- apply(counts.amp.m, 2, sum);
print(marginals.amp)

counts.del.m <- count_classes(del.genes, direction = -1);
marginals.del <- apply(counts.del.m, 2, sum);
print(marginals.del)

qwrite(counts.amp.drivers.m, insert(rds.fname, "amp-drivers"));
qwrite(counts.del.drivers.m, insert(rds.fname, "del-drivers"));
qwrite(counts.amp.m, insert(rds.fname, "amp-genes"));
qwrite(counts.del.m, insert(rds.fname, "del-genes"));

qwrite(amp.genes, insert(out.fname, tag="amp-genes", ext="vtr"));
qwrite(del.genes, insert(out.fname, tag="del-genes", ext="vtr"));


####

# Amplifications

fisher.test(matrix(c(4 + 2 + 4, 12300, 1 + 2 + 0, 7694), nrow=2))
fisher.test(matrix(c(4, 12300, 0, 7694), nrow=2))
fisher.test(matrix(c(4 + 2 + 4, 12300, 1 + 2 + 0, 7694), nrow=2), alternative="greater")
fisher.test(matrix(c(4, 12300, 0, 7694), nrow=2), alternative="greater")

fisher.test(matrix(c(4 + 2 + 4, 12300, 3 + 3 + 5, 12098), nrow=2))
fisher.test(matrix(c(1 + 2 + 0, 12300, 3 + 3 + 5, 12098), nrow=2))

fisher.test(matrix(c(4 + 2 + 4, 5164, 1 + 2 + 0, 3189), nrow=2))
fisher.test(matrix(c(4 + 2 + 4, 5164, 3 + 3 + 5, 5318), nrow=2))
fisher.test(matrix(c(1 + 2 + 0, 5164, 3 + 3 + 5, 5318), nrow=2))

fisher.test(matrix(c(4 + 2 + 4, 2177, 1 + 2 + 0, 1449), nrow=2))
fisher.test(matrix(c(4 + 2 + 4, 2177, 3 + 3 + 5, 10045), nrow=2))
fisher.test(matrix(c(1 + 2 + 0, 1449, 3 + 3 + 5, 10045), nrow=2))

# B P S
# 2177 1449 10045
fisher.test(matrix(c(2 + 4 + 2, 2177, 1 + 1 + 0, 1449), nrow=2))
fisher.test(matrix(c(2 + 4 + 2, 2177, 3 + 3 + 6, 10045), nrow=2))
fisher.test(matrix(c(1 + 1 + 0, 1449, 3 + 3 + 6, 10045), nrow=2))

# B P S
# 5387 3086 23619
fisher.test(matrix(c(2 + 4 + 2, 5387, 1 + 1 + 0, 3086), nrow=2))
fisher.test(matrix(c(2 + 4 + 2, 5387, 3 + 3 + 6, 23619), nrow=2))
fisher.test(matrix(c(1 + 1 + 0, 3086, 3 + 3 + 6, 23619), nrow=2))

# B P S
# 1829 896 8548
fisher.test(matrix(c(2 + 4 + 2, 1829, 1 + 1 + 0, 896), nrow=2))
fisher.test(matrix(c(2 + 4 + 2, 1829, 3 + 3 + 6, 8548), nrow=2))
fisher.test(matrix(c(1 + 1 + 0, 896, 3 + 3 + 6, 8548), nrow=2))

# B P S
# 2399 1127 8542
fisher.test(matrix(c(2 + 4 + 3, 2399, 1 + 2 + 0, 1127), nrow=2))
fisher.test(matrix(c(2 + 4 + 3, 2399, 3 + 3 + 6, 8542), nrow=2))
fisher.test(matrix(c(1 + 2 + 0, 1127, 3 + 3 + 6, 8542), nrow=2))

# Deletions

fisher.test(matrix(c(3, 394, 1, 1301), nrow=2))
fisher.test(matrix(c(3, 394, 13, 956), nrow=2))
fisher.test(matrix(c(1, 1301, 13, 956), nrow=2))

fisher.test(matrix(c(3, 132, 1, 626), nrow=2))
fisher.test(matrix(c(3, 132, 13, 125), nrow=2))
fisher.test(matrix(c(1, 626, 13, 125), nrow=2))

# B P S
# 93 310 123
fisher.test(matrix(c(3, 93, 1, 310), nrow=2))
fisher.test(matrix(c(3, 93, 13, 123), nrow=2))
fisher.test(matrix(c(1, 310, 13, 123), nrow=2))

# B P S
# 326 572 611
fisher.test(matrix(c(3, 326, 1, 572), nrow=2))
fisher.test(matrix(c(3, 326, 13, 611), nrow=2))
fisher.test(matrix(c(1, 572, 13, 611), nrow=2))

# B P S
# 154 222 306
fisher.test(matrix(c(3, 154, 1, 222), nrow=2))
fisher.test(matrix(c(3, 154, 13, 306), nrow=2))
fisher.test(matrix(c(1, 222, 13, 306), nrow=2))

# B P S
# 358 313 390
fisher.test(matrix(c(3, 358, 1, 313), nrow=2))
fisher.test(matrix(c(3, 358, 13, 390), nrow=2))
fisher.test(matrix(c(1, 313, 13, 390), nrow=2))

