library(io);
library(dplyr);
library(ggplot2);
library(reshape2);

prepare <- function(d) {
	d$position <- ceiling((d$start + d$end) / 2);
	d
}

combine_cn <- function(xs) {
	do.call(rbind, mapply(
		function(x, id) {
			if (is.null(x)) {
				NULL
			} else {
				data.frame(x, id = id)
			}
		},
		xs,
		names(xs),
		SIMPLIFY=FALSE
	))
}

mark_region <- function(x, region.start, region.end) {
	mutate(x, marked = position >= region.start & position <= region.end)
}

region_plot <- function(cn.l, patient, chrom, window.lim, region.lim, direction=1) {
	window.start <- window.lim[1];
	window.end <- window.lim[2];
	region.start <- region.lim[1];
	region.end <- region.lim[2];

	patient.d <- combine_cn(cn.l) %>% 
		filter(chromosome == chrom, position >= window.start, position <= window.end) %>%
		mark_region(region.start, region.end);

	g <- ggplot(patient.d, aes(x = position / 1e6, y = value, colour=marked)) +
		geom_point(alpha=0.3) + theme_bw() + 
		facet_wrap(~ id, ncol=1, strip.position="top") +
		xlab("Genomic position (Mbp)") +
		ggtitle(sprintf("Patient %s on chr%s near", patient, chrom)) +
		ggtitle(substitute("Patient " * {p} * " on chr " * {chr} * " near " * {italic(g)},
			list(p = patient, chr = chrom, g = gene))) +
		theme(
			legend.position = "none",
			strip.background = element_blank(),
			strip.text = element_text(hjust = 0)
		)

	if (direction > 0) {
		g  <- g + scale_colour_manual(values=c("black", "red"));
	} else {
		g <- g + scale_colour_manual(values=c("black", "royalblue"));
	}

	g
}


# project specific function
sort_samples <- function(samples) {
	name_to_type <- function(name) {
		ifelse(grepl("-N", name), 0,
			ifelse(grepl("-P", name), 1,
				ifelse(grepl("-M", name), 3, 2)))
	}
	samples[ order(name_to_type(samples)) ]
}

#' @param r copy-ratios
#' @param ploidy  tumour ploidy
#' @param purity  tumour purity
rescaled_cn <- function(r, ploidy, purity) {
	ploidy * r  +  2*(1 - purity)/purity * (r - 1)
}

#' @param s  copy-number state
#' @param ploidy  tumour ploidy
#' @param purity  tumour purity
copy_ratio <- function(s, ploidy, purity) {
	v <- 2 * (1 - purity);
	(s * purity + v) / (ploidy * purity + v)
}


pheno <- qread("~/exigens/brain-mets/annot/sample-info_wes_stage2.tsv");

in.dir <- "rds";
out.fname <- filename("cnplot", path="png", ext="png", date=NA);
#out.fname <- filename("cnplot", path="png", ext="pdf", date=NA);

params <- list(
	YAP1 = list(
		chromosome = "11",
		window = c(90e6, 110e6),
		region = c(101981192, 102826463),
		direction = 1
	),
	MYC = list(
		chromosome = "8",
		window = c(120e6, 140e6),
		region = c(128748315, 128753680),
		direction = 1
	),
	CDKN2A = list(
		chromosome = "9",
		window = c(10e6, 30e6),
		region = c(21967751, 22009312),
		direction = -1
	)
);

#cna.samples <- qread("cna-samples.rds");
sample.sets <- qread("cna-patients.rds");

for (gene in names(sample.sets)) {
	message(gene)

	set <- sample.sets[[gene]];
	p <- params[[gene]];


	for (patient in names(set)) {
		message(patient)

		if (length(set[[patient]]) == 0) next;

		samples <- sort_samples(set[[patient]]);
		
		cn.l <- lapply(samples, function(s) {
			fname <- 	filename(s, ext=c("tn", "rds"), path=in.dir, date=NA);
			x <- qread(fname) %>% prepare();

			idx <- which(pheno$sample_id == s);
			if (pheno$sample_type[idx] == "Normal") {
				ploidy <- 2;
				purity <- 1;
			} else {
				ploidy <- pheno$tau[idx];
				purity <- pheno$purity[idx];
			}
			if (is.na(ploidy) || is.na(purity)) {
				warning("Sample ", s, " has no ploidy or purity estimate");
				x <- NULL;
			} else {
				s <- rescaled_cn(2^x$value, ploidy, purity);
				x$value <- s;	
			}

			x
		});
		names(cn.l) <- samples;
		
		g <- region_plot(cn.l, patient, p$chromosome, p$window, p$region, p$direction) + 
			geom_hline(yintercept=2, linetype="dotted") +
			ylab("Rescaled copy-number");

		if (p$direction < 0) {
			g <- g + ylim(-0.75, 4);
		}

		qdraw(
			g,
			width = 6, height = 10,
			device = NULL,
			file = insert(out.fname, c(patient, tolower(gene)))
		);
	}
}

