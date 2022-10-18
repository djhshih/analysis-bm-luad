library(io);
library(gpldiff);
library(ggplot2);
library(RColorBrewer);

load_all("~/projects/r/gpldiff");

gscores.fname <- "../luad_absolute_bmet-only_cnvf/scores.gistic";
gscores <- read_gistic(gscores.fname);

gscores.control.fname <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_absolute_cnvf_wt/scores.gistic";
gscores.control <- read_gistic(gscores.control.fname);

out.fname <- filename("pkb-tcga-luad", tag="gistic");
pdf.fname <- insert(out.fname, ext="pdf");

split_gscores <- function(gscores) {
	lapply(
		split(gscores, list(type=gscores$type)),
		function(d) split(d, list(chromosome=d$chromosome))
	);
}


insert_gap <- function(track, gap.cut = 0.1) {
	x.diff <- diff(track$x);
	domain <- track$x[length(track$x)] - track$x[1];
	
	if (max(x.diff) / domain > gap.cut) {
		# insert zero at the largest gap
		i <- which.max(x.diff);
		track$x <- c(track$x[1:i], track$x[i], track$x[i+1], track$x[(i+1):length(track$x)])
		track$y <- c(track$y[1:i], 0, 0, track$y[(i+1):length(track$y)])
	}

	track
}

make_chromosome_track <- function(d) {
	pos <- c(d$start, d$end);
	scores <- c(d$g_score, d$g_score);
	idx <- order(pos);
	track <- list(x = pos[idx], y = scores[idx]);

	data.frame(chromosome = d$chromosome[1], insert_gap(track))
}

int_to_chrom <- function(x, species="hsa") {
	if (species == "hsa") {
		factor(x, levels=1:24, labels=c(1:22, "X", "Y"))
	} else {
		x
	}
}

make_chromosome_tracks <- function(d) {
	tracks <- do.call(rbind, lapply(d, make_chromosome_track));
	tracks$chromosome <- int_to_chrom(tracks$chromosome);
	rownames(tracks) <- NULL;
	tracks
}

#' @param xs a named list of data.frame
combine_cohorts <- function(xs) {
	cohorts <- names(xs);
	# TODO extend this to multiple cohorts
	d <- rbind(
		data.frame(
			cohort = cohorts[1],
			xs[[1]]
		),
		data.frame(
			cohort = cohorts[2],
			xs[[2]]
		)
	);
	d$cohort <- factor(d$cohort, levels=cohorts);
	d
}


cohorts <- c("TCGA", "Present");

gscores.split <- split_gscores(gscores);
amp.tracks <- make_chromosome_tracks(gscores.split$Amp);

gscores.control.split <- split_gscores(gscores.control);
amp.tracks.control <- make_chromosome_tracks(gscores.control.split$Amp);

amp.cols <- brewer.pal(3, "Reds")[2:3];

amp.tracks.combined <- combine_cohorts(list(TCGA = amp.tracks.control, Present = amp.tracks));


qdraw(
	{
		ggplot(amp.tracks.combined, aes(x=x, y=y, colour=cohort)) + geom_line() + theme_bw() + 
			facet_grid(. ~ chromosome, scales="free_x", space="free_x") +
			xlab("Genomic position") + ylab("GISTIC score") +
			scale_colour_manual(values=amp.cols) +
			theme(
				panel.spacing = unit(0.2, "lines"),
				strip.background = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank()
			)
	},
	width = 20,
	file = insert(pdf.fname, c("genome", "amp"))
);


gscores.split <- split_gscores(gscores);
del.tracks <- make_chromosome_tracks(gscores.split$Del);

gscores.control.split <- split_gscores(gscores.control);
del.tracks.control <- make_chromosome_tracks(gscores.control.split$Del);

del.cols <- brewer.pal(3, "Blues")[2:3];

del.tracks.combined <- combine_cohorts(list(TCGA = del.tracks.control, Present = del.tracks));

qdraw(
	{
		ggplot(del.tracks.combined, aes(x=x, y=y, colour=cohort)) + geom_line() + theme_bw() + 
			facet_grid(. ~ chromosome, scales="free_x", space="free_x") +
			xlab("Genomic position") + ylab("GISTIC score") +
			scale_colour_manual(values=del.cols) +
			theme(
				panel.spacing = unit(0.2, "lines"),
				strip.background = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank()
			)
	},
	width = 20,
	file = insert(pdf.fname, c("genome", "del"))
);

