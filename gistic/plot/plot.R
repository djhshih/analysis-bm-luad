library(io);
library(gpldiff);
library(ggplot2);
library(RColorBrewer);
library(dplyr);

load_all("~/projects/r/gpldiff");

gscores.fname <- "../luad_absolute_bmet-only_cnvf/scores.gistic";
gscores <- read_gistic(gscores.fname);

gscores.control.fname <- "~/share/exigens/tcga/tcga-luad/gistic/absolute/tcga-luad_absolute_cnvf_wt/scores.gistic";
gscores.control <- read_gistic(gscores.control.fname);

out.fname <- filename("pkb-tcga-luad", tag="gistic");
pdf.fname <- insert(out.fname, ext="pdf");

source("~/exigens/brain-mets/common/params.R");

#amp.cols <- c("grey30", brewer.pal(3, "Reds")[3]);
#del.cols <- c("grey30", brewer.pal(3, "Blues")[3]);

amp.cols <- c("grey60", brewer.pal(3, "Reds")[3]);
del.cols <- c("grey60", brewer.pal(3, "Blues")[3]);

amp.fill.cols <- c("green4", brewer.pal(3, "Reds")[3]);
del.fill.cols <- c("redorange2", brewer.pal(3, "Blues")[3]);

split_gscores <- function(gscores) {
	lapply(
		split(gscores, list(type=gscores$type)),
		function(d) split(d, list(chromosome=d$chromosome))
	);
}


insert_gap <- function(track, gap.cut = 0) {
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

gscores.split <- split_gscores(gscores);
amp.tracks <- make_chromosome_tracks(gscores.split$Amp);

gscores.control.split <- split_gscores(gscores.control);
amp.tracks.control <- make_chromosome_tracks(gscores.control.split$Amp);


xs <- list(amp.tracks.control, amp.tracks);
names(xs) <- cohorts;
amp.tracks.combined <- combine_cohorts(xs);

gistic_plot <- function(tracks) {
	ggplot(tracks, aes(x=x, y=y, colour=cohort)) + 
		geom_line(alpha=0.7) + theme_bw() + 
		facet_grid(. ~ chromosome, scales="free_x", space="free_x") +
		xlab("Genomic position") + ylab("GISTIC score") +
		theme(
			panel.spacing = unit(0.2, "lines"),
			strip.background = element_blank(),
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.border = element_rect(colour="grey80"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		);
}


gistic_fill_plot <- function(tracks) {
	ggplot(tracks, aes(x=x, ymin=0, ymax=y, fill=cohort)) + 
		geom_ribbon(alpha=0.6) + theme_bw() + 
		facet_grid(. ~ chromosome, scales="free_x", space="free_x") +
		xlab("Genomic position") + ylab("GISTIC score") +
		theme(
			panel.spacing = unit(0.2, "lines"),
			strip.background = element_blank(),
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.border = element_rect(colour="grey80"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		);
}

gistic_region_plot <- function(tracks) {
	tracks$x <- tracks$x / 1e6;
	ggplot(tracks, aes(x=x, y=y, colour=cohort)) + 
		geom_line(alpha=0.7) + theme_bw() + 
		xlab("Genomic position (Mb)") + ylab("GISTIC score") +
		ylim(min(tracks$y), max(tracks$y) * 1.25) +
		theme(
			panel.border = element_rect(colour="grey80"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			legend.position = "bottom"
		);
}

gistic_fill_region_plot <- function(tracks) {
	tracks$x <- tracks$x / 1e6;
	ggplot(tracks, aes(x=x, ymin=0, ymax=y, fill=cohort)) + 
		geom_ribbon(alpha=0.6) + theme_bw() + 
		xlab("Genomic position (Mb)") + ylab("GISTIC score") +
		theme(
			panel.border = element_rect(colour="grey80"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			legend.position = "bottom"
		);
}


g.amp <- gistic_plot(amp.tracks.combined) + scale_colour_manual(values=amp.cols);

qdraw(
	g.amp,
	width = 20,
	file = insert(pdf.fname, c("genome", "amp"))
);

qdraw(
	g.amp,
	width = 10,
	height = 2,
	file = insert(pdf.fname, c("genome", "amp", "small"))
);


g.fill.amp <- gistic_fill_plot(amp.tracks.combined) + scale_fill_manual(values=amp.cols);

qdraw(
	g.fill.amp,
	width = 20,
	file = insert(pdf.fname, c("genome", "fill", "amp"))
);

qdraw(
	g.fill.amp,
	width = 10,
	height = 2,
	file = insert(pdf.fname, c("genome", "fill", "amp", "small"))
);


gscores.split <- split_gscores(gscores);
del.tracks <- make_chromosome_tracks(gscores.split$Del);

gscores.control.split <- split_gscores(gscores.control);
del.tracks.control <- make_chromosome_tracks(gscores.control.split$Del);


xs <- list(del.tracks.control, del.tracks);
names(xs) <- cohorts;
del.tracks.combined <- combine_cohorts(xs);

g.del <- gistic_plot(del.tracks.combined) + scale_colour_manual(values=del.cols);

qdraw(
	g.del,
	width = 20,
	file = insert(pdf.fname, c("genome", "del"))
);

qdraw(
	g.del,
	width = 10,
	height = 2,
	file = insert(pdf.fname, c("genome", "del", "small"))
);


g.fill.del <- gistic_fill_plot(del.tracks.combined) + scale_fill_manual(values=del.cols);

qdraw(
	g.fill.del,
	width = 20,
	file = insert(pdf.fname, c("genome", "fill", "del"))
);

qdraw(
	g.fill.del,
	width = 10,
	height = 2,
	file = insert(pdf.fname, c("genome", "fill", "del", "small"))
);

####

# chr8:126445706-130762261
r <- filter(amp.tracks.combined, chromosome == 8, x >= 100e6, x <= 150e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("myc", "amp"), ext="pdf")
);

qdraw(
	gistic_fill_plot(r) + scale_fill_manual(values=amp.fill.cols) + theme(legend.position="bottom"),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("fill", "myc", "amp"), ext="pdf")
);

# chr11:101863542-102585362
r <- filter(amp.tracks.combined, chromosome == 11, x >= 80e6, x <= 120e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("yap1", "amp"), ext="pdf")
);

# chr4:148893121-148899092
r <- filter(amp.tracks.combined, chromosome == 4, x >= 130e6, x <= 170e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("ednra", "amp"), ext="pdf")
);

# chr14:36334948-39510012
r <- filter(amp.tracks.combined, chromosome == 14, x >= 10e6, x <= 80e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("nkx2-1", "amp"), ext="pdf")
);

# chr19:38189869-40883989
r <- filter(amp.tracks.combined, chromosome == 19, x >= 35e6, x <= 60e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("ccne1", "amp"), ext="pdf")
);

# chr1:149897699-153391702 (1q21)
r <- filter(amp.tracks.combined, chromosome == 1, x >= 140e6, x <= 170e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("1q21", "amp"), ext="pdf")
);

# chr12:66935663-69645890
r <- filter(amp.tracks.combined, chromosome == 12, x >= 45e6, x <= 90e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=amp.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("mdm2", "amp"), ext="pdf")
);

# chr9:21227888-23693464
r <- filter(del.tracks.combined, chromosome == 9, x >= 0e6, x <= 40e6);

qdraw(
	gistic_region_plot(r) + scale_colour_manual(values=del.cols),
	width = 3, height = 3,
	file = insert(out.fname, tag=c("cdkn2a-cdkn2b", "del"), ext="pdf")
);

