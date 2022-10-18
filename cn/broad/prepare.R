# ==============================================================================
# PURPOSE
# To determine broad events from segmentation profiles
#
# @Author:   David JH Shih  (djh.shih@gmail.com)
# @License:  GNU General Public License v3 
# @Created:  2011-11-30
# @Input:    segmentation file
# @Output:   matrix of broad events

# ==============================================================================
# HISTORY
#
# Version:   0.1
# Date:      2011-11-30
# Comment:   Initial write

# ==============================================================================
# PREAMBLE
#
library(bioinf);
library(io);
library(filenamer);
library(dplyr);

cytoband <- qread("~/data/ucsc/hg19/cytoBand.rds");

input.fname <- "~/exigens/brain-mets/gistic/brain-mets_pass_luad_absolute_bmet-only.lgr.cntf.cenf.cnvf.seg";
output.fstem <- "pkb-bm-luad_brain-mets";

#input.fname <- "~/exigens/brain-mets/gistic/tcga-luad_pass_absolute.lgr.cntf.cenf.cnvf.seg";
#output.fstem <- "tcga-luad";


#state.cut <- 0.2;
#state.cut <- 0.15;
state.cut <- 0.1;

remove.chry <- TRUE;

# ==============================================================================
# FUNCTIONS
# 

get.bounded.region <- function(region, bound) {
	if (region$chromosome != bound$chromosome) return (NULL);
	s <- max(region$start, bound$start);
	e <- min(region$end, bound$end);
	if (e < s) return (NULL);
	region$start <- s;
	region$end <- e;
	return ( region );
}

get.frequency <- function(x, state) apply(x == state, 1, sum) / ncol(x);

# ==============================================================================
# INPUT
# 

chromosomes.levels <- unique(cytoband$cytobands$chromosome);
arms.levels <- unique(cytoband$chrom.arms$arm);

chromosomes <- rep(chromosomes.levels, each=length(arms.levels));
arms <- rep(arms.levels, times=length(chromosomes.levels));

seg <- qread(input.fname);

# make seg play nice
colnames(seg) <- c("sample", "chromosome", "start", "end", "count", "state");
seg$chromosome <- sub("chr", "", seg$chromosome);


# ==============================================================================
# PROCESS
# 

x <- cytoband$chrom.arms;
x$chromosome <- as.character(x$chromosome);

chrom.arms <- as.list( as.data.frame( mapply( function(chrom, arm) {
	return ( x[x$chromosome == chrom & x$arm == arm, c("chromosome", "start", "end")] );
}, chromosomes, arms) ) );

names(chrom.arms) <- paste(chromosomes, arms, sep="");

if (remove.chry) {
	# remove Y chromosome and arms
	chrom.arms <- chrom.arms[ !( names(chrom.arms) %in% c("Y", "Yp", "Yq") ) ];
}


seg.samp <- split(seg, list(seg$sample));

mat <- matrix(nrow=length(chrom.arms), ncol=length(seg.samp));
rownames(mat) <- names(chrom.arms);
colnames(mat) <- names(seg.samp);


for (chrom in rownames(mat)) {
	message("Processing ", chrom);
	for (samp in colnames(mat)) {
		x <- seg.samp[[samp]][ seg.samp[[samp]]$chromosome == chrom.arms[[chrom]]$chromosome, ];
		if (nrow(x) > 0) {
			regions <- lapply(df.to.row.list(x), function(r) {
					return (get.bounded.region(r, chrom.arms[[chrom]]));
			});
			regions.df <- row.list.to.df(regions);
			total.length <- sum(regions.df$end - regions.df$start + 1);
			# sum copy number state from segments, weighted by size of segment
			state <- sum( unlist( lapply(regions, function(r) {
					(r$end - r$start + 1)/total.length * r$state
			}) ) );

			mat[chrom, samp] <- state;
		}
	}
}

mat.cna <- matrix(0, nrow=nrow(mat), ncol=ncol(mat));
rownames(mat.cna) <- rownames(mat);
colnames(mat.cna) <- colnames(mat);
mat.cna[mat > state.cut] <- 1;
mat.cna[mat < -state.cut] <- -1;

freq.gains <- get.frequency(mat.cna, 1);
freq.losses <- get.frequency(mat.cna, -1);


# ==============================================================================
# OUTPUT
# 

qwrite(mat, filename(output.fstem, tag="state", ext="rda"));
qwrite(mat.cna, filename(output.fstem, tag="cna", ext="mtx"));

