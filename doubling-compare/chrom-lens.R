library(io)

cytoband <- read.table("cytoBand.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE);
colnames(cytoband) <- c("chromosome", "start", "end", "name", "stain");
cytoband$chromosome <- sub("chr", "", cytoband$chromosome);

cytoband.nocen <- cytoband[ cytoband$stain != "acen", ];
cytoband.nocen$length <- cytoband.nocen$end - cytoband.nocen$start;

chrom.lens <- tapply(cytoband.nocen$length, cytoband.nocen$chromosome, sum);

chroms <- names(chrom.lens);
chroms[chroms == "X"] <- 23;
chroms[chroms == "Y"] <- 24;
chroms <- as.integer(chroms);

chrom.lens <- chrom.lens[order(chroms)];

qwrite(chrom.lens, "chrom-lens_no-cen.vtr");

