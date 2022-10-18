library(reshape2);
library(dplyr);
library(ggplot2);
library(io);

# calculate median from binned data
median_binned <- function(values, counts) {
	n <- sum(counts);
	m <- ceiling(n/2);
	sums <- cumsum(counts);
	values[which(sums >= m)[1]]
}

pdf.fname <- filename("insert-qc", ext="pdf");

qc <- qread("../luad-bm_qc.rds");

files <- list.files(pattern = "\\.insert_size_metrics$");
names(files) <- gsub(".insert_size_metrics", "", files, fixed=TRUE);

hs <- lapply(files, function(f) {
	d <- qread(f, type="tsv", skip=11)[, 1:2];
	colnames(d) <- c("insert_size", "count");
	d
});

summary.files <- list.files(pattern = "\\.insert_size_metrics.head$");
names(summary.files) <- names(files);
ss <- lapply(summary.files, function(f) {
	qread(f, type="tsv")[1,]
});

reported_insert_medians <- unlist(lapply(ss, function(x) x$MEDIAN_INSERT_SIZE));

hs.m <- melt(hs, id.vars="insert_size") %>% select(sample = L1, insert_size=insert_size, count=value);

qc$prop_zero_pass <- qc$prop_zero < 0.1;

hs.m.qc <- left_join(hs.m, select(qc, sample=rnaseq_id, prop_zero_pass=prop_zero_pass)) %>%
	filter(sample %in% qc$rnaseq_id[!qc$fresh_frozen]);

g <- ggplot(hs.m.qc, aes(x=insert_size, y=count, fill=sample)) + geom_bar(stat="identity") +
	facet_grid(prop_zero_pass ~ ., scale="free_y");
qdraw(g, insert(pdf.fname, "hist"), width=8);

g <- ggplot(hs.m.qc, aes(x=insert_size, y=count, fill=sample)) + 
	geom_bar(stat="identity", position="dodge") +
	xlim(0, 200) + facet_grid(prop_zero_pass ~ ., scale="free_y");
qdraw(g, insert(pdf.fname, c("hist", "dodge")), width=8);

g <- ggplot(hs.m.qc, aes(x=insert_size, y=count, fill=sample)) + 
	geom_density(stat="identity", alpha=0.2) +
	xlim(0, 200) + facet_grid(prop_zero_pass ~ ., scale="free_y");
qdraw(g, insert(pdf.fname, c("density")), width=8);

insert_medians <- unlist(lapply(hs,
	function(h) {
		median_binned(h$insert_size, h$count)
	}
));

insert_counts <- unlist(lapply(hs,
	function(h) sum(h$count)
));

stopifnot(names(insert_medians) == names(reported_insert_medians));
insert.df <- data.frame(
	rnaseq_id=names(insert_medians),
	insert_median=insert_medians,
	insert_count=insert_counts,
	reported_insert_median=reported_insert_medians
);

insert.qc <- left_join(qc, insert.df) %>% filter(!fresh_frozen);

ggplot(insert.qc, aes(x=insert_count, y=prop_zero)) + geom_point();
ggplot(insert.qc, aes(x=insert_count, y=insert_median)) + geom_point();
ggplot(insert.qc, aes(x=reported_insert_median, y=insert_median)) + geom_point();

ggplot(insert.qc, aes(x=insert_median, y=prop_zero)) + geom_point();

g <- ggplot(insert.qc, aes(x=reported_insert_median, y=prop_zero, colour=dv200)) +
	geom_point() + scale_x_continuous(breaks=seq(100, 600, by=50));
qdraw(g, insert(pdf.fname, "rmedian-vs-propzero"));

g <- ggplot(insert.qc, aes(x=reported_insert_median, y=dv200, colour=prop_zero > 0.1)) + 
	geom_point();
qdraw(g, insert(pdf.fname, "rmedian-vs-dv200"));

