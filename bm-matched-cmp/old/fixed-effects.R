library(io);
library(ggplot2);
library(reshape2);
library(tidyr);
library(dplyr);

source("mid-p-binom.R");

# Analyze matched BM and primary with candidate CNA drivers
# using the mid p binomial test (e.g. a fixed-effect model)

# 73 brain metastases
# 57 matched primary

out.fname <- filename("pkb-luad", tag="bm-matched-cmp");
pdf.fname <- insert(out.fname, ext="pdf");

compart <- qread("cna-compart-specificity.tsv");
genes <- as.character(unique(compart$gene));

timing <- qread("cna-bm-timing.tsv");

# added shared compartment
shared <- filter(timing, epoch == "early") %>% mutate(compartment = "Truncal") %>% select(gene, compartment, count);

branch <- rbind(compart, shared) %>% rename(branch=compartment)
branch$branch <- factor(branch$branch, levels=c("Primary", "Brain metastasis", "Truncal"))

####

exclusivity.p <- mid_p_binomial_test(
	sum(compart$count[compart$compartment == "Brain metastasis"]),
	sum(compart$count)
);

compart <- mutate(compart,
	exclusivity = ifelse(compartment == "Brain metastasis", "BM-enriched", "Primary-enriched")
);

qdraw(
	ggplot(compart, aes(x=exclusivity, y=count, fill=gene)) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		theme(
			legend.text = element_text(face = "bold.italic"),
			axis.text.x = element_text(angle=30, hjust=1)
		) +
		xlab("") + ylab("# patients") +
		ylim(0, 15) +
		annotate(x = 1.5, y = 14, geom="text", label = sprintf("p = %s", format(exclusivity.p, digits=3)))
	,
	width = 2.5, height=3,
	file = insert(pdf.fname, c("bm-exclusivity", "bar"))
);

qdraw(
	ggplot(branch, aes(x=branch, y=count, fill=gene)) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		theme(
			legend.text = element_text(face = "bold.italic"),
			axis.text.x = element_text(angle=30, hjust=1)
		) +
		xlab("Phylogenetic branch") + ylab("# patients") +
		annotate(x = 1.5, y = 14, geom="text", label = sprintf("p = %s", format(exclusivity.p, digits=3)))
	,
	width = 3, height=4,
	file = insert(pdf.fname, c("branch", "bar"))
);

qdraw(
	ggplot(branch, aes(x=factor(branch, levels=rev(levels(branch))), y=count, fill=gene)) +
		coord_flip() +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		labs(fill="") +
		theme(
			legend.text = element_text(face = "bold.italic")
		) +
		xlab("") + ylab("# patients") +
		annotate(x = 2.5, y = 18, geom="text", label = sprintf("p = %s", format(exclusivity.p, digits=3)))
	,
	width = 5, height = 1.5,
	file = insert(pdf.fname, c("branch", "hbar"))
);

####

early.p <- mid_p_binomial_test(
	sum(timing$count[timing$epoch == "early"]),
	sum(timing$count)
);

qdraw(
	ggplot(timing, aes(x=epoch, y=count, fill=gene)) +
		geom_bar(stat="identity") + theme_bw() +
		scale_fill_brewer(palette="Set2") +
		theme(
			legend.text = element_text(face = "bold.italic"),
			axis.text.x = element_text(angle=30, hjust=1)
		) +
		xlab("") + ylab("# patients") +
		ylim(0, 26) +
		annotate(x = 1.5, y = 25.5, geom="text", label = sprintf("p = %s", format(early.p, digits=2)))
	,
	width = 2.5, height=3,
	file = insert(pdf.fname, c("bm-timing", "bar"))
);

