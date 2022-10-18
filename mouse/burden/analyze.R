library(io)
library(ggplot2)

x <- qread("tumor-burden_2019-12-10.csv");

out.fname <- filename("tumor-burden");


groups <- colnames(x);

d <- do.call(rbind, lapply(groups, function(g) {
	y <- x[[g]];
	idx <- !is.na(y);
	y <- y[idx];
	data.frame(value = y, group = g)
}));

d2 <- d;
d2$value[d2$value < 1] <- 1;

h <- kruskal.test(value ~ group, data = d2);

p_to_str <- function(x) {
	sprintf("p = %.2f", x)
}

p.str <- p_to_str(h$p.value);

ns <- table(d2$group);
d2$group_label <- factor(d2$group, levels=names(ns), labels=paste0(names(ns), "\n(", ns, ")"));

qdraw(
	ggplot(d2, aes(x=group_label, y=value)) +
		geom_boxplot(outlier.shape=NA) + theme_bw() +
		#geom_jitter(alpha=0.3, colour="darkorange", width=0.1) +
		scale_y_log10() +
		theme(legend.position = "none", panel.grid = element_blank()) +
		annotate(x = 0.75, y = 1e1, geom="text", label = p_to_str(h$p.value), hjust=0) +
		xlab("") + ylab("Tumor burden (photon flux)")
	,
	width = 3, height = 3,
	file = insert(out.fname, tag="box", ext="pdf")
)

