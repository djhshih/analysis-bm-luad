library(io);
library(dplyr);

x <- qread("An_ALL_Agilent_BrainMets_individuals.tsv");

dup <- duplicated(x$collaborator_id) | duplicated(x$collaborator_id, fromLast=TRUE);

y <- x[dup, ] %>%
	select(individual_id, collaborator_id);

y2 <- y[order(y$collaborator_id), ];


qwrite(y2, "agilent_duplicated-individuals.tsv");

# cut -f 1 agilent_duplicated-individuals.tsv | sed '1,3d' > agilent_duplicated-individuals.txt
