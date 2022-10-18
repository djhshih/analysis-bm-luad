library(io);
library(dplyr);

x <- qread("../agilent_duplicated-individuals.tsv") %>% 
	filter(collaborator_id != "");

individuals <- as.character(unique(x$collaborator_id));

for (individual in individuals) {
	individual_ids <- as.character(x$individual_id[x$collaborator_id == individual]);
	samples <- lapply(individual_ids,
		function(individual_id) {
			qread(sprintf("samples_%s.tsv", individual_id))
		}
	);
}

