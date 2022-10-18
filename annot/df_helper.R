library(magrittr);

normalize_field_names <- function(x) {
	tolower(x) %>%
	# remove everything after and including a colon
	#sub(":.*$", "", .) %>%
	# remove all parenthesized contents
	gsub(" ?\\([^)]*\\)", "", .) %>%
	gsub(" ?\\[[^]]*\\]", "", .) %>%
	# replace spaces and dashes with underscores
	gsub("[ -]", "_", .) %>%
	# replace repeated underscores
	gsub("__+", "_", .) %>%
	gsub("%", "pct", .) %>%
	abbreviate_terms %>%
	prioritize_word(., "date")
}

is_na_field <- function(x) {
	all(is.na(x)) && is.logical(x)
}

convert_df_na_fields_to_type <- function(d, type="character") {
	for (j in 1:ncol(d)) {
		if (is_na_field(d[, j])) {
			d[, j] <- as(d[, j], type);
		}
	}
	d
}

# field names should already be normalized
abbreviate_terms <- function(x) {
	x %>%
	gsub("diagnosis", "dx", .) %>%
	gsub("history", "hx", .) %>%
	gsub("surgery", "sx", .) %>%
	gsub("biopsy", "bx", .) %>%
	gsub("resection", "sx", .) %>%
	gsub("resxn", "sx", .) %>%
	gsub("treatment", "tx", .) %>%
	gsub("therapy", "tx", .) %>%
	gsub("chemotherapy", "ct", .) %>%
	gsub("chemo", "ct", .) %>%
	gsub("radiotherapy", "rt", .) %>%
	gsub("radiation", "rt", .) %>%
	gsub("primary_tumor", "primary", .) %>%
	gsub("metastasis", "met", .) %>%
	gsub("metastases", "met", .) %>%
	gsub("metastatic", "met", .) %>%
	gsub("progression_free_survival", "pfs", .) %>%
	gsub("progression", "pd", .) %>%
	gsub("brain_met", "bm", .)
}

# move word to the front if it appears at the end
prioritize_word <- function(x, word) {
	gsub(sprintf("(.+)_(%ss?)$", word), "\\2_of_\\1", x)
}

normalize_empty<- function(x) {
	x[x == ""] <- NA
	x
}

normalize_df_field_names <- function(d) {
	stopifnot(is.data.frame(d) || is.matrix(d));
	colnames(d) <- normalize_field_names(colnames(d));
	d
}

rename_df_field <- function(d, from, to) {
	i <- match(from, colnames(d));
	if (length(i) == 1) {
		colnames(d)[i] <- to;
	} else if (length(i) == 0) {
		stop("Field name cannot be found");
	} else {
		stop("Field name is not unique");
	}
	d
}

remove_df_field <- function(d, name) {
	i <- match(name, colnames(d));
	if (length(i) > 0) {
		d[, -i]
	} else {
		stop("field name ", name, " is not found");
	}
}


censor_survival <- function(x, max.time) {
  idx <- x$os_year > max.time;
  x$os_years[idx] <- max.time;
  x$death[idx] <- 0;
  x
}

remove_rare_levels <- function(x, n) {
  d <- table(x);
  factor(x, levels=names(d)[d >= n])
}

mutate_df_field <- function(d, field, f) {
	d[[field]] %>% f	
	d
}

capitalize <- function(x) {
	paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}

normalize_date <- function(x) {
	x <- as.character(x) %>%
		# strip anything in parantheses
		gsub("\\([^)]+\\)", "", .) %>%
		# strip non-legal characters (allow some characters for delimited values)
    gsub("[^0-9 ,/-]", "", .) %>%
		# strip leading and trailing whitespaces
		trimws;

	patterns <- c(
		"%Y-%m-%d" = "\\d{4}-\\d{1,2}-\\d{1,2}",
		"%m/%d/%Y" = "\\d{1,2}/\\d{1,2}/\\d{4}",
		"%m/%d/%y" = "\\d{1,2}/\\d{1,2}/\\d{2}"
	);

	pattern.idx <- rep(0, length(x));
	for (i in 1:length(patterns)) {
		idx <- pattern.idx == 0;
		if (sum(idx) > 0) {
			pattern.idx[idx][ grepl(patterns[i], x[idx]) ] <- i;
		}
	}
	# default to first pattern (0 indices causes hell occurs later)
	pattern.idx[pattern.idx == 0] <- 1;
		
  as.Date(x, names(patterns)[pattern.idx])
}

# TODO
normalize_dates <- function(x) {
	x <- as.character(x) %>%
		# strip anything in parantheses
		gsub("\\([^)]+\\)", "", .) %>%
		# strip non-legal characters (allow some characters for delimited values)
    gsub("[^0-9 ,/-]", "", .)

	# NB this assumes that dates are separate by - and ,
	#    so it will not work for POSIX format, e.g. 2016-01-01
	xs <- strsplit(x, "-|,");

	patterns <- c(
		"%Y-%m-%d" = "\\d{4}-\\d{1,2}-\\d{1,2}",
		"%m/%d/%Y" = "\\d{1,2}/\\d{1,2}/\\d{4}",
		"%m/%d/%y" = "\\d{1,2}/\\d{1,2}/\\d{2}"
	);

	# "surgery: 12/30/10, RT: 4/4/11-5/17/11"
  # "RT: 12/2010-1/21/11"                  
  # "RT 4/4/05-5/23/05"                    
	# "RT: 3/29/10-5/4/10"    
}

#' An R function for filling in missing values of a variable from one data frame with the values from another variable.
#'
#' \code{FillIn} uses values of a variable from one data set to fill in missing values in another.
#' 
#' @param x the data frame with the variable you would like to fill in.
#' @param y the data frame with the variable you would like to use to fill in \code{x}.
#' @param target a character string of the name of the variable in \code{x} you want to fill in.
#' @param src an optional character string of variable name in \code{y} that you would like to use to fill in.
#' @param by a character vector of variable names that are shared by \code{x} and \code{y} that can be used to join the data frames.
#' 
#' @examples 
#' # Create data set with missing values
#' naDF <- data.frame(a = sample(c(1,2), 100, rep=TRUE), 
#'                    b = sample(c(3,4), 100, rep=TRUE), 
#'                    fNA = sample(c(100, 200, 300, 400, NA), 100, rep=TRUE))
#'
#' # Created full data set
#' fillDF <- data.frame(a = c(1,2,1,2), 
#'                      b = c(3,3,4,4),
#'                      fFull = c(100, 200, 300, 400))
#'
#' # Fill in missing f's from naDF with values from fillDF
#' FilledInData <- FillIn(naDF, fillDF, target = "fNA", src = "fFull", by = c("a", "b"))
#'
fill_in <- function(x, y, target, src = NULL, by = intersect(colnames(x), colnames(y))) {
  # Give src the same name as var1 if src is NULL
  if (is.null(src)){
    src <- target
  } else {
    src <- src
  }
  
  # Give var a generic name
  names(x)[match(target, names(x))] <- "VarGen"
  names(y)[match(src, names(y))] <- "VarGen.1"
  
  # Convert data frames to data.table type objects
  xTemp <- data.table::data.table(x, key = by)
  yTemp <- data.table::data.table(y, key = by)
  
  # Merge data.tables
  OutDT <- yTemp[xTemp]
  
  # Tell the user how many values will be filled in
  SubNA <- OutDT[, list(VarGen, VarGen.1)]
  SubNA <- subset(SubNA, is.na(VarGen) & !is.na(VarGen.1))
  print(paste(nrow(SubNA), "NAs were replaced."))
  
  # Fill in missing values from x with values from y
  OutDT <- OutDT[is.na(VarGen), VarGen := VarGen.1]
  
  # Convert back to data frame
  OutDF <- data.frame(OutDT)
  
  # Tell the user what the correlation coefficient is between the variables
  SubNoNA <- subset(OutDF, !is.na(VarGen) & !is.na(VarGen.1))
  HowMany <- nrow(SubNoNA)
  CORR <- cor(SubNoNA$VarGen, SubNoNA$VarGen.1, use = "complete.obs")
  print(paste("The correlation between", target, "and", src, "is", round(CORR, digits = 3), "based on", HowMany, "shared observations." ))
  
  # Remove uncombined variable and return main variable's name
  names(OutDF)[match("VarGen", names(OutDF))] <- target
  Keepers <- setdiff(names(OutDF), "VarGen.1")
  OutDF <- OutDF[, Keepers]
  OutDF
}
