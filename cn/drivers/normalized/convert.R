#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)
library(io, quietly=TRUE)

pr <- arg_parser("Convert copy-number in TSV format to RDS format with removal of name field");
pr <- add_argument(pr, "input", help="input TSV file");
pr <- add_argument(pr, "--outdir", default=".", help="output directory");

argv <- parse_args(pr);

in.fname <- as.filename(argv$input);
out.fname <- set_fext(set_fpath(in.fname, argv$outdir), "rds");

x <- qread(in.fname);

# delete name field
x$name <- NULL;
names(x) <- c("chromosome", "start", "end", "value");

table(x$chromosome)

chroms <- c(1:22, "X", "Y");
x$chromosome <- factor(x$chromosome, levels=chroms);

table(x$chromosome)

qwrite(x, out.fname);
