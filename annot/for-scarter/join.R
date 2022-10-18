library(io);
library(dplyr);

x <- qread("../sample-info_wes_stage2_pass_luad.tsv");
y <- qread("../patient-info.rds");

z <- left_join(x, y, by="clinical_id");

qwrite(z, "sample-info_wes2_stage2_pass_luad_w-patient-info.tsv");

