library(io);
library(rstan);

data <- qread("gistic-chr8q-amp_bm-luad_tcga-luad.rds");

# if y ~ normal(\mu, \sigma), then
# var(diff(y)) should have 2 * \sigma^2 variance;
# therefore,
sigma <- sqrt(var(diff(data$y))/2);
data$sigma <- sigma * 2;
# multiple by 2 since sigma is the std of the difference
# assume homoscesdascity of the two groups

# convert g from {0, 1} to {-0.5, 0.5}
data$g <- data$g - 0.5;

options(mc.cores=4);
fit <- stan(file="gp-compare_model3.stan", data=data, iter=500, chains=4);

out.fname <- filename("gp-compare", tag="model3");
qwrite(fit, insert(out.fname, ext="rds"));

