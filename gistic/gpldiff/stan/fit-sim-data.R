library(rstan);

data <- readRDS("gp-compare_sim-data.rds");

# do not use Nystr\"om approximation
data$M <- data$J;

options(mc.cores=4);
fit <- stan(
	file="gp-compare.stan",
	data=data, iter=500, chains=4
);

f.mcmc <- extract(fit, "f")[[1]];
f.mcmc.mean <- apply(f.mcmc, 2, mean);

cor(data$f, f.mcmc.mean)

