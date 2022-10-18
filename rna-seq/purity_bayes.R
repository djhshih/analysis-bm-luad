library(io);
library(RColorBrewer);
library(magrittr);
library(limma);
library(dplyr);
library(ggplot2);
library(MASS);
library(reshape2);
library(cluster);

library(rstan);
library(tntrna);
load_all("~/projects/r/tntrna/");

x <- qread("luad-bm_cn-expr-snv-pheno-gsva.rds");

met.idx <- x$pheno$met_type == "BM";
y <- x$expr["GFAP", met.idx];
y <- x$expr["NKX2-1", met.idx];
y <- x$expr["AGO2", met.idx];
y <- x$expr["MMP13", met.idx];
y <- x$expr["MMP2", met.idx];
y <- x$expr["MMP1", met.idx];
y <- x$expr["ACTB", met.idx];
y <- x$expr["GAPDH", met.idx];
purity <- x$pheno$purity[met.idx];

pri.idx <- x$pheno$met_type == "Primary";
y <- x$expr["NKX2-1", pri.idx];
y <- x$expr["MMP13", pri.idx];
y <- x$expr["GFAP", pri.idx];
purity <- x$pheno$purity[pri.idx];

y <- x$expr["NKX2-1", ];
y <- x$expr["MMP13", ];
y <- x$expr["MMP2", ];
y <- x$expr["MMP1", ];
y <- x$expr["MYC", ];
y <- x$expr["YAP1", ];
purity <- x$pheno$purity;


lfit <- lm(y ~ purity);
lmbeta <- c(coef(lfit)[1], coef(lfit)[2]);

ffit <- rlm(y ~ purity);
fbeta <- c(coef(ffit)[1], coef(ffit)[2]);
summary(ffit);

data <- list(
	N = length(y),
  y = y,
	p = purity,
	nu = 1,
	mu = 5,
	tau = 2,
	kappa = 2,
	omega = 0.3,
	psi = 1
);

bfit <- tntrna_sample(data, "norm-exp", iter=2e4);

# expexp tends to give a lot of local minima in p(theta0, theta1)

theta0s <- extract(bfit, "theta0")[[1]];
theta1s <- extract(bfit, "theta1")[[1]];
sigmas <- extract(bfit, "sigma")[[1]];
thetas <- c(median(theta0s), median(theta1s));
bbeta <- c(thetas[1], thetas[2] - thetas[1]);

hist(theta0s, breaks=100);
hist(theta1s, breaks=100);
smoothScatter(theta0s, theta1s);
abline(a=0, b=1, col="gray30");
hist(c(theta0s, theta1s), breaks=100);

hist(sigmas, breaks=100);

theta.mat <- matrix(c(theta0s, theta1s), ncol=2);
idx <- 1:500;
k <- pam(theta.mat[idx, ], 3);
plot(theta0s[idx], theta1s[idx], pch=".", col=k$cluster);

plot(purity, y, xlim=c(0, 1));
abline(a=lmbeta[1], b=lmbeta[2], col="orange");
abline(a=fbeta[1], b=fbeta[2], col="firebrick");
abline(a=bbeta[1], b=bbeta[2], col="royalblue");


# Visual distributions

n <- 1e4;

hist(rgamma(n, 1, 5), breaks=100);
hist(rgamma(n, 5, 1), breaks=100);
hist(c(rgamma(n*0.3, 1, 5), rgamma(n*0.7, 5, 1)), breaks=100);
hist(c(rexp(n*0.3, 100), rgamma(n*0.7, 5, 1)), breaks=100);

z <- c(rexp(n*0.3, 100), rnorm(n*0.7, 5, 2));
hist(z[z >= 0], breaks=100);

z <- c(rexp(n*0.3, 10), rnorm(n*0.7, 5, 2));
hist(z[z >= 0], breaks=100);

rb0 <- c(rgamma(n*0.3, 1, 5), rgamma(n*0.7, 5, 1));
rb1 <- c(rgamma(n*0.3, 1, 5), rgamma(n*0.7, 5, 1));

smoothScatter(rb0, rb1);
smoothScatter(rb0, rb1 - rb0);

plot(rb0, rb1, pch=".");
plot(rb0, rb1 - rb0, pch=".");
abline(a=0, b=-1, col="royalblue");

hist(rexp(n, 0.5), breaks=100);


load_all("~/projects/r/tntrna/");
params0.init <- list(theta=3, sigma=1);
params1.init <- list(theta=c(3, 3), sigma=1);


# prior of theta pulls it aggressive to 0: need to avoid initializing at 0
gemfit0 <- tntrna_gem_normexp(data, init=params0.init);
gemfit <- tntrna_gem_normexp(data, init=params1.init, max.iter=2e3);
nfit <- laplace(tntrna_joint_distrib, gemfit$parameters, data=data)
nfit0 <- laplace(tntrna_joint_distrib, gemfit0$parameters, data=data)

tntrna_hessian(gemfit$parameters, data);
nfit$parameters$Sigma


tntrna_joint_distrib(gemfit$parameters, data);
tntrna_joint_distrib(gemfit0$parameters, data);

gemfit1 <- tntrna_gem_normexp(data, init=list(theta=c(0, 3), sigma=3), max.iter=2e3);
tntrna_joint_distrib(gemfit1$parameters, data);
gemfit2 <- tntrna_gem_normexp(data, init=list(theta=c(3, 0), sigma=3), max.iter=2e3);
tntrna_joint_distrib(gemfit2$parameters, data);
gemfit3 <- tntrna_gem_normexp(data, init=list(theta=c(0, 0), sigma=3), max.iter=2e3);
tntrna_joint_distrib(gemfit3$parameters, data);


tntrna_gradient_plot(gemfit$parameters, data);

tntrna_joint_distrib(list(theta=c(-1, 0), sigma=1), data);

theta_joint_distrib <- function(sigma) {
	function(theta, data) {
		tntrna_joint_distrib(list(theta=theta, sigma=sigma), data)
	}
}

f <- theta_joint_distrib(gemfit$parameters$sigma);
f(gemfit$parameters$theta, data)
f(c(-1, 0), data)

d <- seq(0, 10, by=0.1);
l <- unlist(lapply(d, function(x) f(c(x, gemfit$parameters$theta[2]), data)));
r <- unlist(lapply(d, function(x) f(c(gemfit$parameters$theta[1], x), data)));
dd <- seq(0, 10, 0.01);
s <- unlist(lapply(dd,
	function(x) tntrna_joint_distrib(list(theta = gemfit$parameters$theta, sigma=x), data)));
plot(d, l, type="l");
plot(d, r, type="l");
plot(dd, s, type="l");

mycontour(f, c(0, 8, 0, 8), data);


mycontour(lbinorm, c(0, 8, 0, 8), list(m=nfit$mode, v=nfit$var));


model_posterior(log(c(0.5, 0.5)), c(nfit0$lp.data, nfit$lp.data));
compare_bayes_models(list(nfit0, nfit), lods=TRUE);

plot(purity, y, xlim=c(0, 1));
abline(a=lmbeta[1], b=lmbeta[2], col="orange");
abline(a=fbeta[1], b=fbeta[2], col="firebrick");
abline(a=bbeta[1], b=bbeta[2], col="royalblue");
abline(a=gemfit$parameters$theta[1], b=diff(gemfit$parameters$theta), col="green4");
abline(a=gemfit0$parameters$theta, b=0, col="yellowgreen");


nfit.mean <- nfit$parameters$mu;
nfit.sd <- sqrt(diag(nfit$parameters$Sigma));

old.par <- par(mfrow=c(3,1));

hist(theta0s, breaks=100, freq=FALSE);
curve(dnorm(x, nfit.mean[1], nfit.sd[1]), add=TRUE, col="firebrick", n=201);

hist(theta1s, breaks=100, freq=FALSE);
curve(dnorm(x, nfit.mean[2], nfit.sd[2]), add=TRUE, col="firebrick", n=201);

hist(sigmas, breaks=100, freq=FALSE);
curve(dnorm(x, nfit.mean[3], nfit.sd[3]), add=TRUE, col="firebrick", n=201);

par(old.par);


gemfit0 <- tntrna_gem_normexp(data, init=list(theta=3, sigma=1));
gemfit <- tntrna_gem_normexp(data, init=list(theta=c(3, 3), sigma=1));
nfit0 <- laplace(tntrna_joint_distrib, gemfit0$parameters, data=data);
nfit <- laplace(tntrna_joint_distrib, gemfit$parameters, data=data);

plot(purity, y, xlim=c(0, 1));
abline(a=gemfit$parameters$theta[1], b=diff(gemfit$parameters$theta), col="green4");
abline(a=gemfit0$parameters$theta, b=0, col="yellowgreen");

