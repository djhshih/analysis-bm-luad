
N <- 100;
sigma <- 0.1;

x <- seq(0, 4*pi, length.out=N) + rnorm(N, sd = 0.1);

y_a <- rnorm(N, mean = sin(x), sd = sigma);
y_b <- rnorm(N, mean = sin(x) + pi/3, sd = sigma);

g <- sample(c(-0.5, 0.5), N, replace=TRUE);

f <- y_b - y_a;

data <- list(
	# inputs
	J = N,
	x = x,
	y = ifelse(g > 0, y_b, y_a),
	g = g,
	# answers
	f = f,
	y_a = y_a,
	y_b = y_b,
	sigma = sigma
);

saveRDS(data, "gp-compare_sim-data.rds");

