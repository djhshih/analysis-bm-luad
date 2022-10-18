mid_p_binomial_test <- function(x, n, prob=0.5, alternative=c("two.sided", "less", "greater"), mid.p=TRUE) {
	alternative <- match.arg(alternative);

	# alternative
	# less: x < n - x
	# greater: x > n - x
	
	if (alternative == "less") {
		# integrate over i \forall i <= x
		exact.p <- pbinom(x, n, prob);
	} else if (alternative == "greater") {
		# integrate over i \forall i >= x
		exact.p <- pbinom(x - 1, n, prob, lower.tail = FALSE);
	} else if (alternative == "two.sided") {
		exact.p <- 2 * min(
			pbinom(x, n, prob),
			pbinom(x - 1, n, prob, lower.tail = FALSE)
		);
	}
	
	if (mid.p) {
		# calculate the mid p-value by subtracting half of the 
		# probability mass at x
		max(0, exact.p - 0.5 * dbinom(x, n, prob))
	} else {
		exact.p
	}
}

