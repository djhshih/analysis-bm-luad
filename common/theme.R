theme_bar <- function() {
	theme_bw() +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(face="bold.italic"),
		#axis.text.x = element_text(angle = 90, hjust = 1),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		#panel.grid.major.x = element_blank(),
		#panel.grid.minor.x = element_blank(),
		panel.grid = element_blank(),
		legend.position = "bottom"
	)
}

theme_vdot <- function() {
	theme_bw() +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(face="bold.italic"),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		panel.grid = element_blank(),
		legend.position = "bottom"
	)
}

theme_dot <- function() {
	theme_bw() +
	theme(
		strip.background = element_blank(),
		strip.text = element_text(face="bold.italic"),
		panel.grid = element_blank(),
		legend.position = "none"
	)
}

theme_line <- function() {
	theme_bw() +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank()
	)
}

theme_box <- function() {
	theme_bw() +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank()
	)
}
