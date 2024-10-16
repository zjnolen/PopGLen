sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(data.table)
library(reshape2)
library(ggplot2)
library(ggtext)
library(cowplot)
library(dplyr)

plot_admix <- function(kvalues, qoptlist, pop, optsumm) {

	admix <- list()

	# Generate a plot for each K
	for (i in c(1:length(kvalues))) {
		k <- kvalues[i]
		qopt <- read.table(qoptlist[i])
		qopt$sample <- pop$sample
		qopt <- reshape2::melt(as.data.table(inner_join(qopt, pop)))
		qopt$time <- factor(
			qopt$time,
			levels = c("historical", "modern"),
			labels = c("Historical", "Modern")
		)
		qopt$population <- as.factor(qopt$population)

		# Show strip only on first plot
		if (i == 1) {
			striptext <- element_text(size = 8)
		} else {
			striptext <- element_text(size = 0)
		}

		# Show sample names only on last
		if (i == length(kvalues)) {
			xtext <- element_text(angle = 90)
		} else {
			xtext <- element_blank()
		}

		# Bold K values of converged Ks
		if (optsumm[optsumm$k == k, "converged"] == "Y") {
			ytext <- element_markdown(lineheight = 1.2, face = "bold")
		} else {
			ytext <- element_markdown(lineheight = 1.2)
		}

		# Make plot for given K
		admix[[i]] <- ggplot(qopt, aes(value, fill = variable, x = sample)) +
			geom_bar(
				position = "fill",
				stat = "identity",
				width = 1,
				color = "black"
			) +
			facet_grid(
				~ population,
				scales = "free",
				space = "free"
			) +
			ylab(paste0("K = ", k)) +
			guides(fill = "none") +
			scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.25)) +
			theme_classic() +
			theme(
				axis.title.y = ytext,
        strip.text.x = striptext,
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = xtext,
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0, 0, 0.1, 0), "cm")
			)
	}

	# Arrange plots in a grid with some rough evenness hopefully. Will probably be
	# a bit off due to variation in sample name length.
	plot_grid(
		plotlist = admix,
		nrow = length(admix),
		ncol = 1,
		rel_heights = c(1, rep(0.88, length(admix)-2), 1.25)
	)
}

# Read in necessary population values
pop <- as.data.frame(
	fread(
		snakemake@input[["poplist"]],
		header = TRUE,
		sep = "\t",
		select = c("sample", "population", "time")
	)
)


# Summarize convergence output to table, use it for plotting
optwraps <- c()

for (i in c(1:length(snakemake@params[["kvals"]]))) {
	k <- snakemake@params[["kvals"]][i]
	optwrap <- read.table(snakemake@input[["optwrap"]][i], header = FALSE)
	optwrap$k <- k
	optwraps <- rbind(optwraps, optwrap)
}

optsumm <- optwraps %>%
	group_by(k) %>%
	summarize(
		niter = max(V1),
		deltaLike = diff(
			range(
				sort(
					V3,
					decreasing = TRUE
				)[1:snakemake@params[["conv"]]]
			)
		)
	)

optsumm$converged <- ifelse(
	optsumm$deltaLike <= snakemake@params[["thresh"]],
	"Y",
	"N"
)

write.table(optsumm, snakemake@output[["convsumm"]], quote = FALSE, sep = "\t",
	row.names = FALSE)

# Run the plotting function
plot_admix(
	snakemake@params[["kvals"]],
	snakemake@input[["qopts"]],
	pop,
	optsumm
)

# Set reasonable dimensions based on number of pops and Ks
if (nrow(pop) <= 50) {
	plotwidth = 6
} else {
	plotwidth = (6/45)*length()
}

plotheight = 1.65 + (1.44 * (length(snakemake@params[["kvals"]])-2)) + 2

# Save plot
ggsave(
	snakemake@output[[1]],
	width = plotwidth,
	height = plotheight,
	units = "in"
)