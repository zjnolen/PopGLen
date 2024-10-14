sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(reactable)
library(data.table)
library(htmlwidgets)

dt <- fread(snakemake@input[[1]], header = TRUE)
table <- reactable(dt, 
	searchable = TRUE, 
	filterable = TRUE, 
	showPageSizeOptions = TRUE,
	defaultPageSize = 25,
	paginationType = "jump",
	resizable = TRUE,
	wrap = FALSE)

saveWidget(widget = table, file = snakemake@output[[1]], selfcontained = TRUE)