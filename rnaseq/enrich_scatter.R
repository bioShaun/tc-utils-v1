suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggthemes))
suppressMessages(library(argparser))
suppressMessages(library(stringr))


options(stringsAsFactors = F)
p <- arg_parser("enrichment scatter plot")
p <- add_argument(p, "--enrich_file", help = "enrichment table file.")
argv <- parse_args(p)



wrap_long_name <- function(name, width=80) {
  return(paste(strwrap(name, width = width), collapse="\n"))
}

clean_enrich_table <- function(enrich_file) {
  enrich_file_name <- basename(enrich_file)
  enrich_df <- read.csv(enrich_file)
  enrich_df <- enrich_df[, c('p.adjust', 'Description', 'GeneRatio')]
  enrich_df$term <- enrich_df$Description
  enrich_df$numInCat <- as.numeric(str_extract(enrich_df$GeneRatio, "\\d+"))
  ylab_title <- '-log10(p.adjust)'
  
  enrich_plot_data_list <- list(table=enrich_df, title=ylab_title)
  return(enrich_plot_data_list)
}

save_ggplot <- function(ggplot_out, output,
                        width=8, height=6) {
  ggsave(paste(output, "png", sep = "."),
         plot = ggplot_out,
         width = width,
         height = height,
         dpi = 300, type = "cairo")
  ggsave(paste(output, "pdf", sep = "."),
         plot = ggplot_out,
         width = width,
         height = height,
         device = cairo_pdf)
}


cell_cols <- RColorBrewer::brewer.pal(10,"RdYlGn")
gradiant_cell_cols <- rev(colorRampPalette(cell_cols)(100))


enrich_file <- argv$enrich_file
out_prefix <- str_replace(enrich_file, 'enrichment.csv', 'enrichment.scatter_plot')
enrich_data_list <- clean_enrich_table(enrich_file)
enrich.table <- enrich_data_list$table
enrich.table$log_pval <- -log10(enrich.table$p.adjust)


enrich.table$wrap_term <- sapply(enrich.table$term, wrap_long_name)
enrich.table$wrap_term <- factor(enrich.table$wrap_term,
                              levels = rev(enrich.table$wrap_term))


plot_df <- head(enrich.table, 30)
max_term_len <- max(nchar(as.character(plot_df$term)))
max_term_len <- ifelse(max_term_len > 80, 80, max_term_len)

plot_witdh <- 6 + max_term_len *0.05
plot_height <- dim(plot_df)[1] / 3

p <- ggplot(plot_df, aes(log_pval, wrap_term, size=numInCat, color=log_pval)) +
  geom_point() +
  scale_color_gradientn(colours = gradiant_cell_cols)+
  theme_calc() +
  guides(size=guide_legend(title = 'Number of Genes'),
         color=guide_colorbar(title='-log10(Adjusted Pvalue)')) +
  xlab('-log10(Adjusted Pvalue)') + ylab('') +
  theme(axis.text.x = element_text(size = rel(1.2)),
        axis.title.x = element_text(size = rel(1),
                                    face = 'bold'),
        axis.text.y = element_text(size = rel(1.2)))
p

save_ggplot(p, out_prefix, width = plot_witdh, height = plot_height)

