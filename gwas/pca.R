suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(argparser))

options(stringsAsFactors = F)

p <- arg_parser("vcf pca")
p <- add_argument(
  p, "--plink_prefix",
  help = "plink prefix"
)
p <- add_argument(
  p, "--sample_group",
  help = "sample group"
)
p <- add_argument(
  p, "--outdir",
  help = "outdir"
)

argv <- parse_args(p)

df <- read.delim(paste(argv$plink_prefix, 'eigenvec', sep="."), sep=" ", header=F, col.names = c("id", "sample_id", "PC1", "PC2", "PC3", "PC4", "PC5"))
eigenval <- read.delim(paste(argv$plink_prefix, 'eigenval', sep="."), header = F)

pca.explained <- eigenval$V1 / sum(eigenval$V1)

group.df <- read.delim(argv$sample_group)

plot_df <- merge(df, group.df, by = "sample_id")
groupNumber <- length(unique(group.df$group))

plot.colors <- brewer.pal(9, 'Set1')
group.colors <- colorRampPalette(plot.colors)(groupNumber)

xlab_text <- str_c('PCA1 ', '(', round(pca.explained[1] * 100,2), '%)')
ylab_text <- str_c('PCA2 ', '(', round(pca.explained[2] * 100,2), '%)')


p <- ggplot(plot_df, aes(x=PC1, y=PC2, color=group)) + geom_point() + theme_classic() + 
  theme(axis.title = element_text(size=rel(1.2)), axis.text = element_text(size = rel(1))) +
  scale_color_manual(values = group.colors) +
  guides(color=guide_legend(title="")) + xlab(xlab_text) + ylab(ylab_text)
p

ggsave(file.path(argv$outdir, 'pca.pdf'), width = 10, height = 8)
ggsave(file.path(argv$outdir, 'pca.png'), width = 10, height = 8)
