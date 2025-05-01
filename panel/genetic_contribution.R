suppressMessages(library(ggplot2))
suppressMessages(library(argparser))

p <- arg_parser("for vertical origin plot")
p <- add_argument(p, "--plot_file", help = "plot file")
p <- add_argument(p, "--chr_size", help = "chr size file")
p <- add_argument(p, "--p1_name", help = "p1 name")
p <- add_argument(p, "--p2_name", help = "p2 name")
p <- add_argument(p, "--offspring_name", help = "offspring name")
p <- add_argument(p, "--same_name", help = "same name")
p <- add_argument(p, "--p1_color", help = "p1 color")
p <- add_argument(p, "--p2_color", help = "p2 color")
p <- add_argument(p, "--offspring_color", help = "offspring color")
p <- add_argument(p, "--same_color", help = "same color")
p <- add_argument(p, "--out_prefix", help = "output plot prefix")
p <- add_argument(p, "--bin_size", help = "pos bin size", type="numeric", default=0.5)
argv <- parse_args(p)

raw.data <- read.delim(argv$plot_file, check.names = F)
chr_df <- read.delim(argv$chr_size)

plot_color <- c(argv$p1_color, argv$p2_color, argv$offspring_color, argv$same_color)
names(plot_color) <- c(argv$p1_name, argv$p2_name, argv$offspring_name, argv$same_name)
print(plot_color)

raw.data$POS_MB <- raw.data$pos / 1000000
raw.data$END_MB <- raw.data$POS_MB + argv$bin_size

raw.data$chrom <- factor(raw.data$chrom, levels = chr_df$chrom)
chr_df$chrom<- factor(chr_df$chrom, levels = chr_df$chrom)


chr_df$chrom_length_mb <- chr_df$chrom_size / 1000000


out_prefix <- argv$out_prefix
raw.data$origin <- factor(raw.data$origin, levels = c(argv$p1_name, argv$p2_name, argv$offspring_name, argv$same_name))


p <- ggplot(raw.data) +
  geom_rect(aes(
    xmin = POS_MB, xmax = END_MB, ymin = 0,
    ymax = 1, fill = origin
  )) +  
  geom_rect(data = chr_df, aes(
    xmin = 0, xmax = chrom_length_mb, ymin = 0,
    ymax = 1
  ), fill = "white", color='black', alpha=0) +    
  scale_fill_manual(values = plot_color) +
  facet_wrap(chrom~., nrow = 1) +
  theme(strip.text.y = element_text(size=rel(.8), 
                                    face="bold",
                                    angle = 0,
                                    vjust=100
                                    ),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(colour = "white", fill = "white")) +
  coord_flip()+ scale_x_reverse(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.5)) +
  guides(fill=guide_legend(title="")) +
  xlab("Location (Mb)")

#
# axis.text.y = element_blank(),

chrom_num <- length(unique(chr_df$chrom))

p_width = 5 * chrom_num  / 7
ggsave(paste(out_prefix, 'png', sep='.'),
        plot = p, width = p_width, height = 12,
        dpi = 300, type = "cairo")
ggsave(paste(out_prefix, 'pdf', sep='.'), 
       plot = p, width = p_width, height = 12, limitsize=F)
