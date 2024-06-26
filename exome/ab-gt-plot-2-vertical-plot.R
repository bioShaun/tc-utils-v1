suppressMessages(library(ggplot2))
suppressMessages(library(argparser))

p <- arg_parser("for vertical origin plot")
p <- add_argument(p, "--plot_file", help = "plot file")
p <- add_argument(p, "--chr_size", help = "chr size file")
p <- add_argument(p, "--p1_name", help = "p1 name")
p <- add_argument(p, "--p2_name", help = "p2 name")
p <- add_argument(p, "--p1_color", help = "p1 color")
p <- add_argument(p, "--p2_color", help = "p2 color")
p <- add_argument(p, "--both_color", help = "both color")
p <- add_argument(p, "--out_prefix", help = "output plot prefix")
argv <- parse_args(p)

raw.data <- read.delim(argv$plot_file, check.names = F)
chr_df <- read.table(argv$chr_size, header = F)
colnames(chr_df) <- c("CHROM", "chrom_length")

plot_color <- c(argv$p1_color, argv$p2_color, argv$both_color)
names(plot_color) <- c(argv$p1_name, argv$p2_name, "Both")

raw.data$END <- raw.data$POS + 100000
raw.data$CHROM <- factor(raw.data$CHROM, levels = chr_df$CHROM)
chr_df$CHROM<- factor(chr_df$CHROM, levels = chr_df$CHROM)
out_prefix <- argv$out_prefix
raw.data$origin <- factor(raw.data$origin, levels = c(argv$p1_name, argv$p2_name, "Both"))

p <- ggplot(raw.data) +
  geom_rect(aes(
    xmin = POS, xmax = END, ymin = 0,
    ymax = 1, fill = origin
  )) +  
  geom_rect(data = chr_df, aes(
    xmin = 0, xmax = chrom_length, ymin = 0,
    ymax = 1
  ), fill = "white", color='black', alpha=0) +    
  scale_fill_manual(values = plot_color) +
  facet_wrap(CHROM~., nrow = 1) +
  theme(strip.text.y = element_text(size=rel(.8), 
                                    face="bold",
                                    angle = 0,
                                    vjust=100
                                    ),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(colour = "white", fill = "white"),        
        axis.text.x =  element_blank()) +
  coord_flip()+ scale_x_reverse(expand = c(0, 0)) +
  guides(fill=guide_legend(title="")) 

chrom_num <- length(unique(chr_df$CHROM))

p_width = 5 * chrom_num  / 7
ggsave(paste(out_prefix, 'png', sep='.'),
        plot = p, width = p_width, height = 12,
        dpi = 300, type = "cairo")
ggsave(paste(out_prefix, 'pdf', sep='.'), 
       plot = p, width = p_width, height = 12, limitsize=F)
