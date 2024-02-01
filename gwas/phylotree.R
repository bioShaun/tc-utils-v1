suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(argparser))

options(stringsAsFactors = F)

p <- arg_parser("phylo tree plot")
p <- add_argument(
  p, "--tree",
  help = "input tree"
)
p <- add_argument(
  p, "--output",
  help = "output prefix"
)

argv <- parse_args(p)

tree <- read.tree(argv$tree)
# tree <- read.tree('./selected_sample.tree')



linear_width <- (0.2 * tree$Nnode) + 6
linear_height <- linear_width * 0.8

text_size = 5 - ceiling(log10(tree$Nnode))

if (tree$Nnode <= 50) {
  linear_tree <- ggtree(tree, branch.length="none") + geom_tiplab(size=text_size) + geom_rootedge()
  original_x_xmax <- max(linear_tree$data$x)
  linear_x_max = original_x_xmax + original_x_xmax * 2 / (linear_width - 2)
  linear_tree <- ggtree(tree, branch.length="none") + geom_tiplab(size=text_size) + geom_rootedge() + xlim(0, linear_x_max)
  ggsave(linear_tree, filename = paste(argv$output, 'linear.png', sep = '.'), dpi = 300, width = linear_width, height = linear_height)
  ggsave(linear_tree, filename = paste(argv$output, 'linear.pdf', sep = '.'), width = linear_width, height = linear_height)
}
circular_height <- (0.015 * tree$Nnode) + 5
circular_tree <- ggtree(tree, branch.length="none", layout='circular') + geom_tiplab(size=text_size) + geom_rootedge()
original_x_xmax <- max(circular_tree$data$x)
circular_x_max = original_x_xmax + original_x_xmax * 2 / (circular_height - 2)
circular_tree <- ggtree(tree, branch.length="none", layout='circular') + geom_tiplab(size=text_size) + geom_rootedge()+ xlim(0, circular_x_max)
ggsave(circular_tree, filename = paste(argv$output, 'circular.png', sep = '.'), dpi = 300, width = circular_height, height = circular_height)
ggsave(circular_tree, filename = paste(argv$output, 'circular.pdf', sep = '.'), width = circular_height, height = circular_height)



