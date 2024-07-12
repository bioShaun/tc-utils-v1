suppressMessages(library(pophelper))
suppressMessages(library(stringr))
suppressMessages(library(argparser))

options(stringsAsFactors = F)

p <- arg_parser("structure plot")
p <- add_argument(
  p, "--structure_dir",
  help = "structure dir"
)
p <- add_argument(
  p, "--outdir",
  help = "output dir"
)

argv <- parse_args(p)

structureFiles <- list.files(argv$structure_dir, full.names=TRUE)
qFiles <- structureFiles[str_ends(structureFiles, 'meanQ')]
q <- readQ(qFiles)
plotQ(q, exportpath=argv$outdir, sortind="all")
