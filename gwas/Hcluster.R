
suppressMessages(library(gdsfmt))
suppressMessages(library(SNPRelate))
suppressMessages(library(argparser))

options(stringsAsFactors = F)

p <- arg_parser("vcf Hcluster")
p <- add_argument(
  p, "--vcf",
  help = "vcf"
)
p <- add_argument(
  p, "--outdir",
  help = "output dir"
)

argv <- parse_args(p)


#vcf to GDS
snpgdsVCF2GDS(argv$vcf, "my.gds")
genofile <- snpgdsOpen("my.gds")

#dendogram
dissMatrix  <-  snpgdsDiss(genofile,snp.id=NULL, autosome.only=F, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=10, verbose=TRUE)
snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=sample.id[sample.id != DZ55-2018], need.mat=TRUE)
cutTree <- snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,col.outlier=red, col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, verbose=TRUE)

sample_cluster <- data.frame(sample_id=cutTree$sample.id, group=cutTree$samp.group)

group_path = file.path(argv$outdir, "sample.group.tsv")

write.table(sample_cluster, file=group_path, sep="\t", quote=F,row.names=F)

#svg(filename = "diff.snp.hcluster.svg",width = 10,height = 16)
#snpgdsDrawTree(cutTree, main = "",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
#               y.label.kinship=T,leaflab="perpendicular", y.label=0.05)
#dev.off()
