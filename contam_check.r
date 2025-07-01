# input args <- input file, base name of individual
args <- commandArgs(trailingOnly = TRUE)

# no scientific notation
options(scipen=999)

# read in input file
x <- read.table(args[1], sep=",", stringsAsFactors=F)

base_name <- args[2]

output_name <- paste0(base_name, ".pdf")

# ref allele freq
x1 <- x[,1] + x[,2]

# alt allele freq
x2 <- x[,3] + x[,4]

# get allele frequency of allele 1
maf <- x1 / (x1 + x2)

# choose the frequency of the minor allele
maf[maf > 0.5] <- 1 - maf[maf > 0.5]

pdf(output_name, width=5, height=5)

hist(maf, breaks=seq(from=0,to=0.5,by=0.02), main=base_name, xlab="maf coverage proportion")

dev.off()
