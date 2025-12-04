library(ape)

x <- list.files(pattern="*.con.tre")

# read trees as nexus and write as newick
for(a in 1:length(x)) {
	a_rep <- read.nexus(x[a])
	if(a == 1) {
		write.tree(a_rep, file="camponotus_mrbayes.trees")
	} else {
		write.tree(a_rep, file="camponotus_mrbayes.trees", append=T)
	}
}

