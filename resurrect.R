#!/usr/bin/R

args <- commandArgs(TRUE)

library(phangorn)

########################################################
# Arguments order

# args[1] -- sequence file
# args[2] -- newick tree file
# args[3] -- sequence file format

########################################################

tree <- read.tree(args[2])
chrseq <- read.phyDat(args[1], format=args[3])

fit <- pml(tree,chrseq)
temp <- ancestral.pml(fit, type="ml")

ptab = array(dim=c(length(temp), length(temp[[1]])/4))

locmax <- function(list) {
	if (list[1]==list[2] && list[3]==list[4])
	{
		return (18)
	}	
	else
	{
		for (i in seq(1, length(list)))
		{
			if (max(list) == list[i])
			{
				return (i)
			}
		}
	}	
}


for (h in 1:length(temp)) {
	for (j in 1:(length(temp[[1]])%/%4))
	{
		ptab[h, j] <- locmax(temp[[h]][j,])
	}
}

write.table(ptab, "/tmp/rdata.dat", col.name=F, row.names=names(temp))
