
## get data to plot
setwd("/home/thibaut/Dropbox/teaching/MSc-Imperial/phylogenetics/lecture/figs/R")

library(ape)
library(adegenet)


#### WHERE'S THE ROOT? ####
data(woodmouse)
col <- spectral(15)

tre <- nj(dist.dna(woodmouse))

## basic plot
plot(tre, show.tip=FALSE)
tiplabels(col=col, pch=20, cex=5)
tiplabels(1:15, frame="n")


## try plotting with different rootings
pdf("../6trees.pdf", width=21, height=14)

par(mfrow=c(2,3), mar=rep(.1,4))

plot(root(tre,1), show.tip=FALSE)
tiplabels(col=col, pch=20, cex=13)
tiplabels(1:15, frame="n",cex=3)

plot(root(tre,2), show.tip=FALSE)
tiplabels(col=col, pch=20, cex=13)
tiplabels(1:15, frame="n",cex=3)

plot(root(tre,8), show.tip=FALSE)
tiplabels(col=col, pch=20, cex=13)
tiplabels(1:15, frame="n",cex=3)

plot(root(tre,1), "fan", show.tip=FALSE)
tiplabels(col=col, pch=20, cex=13)
tiplabels(1:15, frame="n",cex=3)

plot(root(tre,2), "fan", show.tip=FALSE)
tiplabels(col=col, pch=20, cex=13)
tiplabels(1:15, frame="n",cex=3)

plot(tre, "unrooted", show.tip=FALSE)
tiplabels(col=col, pch=20, cex=13)
tiplabels(1:15, frame="n",cex=3)

dev.off()




#### MOLECULAR CLOCK? ####
set.seed(1)
tree <- rtree(100)
tip.date <- sample(1:30, 100, replace=TRUE)
tree <- ladderize(rtt(tree, tip.date))

library(adephylo)
d.root <- distRoot(tree)
t.root <- round(3.5*d.root + rnorm(100,sd=20))
t.root <- abs(t.root-t.root[which.min(d.root)])



## figure 1
pdf("../molecularClock1.pdf",width=14,height=7)

par(mfrow=c(1,2))

plot(tree)
axisPhylo()

dev.off()

## figure 2
pdf("../molecularClock2.pdf",width=14,height=7)

par(mfrow=c(1,2))

plot(tree)
axisPhylo()

plot(t.root, d.root, col=transp("black"), cex=3,
     xlab="Time to the root", ylab="Mutations to the root", pch=20)
abline(lm(d.root~-1+t.root), col="salmon",lwd=2)

anova(lm(d.root~-1+t.root))

dev.off()




#### UNCERTAINTY? ####
data("hivtree.newick") 
hivtree.newick

## generate file "hivtree.phy" in working directory
cat(hivtree.newick, file = "hivtree.phy", sep = "\n")
tree.hiv <- read.tree("hivtree.phy") # load tree
unlink("hivtree.phy") # delete the file "hivtree.phy"


pdf("../hivtree1.pdf",width=14,height=7)

par(mfrow=c(1,2))

plot(tree.hiv, main="193 HIV-1 sequences from DRC\n(Strimmer & Pybus 2001)", show.tip=FALSE)
axisPhylo()

dev.off()


pdf("../hivtree2.pdf",width=14,height=7)

par(mfrow=c(1,2))

plot(tree.hiv, main="193 HIV-1 sequences from DRC\n(Strimmer & Pybus 2001)", show.tip=FALSE)
axisPhylo()

plot(di2multi(tree.hiv, tol=0.01), main="Collapsed tree (threshold length 0.01)", show.tip=FALSE)
axisPhylo()

dev.off()

