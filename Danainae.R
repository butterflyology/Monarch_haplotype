### Code for Danainae chromosomal evolution paper
## Hamm 2017 

set.seed(89623578)
library("phytools")
library("geiger")
library("phytools")

(sess <- sessionInfo())

setwd("~/Desktop/Projects/Monarch_haplotype/")

# save(list = ls(), file = "Danainae_data.RData")
# load("Danainae_data.RData")


Dan <- read.nexus("Data/Bayes/Danaini_cat.nex.con.tre")
plot(Dan)
is.rooted(Dan)
is.ultrametric(Dan)
sort(branching.times(Dan))

lam <- 10^(-1:6)
cv <- sapply(lam, function(x) sum(attr(chronopl(Dan , lambda = x, CV = TRUE), "D2")))
cv
plot(x = lam, y = cv, pch = 19, ylab = "cross-validation score", xlab = expression(paste(lambda)), las = 1, cex = 1.5) # lowest CV is alpha 1e6, suggests that branches each have different rates


Dan2 <- chronopl(phy = Dan, lambda = 1e6, CV = TRUE, eval.max = 1e3, iter.max = 1e4)
is.ultrametric(Dan2)
Dan2 <- root(phy = Dan2, outgroup = 44, resolve.root = TRUE)
is.rooted(Dan2)

# pdf(file = "Danainae.pdf", bg = "white")
plot(Dan2, type = "phylogram", use.edge.length = TRUE, edge.width = 2, cex = 0.8, no.margin = TRUE)
# dev.off()

Dan_data <- read.csv("Danainae_karyotypes2.csv", header = TRUE, row.names = 1)
str(Dan_data)
Dan_data


Dan <- root(phy = Dan, outgroup = "Libythea_celtis", resolve.root = TRUE)
Dan_norm <- treedata(phy = Dan2, data = Dan_data, sort = TRUE, warnings = TRUE)
is.rooted(Dan_norm$phy)
plot.phylo(Dan_norm$phy)
#ape::write.tree(Dan_norm$phy, file = "Dan_norm.tre", digits = 10, tree.names = FALSE)
nkaryo <- length(Dan_norm$phy$tip.label)
Covs <- Dan_norm$data[Dan_norm$phy$tip.label,, drop = FALSE]

# segments(rep(1.6, nkaryo), 1:nkaryo, rep(1.6, nkaryo) + (0.005 * Covs), 1:nkaryo, lwd = 2, col = "black")

dotTree(Dan_norm$phy, Covs, method = "plotTree", legend = TRUE, standardize = FALSE, data.type = "continuous", fsize = 1.5)

c1 <- colorRampPalette(c("blue", "red"))
c1(9)

Covs_cols <- Covs
Covs_cols[Covs == 22] <- "#0000FF"
Covs_cols[Covs == 28] <- "#1F00DF"
Covs_cols[Covs == 29] <- "#3F00BF"
Covs_cols[Covs == 30] <- "#5F009F"
Covs_cols[Covs == 31] <- "#7F007F"
Covs_cols[Covs == 32] <- "#9F005F"
Covs_cols[Covs == 33] <- "#BF003F"
Covs_cols[Covs == 36] <- "#DF001F"
Covs_cols[Covs == 40] <- "#FF0000"
Covs_cols

dotTree(Dan_norm$phy, Covs, method = "plotTree", legend = TRUE, standardize = FALSE, data.type = "discrete")


# Figure for pre-print
Covs2 <- as.matrix(Covs)[, 1]
seq1 <- seq(from = 20, to = 36)
ancs <- c(31, 31, 31, 29, 31, 31, 31, 30, 30, 30, 30, 28, 33, 30, 30, 31, 31)
sa <- cbind(seq1, ancs)

# pdf(file = "Images/Dan_phylo.pdf", bg = "white")
 tiff(file = "Images/Dan_phylo.tiff", bg= "white", height = 720, width = 520, units = "px", compression = "none")
plot.phylo(Dan_norm$phy, align.tip.label = TRUE, type = "phylogram", use.edge.length = TRUE, edge.width = 2, root.edge = TRUE, no.margin = TRUE, adj = 0)
# tiplabels(text = Covs2, adj = -9.3, frame = "none", bg = "white")
tiplabels(text = Covs2, adj = 1.2, frame = "rect", bg = "white")
nodelabels(text = sa[, 2], frame = "rect", bg = "white")
 dev.off()






phylo.heatmap(Dan_norm$phy, Covs)




plotTree.wBars(Dan_norm$phy, Covs2, type = "phylogram", method = "plotTree", scale = 0.005, tip.labels = TRUE) # not a fan of this look, not enough to distinguish.

plotTree.barplot(Dan_norm$phy, Covs2, args.barplot = list(xlim = c(20, 40)))




Dan_pruned <- treedata(phy = Dan2, data = Dan_data, sort = TRUE, warnings = TRUE)
is.rooted(Dan_pruned$phy)
plot(Dan_pruned$phy)
Dan3 <- root(phy = Dan_pruned$phy, outgroup = "Libythea_celtis", resolve.root = TRUE)
plot(Dan3)

#write.nexus(Dan_pruned$phy, file = "Dan2.nxs")
#ape::write.tree(phy = Dan_pruned$phy, file = "Danaini.tre", digits = 10, tree.names = FALSE)



plotTree(Dan_norm$phy, type = "phylogram", use.edge.length = TRUE, edge.width = 2, cex = 0.8, use.edge.length = TRUE)


par(fig = c(x1 = 0.80, x2 = 0.9, y1 = 0.8, y2 = 0.85), new = TRUE)
hist(x = density(c(30, 30)), xlim = c(22, 40), ylim = c(0, 1), main = "", xlab = "", ylab = "", las = 1, breaks = 30, col = "grey", yaxt = "n", xaxt = "n")


plot.new()
chrys <- 30

segments(x0 = 0.555, y0 = 0.625, x1 = 1.8, y1 = 1)

par(new = TRUE)


#segments(rep(1.6, nkaryo), 1:nkaryo, rep(1.6, nkaryo) + (0.005 * Covs), 1:nkaryo, lwd = 2, col = "black")
#plot.new()
#par(mar = c(1, 1, 0, 0))

par(fig = c(0.9, 1.7, 0.8, 1.2), new = TRUE)


# Run iteRates, get richness data for tips
library("iteRates")
sort(branching.times(Dan_norm$phy))

Dv <- vcvPhylo(Dan_norm$phy)
Dv

Dan.exp <- comp.subs(Dan_norm$phy, mod.id = c(1, 0, 0, 0))
Dan.rv <- comp.subs(Dan_norm$phy, mod.id = c(0, 0, 0, 1))
Dan.Wei <- comp.subs(Dan_norm$phy, mod.id = c(0, 1, 0, 0))

par(mfrow = c(1, 3))
color.tree.plot(Dan.exp, Dan3, p.thres = 0.05, root.edge = TRUE)
color.tree.plot(Dan.Wei, Dan3, p.thres = 0.05, root.edge = TRUE) #entirely due to taxon sampling, methinks. 
color.tree.plot(Dan.rv, Dan3, p.thres = 0.05, root.edge = TRUE)


