setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Code")

# if want population data:
#mutations <- read.csv("../Data/mutations_populations.csv")

# if want clonal data:
mutations=read.csv("../Data/mutations_clones.csv")
colnames(mutations)[1] <- "Gene"


mutations <- mutations[, c("Gene", "Population")]

tested.genes <- c()
for (gene in unique(mutations$Gene)) {
  aa <- subset(mutations, Gene == gene)
  pp <- unique(aa[, "Population"])
  if (length(pp) > 2) {
    tested.genes <- rbind(tested.genes, data.frame(gene, pp))
  }
}

tested.genes$group <- substr(tested.genes$pp, 1, 4)
tested.genes$category <- as.factor(gsub("[FS]", "D", tested.genes$group))

tt <- table(tested.genes[, c("gene", "category")])

pval <- c()
for (i in 1:dim(tt)[1]) {
  a <- tt[i, ]
  b <- c(8, 8, 4, 4) - a
  pval <- c(pval, fisher.test(cbind(a, b))$p.val)
}

p.BH <- p.adjust(pval, method = "BH")

genes <- cbind(tt, pval, p.BH)

library(pander)
pander(genes)
