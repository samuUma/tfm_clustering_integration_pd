#filtering2
wb2=wb[,intersect(rownames(PLA),colnames(wb))]
# expression normalization first
library("edgeR")
# get cpm counts:
cpms <- cpm(wb2)
# filter 80%
genesCounts <- apply(cpms, 1, function(x) {sum(x>1)})
genesProp <- sapply(genesCounts, function(x) {x > ncol(cpms)*0.8})
filteredcpms=cpms[genesProp,]
# function to filter by coefficient of variation
expression_varFilter = apply(filteredcpms, 1, function(x) {sd(x, na.rm = T)/abs(mean(x, na.rm = T)) > 0.5})
expression_filtered = filteredcpms[expression_varFilter,]

filteredWB <- wb2[rownames(expression_filtered),]

y <- DGEList(counts = filteredWB)
# calculate TMM normalization factors:
y <- calcNormFactors(y,method = "TMM")
# get the normalized counts:
cpms <- cpm(y,log = TRUE)
tcpms=t(cpms)
