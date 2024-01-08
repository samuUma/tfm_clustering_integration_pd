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

# combat batch
library("sva")
mat <- ComBat(dat = t(PLA[,-1459]), 
              batch = PLA$Dataset, 
              mod = NULL)

mat=mat[,intersect(colnames(filteredWB),colnames(mat))]
# error in pacient "PP-60024-BLM0T1"
#mat=mat[,-c(294)]
#tcpms=tcpms[-c(294),]
set.seed(123)
# icluster with pla
library(iClusterPlus)
for(k in 3:8){
  plaD1.fit=iClusterPlus(dt1=(tcpms),dt2=t(mat),
                         type=c("gaussian","gaussian"),K=k,
                         #lambda=c(0.02,0.02),
                         #scale.lambda=c(1,1),
                         maxiter=10)
  save(plaD1.fit, file=paste("plaD1.fit",k,".Rdata",sep=""))
}

output=alist()
files=grep("plaD1.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=plaD1.fit
}

bic=c()
dev_ratio=c()
for(i in 1:length(output)){
  bic=c(bic,output[[i]]$BIC)
  dev_ratio=c(dev_ratio,output[[i]]$dev.ratio)
}

# only genes and proteins
genesprot=merge(x=tcpms,y=t(mat),by="row.names")
rownames(genesprot)=genesprot$Row.names
genesprot=genesprot[-1]
genesprot$icluster=output[[2]]$clusters

