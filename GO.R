#### GO enrichment analysis (hypergeometric test) ####
args=commandArgs(trailingOnly = TRUE)
input_and_output_dir <- args[1]
cor.global.file <- file.path(input_and_output_dir, "cor_global_only.txt")
cor.local.file <- file.path(input_and_output_dir, "cor_local_only.txt")

cor.global=read.delim(cor.global.file, header = T,stringsAsFactors = F)
cor.local=read.delim(cor.global.file, header = T,stringsAsFactors = F)
sgd=read.delim("./gene_association.sgd", header = F,stringsAsFactors = F,sep = "\t",quote = "")
goterm=unique(sgd$V5)

setwd(input_and_output_dir)
#enrichment test-global (top20)
GO_result_global=c()
for (i in 1:length(goterm)) {
  m=length(unique(sgd[sgd$V5==goterm[i],3]))
  k=length(cor.global[order(cor.global[,3])[1:20],1])
  x=length(intersect(k,m))
  n=length(union(cor.global[,1],cor.local[,1]))-m
  p.val = phyper(x,m,n,k,lower.tail = F)
  pre=c(goterm[i],p.val)
  GO_result_global=rbind(GO_result_global,pre)
}
colnames(GO_result_global)=c("GO_term","p.value")
write.table(GO_result_global,file = "./GO_result_global.txt",quote = F,row.names = F)


#enrichment test-local (top20)
GO_result_local=c()
for (i in 1:length(goterm)) {
  m=length(unique(sgd[sgd$V5==goterm[i],3]))
  k=length(cor.local[order(cor.local[,3])[1:20],1])
  x=length(intersect(k,m))
  n=length(union(cor.global[,1],cor.local[,1]))-m
  p.val = phyper(x,m,n,k,lower.tail = F)
  pre=c(goterm[i],p.val)
  GO_result_local=rbind(GO_result_local,pre)
}
colnames(GO_result_local)=c("GO_term","p.value")
write.table(GO_result_local,file = "./GO_result_local.txt",quote = F,row.names = F)
