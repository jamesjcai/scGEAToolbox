require(node2vec)
gene_edges<-read.csv("input.txt")
#emb<-node2vecR(gene_edges,p=2,q=1,num_walks=5,walk_length=5,dim=10)
emb<-node2vecR(gene_edges,directed=TRUE)
write.csv(emb,file="output.txt")
