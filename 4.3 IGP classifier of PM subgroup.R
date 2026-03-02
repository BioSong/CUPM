load("./CUPM/data/RNAseq/GSE183202/GSE183202_Exp_Cli.RData")
exp1=log2(fpkm1+1)
#------------------filter out 低表达的gene
index1=rowMeans(exp1<1)
exp1=exp1[which(index1<0.5),]
cli1[cli1==""]=NA
cli11=cli1[!is.na(cli1$groups),]
exp11=exp1[,!is.na(cli1$groups)]





