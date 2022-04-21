#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
#png(argv[2],width = 1280, height = 480, res=326)
png(argv[2],height =5,width = 14,units = "in",res=326)

gwasResults1 <- read.table(argv[1], header=T)
library('qqman')
par(mgp=c(2.2,0.6,0),las=1,cex.lab=1.4,cex.axis=1.2,lend=1,tcl=-0.3,mfrow=c(1,2))

manhattan(gwasResults1,chr="CHR", bp="BP", p="P", snp="SNP",main=argv[4],col = c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62"),suggestiveline = F,genomewideline=-log(as.numeric(argv[3]),10),xaxs='i',yaxs='i',cex=0.8) #annotatePval = 0.0001 最高点标出snp名字 main=标题(name of metabolites)

qq(gwasResults1$P,col='blue',cex=0.8,bty='n',xaxs='i',yaxs='i')
dev.off()
