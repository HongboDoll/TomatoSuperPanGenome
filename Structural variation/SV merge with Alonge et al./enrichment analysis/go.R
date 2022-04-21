#!/usr/bin/env Rscript
library("topGO")

argv<-commandArgs(TRUE)
b <- read.table(argv[2])
geneID2GO <- readMappings(argv[1])
geneNames <- names(geneID2GO)
myInterestingGenes <- as.character(b[,1])
geneList_nscore <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList_nscore) <- geneNames

GOdata_nscore_BP <- new("topGOdata",
                        ontology="BP",
                        allGenes=geneList_nscore,
                        annot=annFUN.gene2GO, gene2GO=geneID2GO,
                        nodeSize=5)


result <- runTest(
    GOdata_nscore_BP, 
    algorithm = "classic", 
    statistic = "fisher")

gtFis <- GenTable(GOdata_nscore_BP,classicFisher=result,orderBy="classic",ranksOf="classicFisher",topNodes = 100)

fdr <- p.adjust(p=gtFis[,"classicFisher"],method="fdr")

all <- paste(as.character(gtFis[,'Annotated']),'/',as.character(length(geneScore(GOdata_nscore_BP,  use.names = FALSE))), ' ', '(',as.character(round(gtFis[,'Annotated']*100/length(geneScore(GOdata_nscore_BP,  use.names = FALSE)),2)) ,'%)', sep='')
dis <- paste(as.character(gtFis[,'Significant']),'/',as.character(numSigGenes(GOdata_nscore_BP)), ' ', '(',as.character(round(gtFis[,'Significant']*100/numSigGenes(GOdata_nscore_BP),2)) ,'%)', sep='')

r <- cbind(gtFis[,1], gtFis[,2], dis, all, gtFis[,6], fdr, deparse.level=0)

write.table(r, file=argv[3], sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)




