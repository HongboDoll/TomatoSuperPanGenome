#!/usr/bin/env Rscript

argv = commandArgs(T)

pdf(argv[2], 12, 5)
par(mfrow=c(3,4), mgp=c(2.7,0.7,0),las=1,mar=c(1,4,1,1),oma=c(1,1,1,1),cex.lab=1.5,cex.axis=1.3)
d1 <- read.table(argv[1],head=F)

for (chr in c('01','02','03','04','05','06','07','08','09','10','11','12')){
#	par(new=T)  
	chrr <- paste('Sly',as.character(chr), sep='')
	sub <- subset(d1,V1==chrr)
	if (chr == '01' || chr == '05' || chr == '09'){
	#if (as.character(chr)=='01'){
		plot(sub[,2],sub[,4],type='h',ylim=c(0,100),xlim=c(0,100000000),ylab='Normalized coverage',col='blue',axes=F,lwd=1, xaxs='i',yaxs='i') # reads per genomic content
	} else{
                 plot(sub[,2],sub[,4],type='h',ylim=c(0,100),xlim=c(0,100000000),ylab='',col='blue',axes=F,lwd=1, xaxs='i',yaxs='i') # reads per genomic content

	}
	axis(1,at=c(0,100000000),labels=c('0','100 Mb'),lwd=1.5,lwd.ticks=1.5,tcl=-0.3)
	axis(2,at=seq(0,100,20),labels=seq(0,100,20),lwd=1.5,lwd.ticks=1.5,tcl=-0.3)
#	par(new=T)
}

