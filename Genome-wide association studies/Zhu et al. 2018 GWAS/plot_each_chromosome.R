#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
pdf(argv[2], 7, 14)
a <- read.table(argv[1],header=F)

par(mgp=c(2.2,0.7,0),oma=c(0,0,0,0),mar=c(0,0,0,0),las=0,lend=1)

chr_col=rgb(255,30,49,max=255)
sub1_col=rgb(255,48,146,max=255)
sub2_col=rgb(76,170,232,max=255)
sub3_col=rgb(255,147,0,max=255)
sub4_col=rgb(165,108,189,max=255)
sub5_col=rgb(255,135,141,max=255)
sub6_col=rgb(95,95,95,max=255)
sub7_col=rgb(97,216,54,max=255)
sub8_col=rgb(22,231,207,max=255)
sub_col = c(sub1_col, sub2_col, sub3_col, sub4_col, sub5_col, sub6_col, sub7_col, sub8_col)

chr_len = c(97000000,56000000,67000000,66000000,69000000,52000000,71000000,66000000,72000000,69000000,55000000,68000000) ## length of each chromosome

plot(1,1,type='n',xlab='',ylab='',main='',xlim=c(0,8),ylim=c(0,97000000),axes=F) ### chr_len

rect(0,0,0.2,chr_len[as.numeric(argv[3])],border=chr_col,col=chr_col)

for (i in 1:length(a[,1]))
{
    y_pos = chr_len[as.numeric(argv[3])] - (a[i,2]+((a[i,3]-a[i,2])/2)) ### chr_len
    for (n in 4:11)
    {
        x_pos = 0.3 + (n-3)*0.25
        c = a[i,n]
        cc = 0
        if (c > 12 && c <= 20)
        {
            cc = 3
        } else if (c > 9 && c <= 12)
        { 
            cc = 2.5
        }else if (c > 6 && c <= 9)
        { 
            cc = 2
        }else if (c > 3 && c <= 6)
        { 
            cc = 1.5
        }else if (c > 0 && c <= 3)
        {
            cc = 1
        } else
	{
	    cc = 0
	}
        points(x_pos,y_pos,type='p',pch=19,col=sub_col[n-3],cex=cc)
        
    }

}

cc_list = c(3,2.5,2,1.5,1)
for (i in 1:5)
{
    points(4,10000000+(i-1)*2500000,type='p',pch=20,col=sub_col[1],cex=cc_list[6-i])
    
}
#for (i in 1:8)
#{
#    points(6,10000000+(i-1)*2500000,type='p',pch=20,col=sub_col[9-i],cex=cc_list[1])
#    
#}
#
