# genewise ml plot
#
# created 2018-03-02

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]

#inputfile = "~/git/heteropecilly/simion2017_per_gene_lnl.tab"

genelnldat = read.table(inputfile, header=TRUE, sep="\t")

colorset = c("#fc8d62aa", "#66c2a5aa", "#377eb888")

t1nln = genelnldat[,2]
t2nln = genelnldat[,3]
t3nln = genelnldat[,4]

rowmaxs = apply( genelnldat[,2:4], 1, which.max)

minlnl = apply( genelnldat[,2:4], 1, min)
maxsitelnl = apply( genelnldat[,2:4], 1, max)
medianlnl = apply( genelnldat[,2:4], 1, median )
first_v_second_lnl = maxsitelnl - medianlnl

xpositions = 1:length(genelnldat[,1])
ymax = max(first_v_second_lnl)

stronggenes = (first_v_second_lnl >= (ymax/2))

outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=28, height=7)
par(mar=c(4.5,4.5,1,1))
plot(xpositions, first_v_second_lnl, type="h", pch=16, col=colorset[rowmaxs], xlim=range(xpositions), ylim=c(0,max(pretty(ymax))), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="Gene", axes=FALSE, cex.lab=1.4)
axis(2, cex.axis=1.3)
axis(1, at=seq(0,1700,100), labels=seq(0,1700,100) , cex.axis=1.3)
abline(h=c(5), lwd=1, lty=2)
text(0,28,"Favors T1", cex=2, col="#fc8d62", pos=4)
text(0,25,"Favors T2", cex=2, col="#66c2a5", pos=4)
text(0,22,"Favors T3", cex=2, col="#377eb8", pos=4)
text( xpositions[stronggenes], first_v_second_lnl[stronggenes]+1, genelnldat[,1][stronggenes])

dev.off()
