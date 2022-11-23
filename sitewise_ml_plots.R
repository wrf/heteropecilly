# sitewise_ml_plots.R

#inputfile = "~/est/shen2017/Metazoan_phylogeny/genewise_lnl/RAxML_perSiteLLs.Whelan_D16_Choanoflagellata_site_lk.tab"
inputfile = "~/est/shen2017/Metazoan_phylogeny/genewise_lnl/RAxML_perSiteLLs.Whelan_D16_Opisthokonta_site_lk.tab"
inputfile = "~/git/heteropecilly/RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab"
inputfile = "~/git/heteropecilly/RAxML_perSiteLLs.simion2017_para-pori.tab"

sitelldat = read.table(inputfile, header=TRUE, sep="\t")

colorset = c("#fc8d62aa", "#66c2a5aa", "#377eb888")

t1nln = sitelldat[,2]
t2nln = sitelldat[,3]
t3nln = sitelldat[,4]

t1vt2diff = t1nln-t2nln

#x = t1vt2diff[10543:11140]

max(abs(t1nln - t2nln))

maxbin = 240
bins = 0:maxbin

rowmaxs = apply( sitelldat[,2:4], 1, which.max)

minlnl = apply( sitelldat[,2:4], 1, min)
maxsitelnl = apply( sitelldat[,2:4], 1, max)
medianlnl = apply( sitelldat[,2:4], 1, median )
first_v_second_lnl = maxsitelnl - medianlnl

hist(first_v_second_lnl, breaks=seq(0,5,0.01), ylim=c(0,200))

strongsites = first_v_second_lnl > 0.25
vstrongsites = first_v_second_lnl > 0.1

sum( t1nln[strongsites] )
sum( t2nln[strongsites] )
sum( t3nln[strongsites] )

sum( t1nln[strongsites][rowmaxs[strongsites]==1] )
sum( t2nln[strongsites][rowmaxs[strongsites]==2] )
sum( t3nln[strongsites][rowmaxs[strongsites]==3] )

outputfile = gsub("([\\w/]+)\\....$","\\1.1v2.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=40, height=7)
#plot(sitelldat[,1], first_v_second_lnl, type="h", lwd=1, col=c("#fc8d62", "#66c2a5", "#377eb8")[rowmaxs], xlim=range(sitelldat[,1]), ylim=c(0,7), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="Alignment position")
plot(sitelldat[,1], first_v_second_lnl, type="p", pch=16, col=colorset[rowmaxs], xlim=range(sitelldat[,1])*.99, ylim=c(0,5), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="Alignment position")
#text( sitelldat[,1][vstrongsites], first_v_second_lnl[vstrongsites], sitelldat[,1][vstrongsites])
#text( sitelldat[,1][vstrongsites], first_v_second_lnl[vstrongsites]+0.1, sitelldat[,1][vstrongsites])
abline(h=c(0.25), lwd=1, lty=2)
text(0,5,"T1-Ctenophora-sister", cex=2, col="#fc8d62", pos=4)
text(0,4.6,"T2-Porifera-sister", cex=2, col="#66c2a5", pos=4)
text(0,4.2,"T3-Porifera-Cteno", cex=2, col="#377eb8", pos=4)

#text(0,5,"T1-Porifera-sister", cex=2, col="#fc8d62", pos=4)
#text(0,4.6,"T2-Polyphyletic-zoa", cex=2, col="#66c2a5", pos=4)
#text(0,4.2,"T3-Paraphyletic-sponge", cex=2, col="#377eb8", pos=4)
dev.off()


outputfile = gsub("([\\w/]+)\\....$","\\1.t1vt3.pdf",inputfile,perl=TRUE)

t1vt3diff = t1nln-t3nln
sitecolor = as.integer(t1vt3diff > 0) + 1
pdf(outputfile, width=40, height=7)
plot(sitelldat[,1], t1vt3diff, type="h", lwd=1, col=c("#fc8d62", "#66c2a5")[sitecolor], xlim=range(sitelldat[,1]), ylim=c(-7,7), ylab=expression(paste(Delta,"ln(L)",sep="")), xlab="Alignment position")
#text( sitelldat[,1][vstrongsites], first_v_second_lnl[vstrongsites], sitelldat[,1][vstrongsites])
#text( sitelldat[,1][vstrongsites], first_v_second_lnl[vstrongsites]+0.1, sitelldat[,1][vstrongsites])
abline(h=c(0.25), lwd=1, lty=2)
text(0,5,"T1-Monophyletic-sponge", cex=2, col="#fc8d62", pos=4)
#text(0,4.6,"T2-Polyphyletic-zoa", cex=2, col="#66c2a5", pos=4)
text(0,4.2,"T2-Paraphyletic-sponge", cex=2, col="#66c2a5", pos=4)
#text(0,4.2,"T3-Paraphyletic-sponge", cex=2, col="#377eb8", pos=4)

dev.off()


dSLS = ( (t1nln - t2nln) + (t1nln - t3nln) + (t2nln - t3nln) ) / 3
absdsls = abs(dSLS)

sum(dSLS[rowmaxs==1])
sum(dSLS[rowmaxs==2])
sum(dSLS[rowmaxs==3])


t1hist = hist(-1 * t1nln, breaks=bins, plot=FALSE)
t2hist = hist(-1 * t2nln, breaks=bins, plot=FALSE)
t3hist = hist(-1 * t3nln, breaks=bins, plot=FALSE)

outputfile = gsub("([\\w/]+)\\....$","\\1.hist.pdf",inputfile,perl=TRUE)

#outputfile = "~/est/shen2017/Metazoan_phylogeny/genewise_lnl/RAxML_perSiteLLs.Whelan_D16_Choanoflagellata_site_lk.hist.pdf"

pdf(outputfile, width=7, height=7)

plot(0,0,type='n', xlim=c(0,maxbin), ylim=c(0,5000), xlab="-ln(L)", ylab="Number of sites")
lines(1:maxbin, t1hist$counts, lwd=3, col=colorset[1])
lines(1:maxbin, t2hist$counts, lwd=3, col=colorset[2])
#lines(0:maxbin, t3hist$counts, lwd=3, col=colorset[3])
legend(100,5000,legend=c("T1-Ctenophora-sister","T2-Porifera-sister","T3-Ahifozoa"), col=colorset, lwd=4)
dev.off()

dnln = t1nln-t2nln
absdnln = abs(dnln)

dnlnhist = hist(dnln, breaks=seq(-7,7,0.1))
absdnlnhist = hist( absdnln, breaks=seq(0,7,0.1))

outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)

pdf(outputfile, width=15, height=7)
plot(sitelldat[,1], dSLS, type="h", lwd=1, col=c("#fc8d62", "#66c2a5", "#377eb8")[rowmins], xlim=range(sitelldat[,1]), ylim=c(-5,5), ylab="-ln(L)", xlab="Alignment position")
text( sitelldat[,1][absdsls>0.5], dSLS[absdsls>0.5], sitelldat[,1][absdsls>0.5])
abline(h=c(0.5,-0.5), lwd=1, lty=2)
text(0,-4,"T1-Ctenophora-sister", cex=2, col="#fc8d62", pos=4)
text(0,4,"T2-Porifera-sister", cex=2, col="#66c2a5", pos=4)
text(0,5,"T3-Porifera-Cteno", cex=2, col="#377eb8", pos=4)
dev.off()


pdf(outputfile, width=15, height=7)
plot(0,0,type='n', xlim=range(sitelldat[,1]), ylim=c(-5,5), ylab="-ln(L)", xlab="Alignment position")
lines(sitelldat[,1], dnln, lwd=1, col="#984ea3")

#lines(sitelldat[,1], t2nln, lwd=1, col=colorset[2])
#lines(0:150, t3hist$counts, lwd=3, col=colorset[3])
#legend(0,150,legend=c("T1-Ctenophora-sister","T2-Porifera-sister","T3-Ahifozoa"), col=colorset, lwd=4)

text(0,-4,"T1-Ctenophora-sister", cex=2, col="#fc8d62", pos=4)
text(0,4,"T2-Porifera-sister", cex=2, col="#66c2a5", pos=4)

text( sitelldat[,1][absdnln>1], dnln[absdnln>1], sitelldat[,1][absdnln>1])

abline(h=c(0.5,-0.5), lwd=1, lty=2)

dev.off()