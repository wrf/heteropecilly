# process raw heteropecilly data

infile = "~/git/heteropecilly/heteropecilly-v2/PAUL-90x341543v-C20.proba.tab"

hpdat1 = read.table(infile, header=TRUE, sep="\t")

sortedlnpip = sort(hpdat1[,3])
maxlnhp = max(hpdat1[,3])

bins = seq(0,45,0.2)

decilebreaks = sortedlnpip[floor(1:9/10*length(sortedlnpip))]
roundedbins = round(decilebreaks*5)/5

lnhphist = hist(hpdat1[,3], breaks=bins, plot=FALSE)
cuts = cut(lnhphist$breaks, c(-100,roundedbins,100))

binheight = lnhphist$counts[match(as.integer(roundedbins*5),as.integer(bins*5) )]

bandcolors = c("#62698d", "#5f6eb9", "#5d72d9", "#838fdb", "#b4b4dc", "#dcb4b4", "#dd7f7f", "#dd5d5d", "#de3737", "#df0000")
cutcolors = bandcolors[match(cuts,levels(cuts))]
# graph color learned from https://stackoverflow.com/questions/21858394/partially-color-histogram-in-r

pdf(file="~/git/heteropecilly/heteropecilly-v2/PAUL-90x341543v-C20.lnPIP.pdf", width=8, height=7)
plot(lnhphist, col=cutcolors, lwd=0.1, xlim=c(0,30), xlab="ln(PIP)", ylab="Number of sites", main="Histogram of heteropecilly scores (341k sites)", cex.axis=1.1, cex.lab=1.3)
#lines( rep(roundedbins,each=2), rep(c(0,5000),9), lwd=2)
segments( roundedbins+0.2, rep(0,9), x1=roundedbins+0.2, y1=binheight, lwd=4)
text(20,10000,"60089 constant sites\nwere not counted",col="#36a800", cex=1.5)
text(20,2000,"heteropecillious",col="#df0000", cex=1.5, pos=4)
text(4,13000,"homopecillious",col="#62698d", cex=1.5, pos=4)
dev.off()

png(file="~/git/heteropecilly/heteropecilly-v2/PAUL-90x341543v-C20.lnPIP.png", width=700, height=600)
#par(mar=c(4.5,4.5,4,2))
plot(lnhphist, col=cutcolors, lwd=0.1, xlim=c(0,30), xlab="ln(PIP)", ylab="Number of sites", main="Histogram of heteropecilly scores (341k sites)", cex.axis=1.2, cex.lab=1.5, cex.main=1.4)
segments( roundedbins+0.2, rep(0,9), x1=roundedbins+0.2, y1=binheight, lwd=4)
text(20,10000,"60089 constant sites\nwere not counted",col="#36a800", cex=1.8)
text(20,2000,"heteropecillious",col="#df0000", cex=1.8, pos=4)
text(4,13000,"homopecillious",col="#62698d", cex=1.8, pos=4)
dev.off()

#C,56642
#c,3447
#1,34154
#0,34157
#3,34154
#2,34154
#5,34154
#4,34154
#7,34154
#6,34154
#9,34154
#8,34154

sitelnldata = read.table("~/git/heteropecilly/RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab", header=TRUE, sep="\t")

t1nln = sitelnldata[,2]
t2nln = sitelnldata[,3]
t3nln = sitelnldata[,4]

#rowmins = apply( sitelnldata[,2:4], 1, which.min)
rowmaxs = apply( sitelnldata[,2:4], 1, which.max)

dSLS = ( (t1nln - t2nln) + (t1nln - t3nln) + (t2nln - t3nln) ) / 3
absdsls = abs(dSLS)

sum(dSLS[rowmins==1])
sum(dSLS[rowmins==2])
sum(dSLS[rowmins==3])

hptable = read.table("~/git/heteropecilly/heteropecilly-v2/hp_by_site_w_const_20CAT.tab", header=FALSE, sep="\t", skip=1)
lnhp = hptable[,3]
length(lnhp)
haslnhp = c(!is.na(lnhp))
length(haslnhp)
sortedlnhp = sort( as.numeric(lnhp[haslnhp]), index.return=TRUE )
length(sortedlnhp$x)

#rowmins_w_hp = rowmins[haslnhp]
rowmaxs_w_hp = rowmaxs[haslnhp]

sites_w_hp = sitelnldata[,1][haslnhp]

dsls_w_hp = dSLS[haslnhp]
absdsls_w_hp = absdsls[haslnhp]

sum(dsls_w_hp[rowmaxs_w_hp==1])
sum(dsls_w_hp[rowmaxs_w_hp==2])
sum(dsls_w_hp[rowmaxs_w_hp==3])

dsls_hp_sorted = dsls_w_hp[sortedlnhp$ix]
absdsls_hp_sorted = absdsls_w_hp[sortedlnhp$ix]

t1vt2 = t1nln - t2nln
t1vt2_w_hp = t1vt2[haslnhp]
t1vt2_hp_sorted = t1vt2_w_hp[sortedlnhp$ix]
abst1vt2 = abs(t1vt2_hp_sorted)

constantsites = is.na(lnhp)
length(constantsites)
absdsls_w_const = absdsls[constantsites]
length(absdsls_w_const)
abst1vt2_const = 

sum(t1nln[constantsites])
sum(t2nln[constantsites])
sum(t3nln[constantsites])

lnhphist = hist(lnhp, breaks=bins, plot=FALSE)
decilebreaks = sortedlnhp$x[floor(1:9/10*length(sortedlnhp$x))]
roundedbins = round(decilebreaks*5)/5
cuts = cut(sortedlnhp$x, c(-100,roundedbins,100))

lastdecile = sortedlnhp$ix[floor(length(sortedlnhp$ix)*0.9):length(sortedlnhp$ix)]
sum(dsls_w_hp[lastdecile][rowmins_w_hp[lastdecile]==1])
sum(dsls_w_hp[lastdecile][rowmins_w_hp[lastdecile]==2])
sum(dsls_w_hp[lastdecile][rowmins_w_hp[lastdecile]==3])

fadedcolors = c("#62698d33", "#5f6eb933", "#5d72d933", "#838fdb33", "#b4b4dc33", "#dcb4b433", "#dd7f7f33", "#dd5d5d33", "#de373733", "#df000033")

#precisecuts = cut(sortedlnhp$x[haslnhp],10)
pcutcolors = fadedcolors[match(cuts,levels(cuts))]

#dslssplits = split(dsls_w_hp[sortedlnhp$ix], ceiling(seq_along(dsls_w_hp[sortedlnhp$ix])/ ceiling(length(sortedlnhp$ix)/10) ))
#absdslssplits = split(absdsls_w_hp[sortedlnhp$ix], ceiling(seq_along(absdsls_w_hp[sortedlnhp$ix])/ ceiling(length(sortedlnhp$ix)/10) ))
#rowminsplits =  split(absdsls[sortedlnhp$ix], ceiling(seq_along(absdsls[sortedlnhp$ix])/ ceiling(length(sortedlnhp$ix)/10) ))
t1vt2_splits = split(t1vt2_hp_sorted, ceiling(seq_along(t1vt2_hp_sorted)/ ceiling(length(sortedlnhp$ix)/10) ))
abst1vt2_splits = split(abst1vt2, ceiling(seq_along(t1vt2_hp_sorted)/ ceiling(length(sortedlnhp$ix)/10) ))


#lapply(absdslssplits, sum)

lapply(t1vt2_splits, sum)
dt1vt2sums = lapply(abst1vt2_splits, sum)

pdf(file="~/git/heteropecilly/simion2017_hp_vs_dln_boxplot.pdf", width=8, height=7)
par(mar=c(4.5,4.5,1,1))
boxplot(abst1vt2_splits, col=bandcolors, pch=16, outcol=fadedcolors, ylim=c(-0.5,8), ylab=expression(paste("|",Delta,"ln(L)|",sep="")), xlab="Heteropecilly interval", cex.lab=1.4, cex.axis=1.4 )
text(c(1:10), rep(c(-0.5),10), round(unlist(dt1vt2sums),0), cex=1.1)
text(0.4,-0.5, expression(paste(Sigma,"=")), cex=1.2)
dev.off()

pdf(file="~/git/heteropecilly/simion2017_hp_vs_absdln_plot.pdf", width=9, height=8)
par(mar=c(4.5,4.5,1,1))
plot(sortedlnhp$x, abst1vt2, type='p', col=pcutcolors, pch=16, xlim=c(0,45), ylim=c(0,8), ylab=expression(paste("|",Delta,"ln(L)|",sep="")), xlab="lnPIP", cex.lab=1.4, cex.axis=1.4 )
hplnl_lm = lm(abst1vt2 ~ sortedlnhp$x )
lines(hplnl_lm)
dev.off()

png(file="~/git/heteropecilly/simion2017_hp_vs_dln_boxplot.png", width=600, height=500)
par(mar=c(4.5,4.5,1,1))
boxplot(abst1vt2_splits, col=bandcolors, pch=16, outcol=fadedcolors, ylim=c(-0.5,8), ylab=expression(paste("|",Delta,"ln(L)|",sep="")), xlab="Heteropecilly interval", cex.lab=1.4, cex.axis=1.4 )
text(c(1:10), rep(c(-0.5),10), round(unlist(dt1vt2sums),0), cex=1.1)
text(0.4,-0.5, expression(paste(Sigma,"=")), cex=1.2)
dev.off()

png(file="~/git/heteropecilly/simion2017_hp_vs_absdln_plot.png", width=800, height=700)
par(mar=c(4.5,4.5,1,1))
plot(sortedlnhp$x, abst1vt2, type='p', col=pcutcolors, pch=16, xlim=c(0,45), ylim=c(0,8), ylab=expression(paste("|",Delta,"ln(L)|",sep="")), xlab="lnPIP", cex.lab=1.4, cex.axis=1.4 )
#lines(hplnl_lm)
dev.off()



dslssums = lapply( dslssplits , function(x) sum(x[x>0.5]) )
absdslssums = lapply( absdslssplits, function(x) sum(x[x>0.5]))
unlist(absdslssums) - unlist(dslssums)

pdf(file="~/git/heteropecilly/simion2017_hp_vs_dsls.pdf", width=11, height=10)
par(mar=c(4.5,4.5,2,2))
plot(sortedlnhp$x, dsls_w_hp[sortedlnhp$ix], xlim=c(0,45), ylim=c(-4,4), type='p', pch=16, col=pcutcolors, xlab="lnPIP", ylab="dSLS", cex.axis=1.4, cex.lab=1.4)
points(rep(-1,length(dSLS[constantsites])), dSLS[constantsites], col="#36a80033",pch=16)
abline(h=c(0.5,-0.5), lwd=1, lty=2)
dev.off()

png(file="~/git/heteropecilly/simion2017_hp_vs_dsls.png", width=800, height=700)
par(mar=c(4.5,4.5,2,2))
plot(sortedlnhp$x, dsls_w_hp[sortedlnhp$ix], xlim=c(0,45), ylim=c(-4,5), type='p', pch=16, col=pcutcolors, xlab="lnPIP", ylab="dSLS", cex.axis=1.4, cex.lab=1.4)
points(rep(-1,length(dSLS[constantsites])), dSLS[constantsites], col="#36a80033",pch=16)
abline(h=c(0.5,-0.5), lwd=1, lty=2)
dev.off()




#
