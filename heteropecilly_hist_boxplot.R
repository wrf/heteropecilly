infile = "~/git/heteropecilly/PAUL-90x341543v-C20.proba.tab"

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



# remove constant sites, need to identify which is which
hptable = read.table("~/git/heteropecilly/simion2017_20CAT_hp_by_site_w_const.tab", header=FALSE, sep="\t", skip=1)
lnhp = hptable[,3]
haslnhp = c(!is.na(lnhp))
sortedlnhp = sort( as.numeric(lnhp[haslnhp]), index.return=TRUE )

sitelnldata = read.table("~/project/phylogeny/RAxML_perSiteLLs.simion2017_97sp_401632pos_1719genes.tab", header=TRUE, sep="\t")
t1nln = sitelnldata[,2]
t2nln = sitelnldata[,3]
t3nln = sitelnldata[,4]
sites_w_hp = sitelnldata[,1][haslnhp]

dSLS = ( (t1nln - t2nln) + (t1nln - t3nln) + (t2nln - t3nln) ) / 3
absdsls = abs(dSLS)
dsls_w_hp = dSLS[haslnhp]
absdsls_w_hp = absdsls[haslnhp]

dsls_hp_sorted = dsls_w_hp[sortedlnhp$ix]
absdsls_hp_sorted = absdsls_w_hp[sortedlnhp$ix]

t1vt2 = t1nln - t2nln
t1vt2_w_hp = t1vt2[haslnhp]
t1vt2_hp_sorted = t1vt2_w_hp[sortedlnhp$ix]
abst1vt2 = abs(t1vt2_hp_sorted)

t1vt2_splits = split(t1vt2_hp_sorted, ceiling(seq_along(t1vt2_hp_sorted)/ ceiling(length(sortedlnhp$ix)/10) ))
abst1vt2_splits = split(abst1vt2, ceiling(seq_along(t1vt2_hp_sorted)/ ceiling(length(sortedlnhp$ix)/10) ))


#lapply(absdslssplits, sum)

lapply(t1vt2_splits, sum)
dt1vt2sums = lapply(abst1vt2_splits, sum)

fadedcolors = c("#62698d33", "#5f6eb933", "#5d72d933", "#838fdb33", "#b4b4dc33", "#dcb4b433", "#dd7f7f33", "#dd5d5d33", "#de373733", "#df000033")

mainlab = "Histogram of heteropecilly scores (341k sites)"

pdf(file="~/git/heteropecilly/heteropecilly_hist_boxplot.pdf", width=10, height=5)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,2,1))
plot(lnhphist, col=cutcolors, lwd=0.1, xlim=c(0,30), xlab="ln(PIP)", ylab="Number of sites", main="", cex.axis=1.2, cex.lab=1.3)
#lines( rep(roundedbins,each=2), rep(c(0,5000),9), lwd=2)
segments( roundedbins+0.2, rep(0,9), x1=roundedbins+0.2, y1=binheight, lwd=2.5)
text(20,10000,"60089 constant sites\nwere not counted",col="#36a800", cex=1.2)
text(17,2000,"heteropecillious",col="#df0000", cex=1.3, pos=4)
text(4,13000,"homopecillious",col="#62698d", cex=1.3, pos=4)

par(mar=c(4.5,2,2,1))
boxplot(abst1vt2_splits, col=bandcolors, pch=16, outcol=fadedcolors, ylim=c(-0.5,6), ylab="", xlab="Heteropecilly decile", cex.lab=1.3, cex.axis=1.2 )
mtext(expression(paste("|",Delta,"ln(L)|",sep="")), side=2, cex=1.3, line=2)
text(c(1:10), rep(c(-0.5),10), round(unlist(dt1vt2sums),0), cex=0.9)
text(0.4,-0.5, expression(paste(Sigma,"=")), cex=1)
dev.off()
