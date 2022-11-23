# plot frequency of constant sites relative to sites overall in an alignment
# using data from simion 2017

inputfile = "~/git/heteropecilly/aa_counts/simion2017_aa_counts.tab"
datasetname = "Simion 2017"
inputfile = "~/git/heteropecilly/aa_counts/cannon2016_aa_counts.tab"
datasetname = "Cannon 2016"
inputfile = "~/git/heteropecilly/aa_counts/whelan2017_aa_counts.tab"
datasetname = "Whelan 2017"
#inputfile = "~/git/heteropecilly/aa_counts/borowiec2015_aa_counts.tab"
#datasetname = "Borowiec 2015"
inputfile = "~/git/heteropecilly/aa_counts/philippe2009_aa_counts.tab"
datasetname = "Philippe 2009"


aadata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )

holozoaconst = aadata[dim(aadata)[1]-2,]
allconst = aadata[dim(aadata)[1]-3,]
alfullhumanprots = aadata[dim(aadata)[1]-1,]
allhuman = aadata[dim(aadata)[1],]

normsites = rowSums(aadata[,1:20]) # excludes gaps and X

normcounts = aadata[,1:20]/normsites

mostcommonaas = colSums(aadata[1:(dim(aadata)[1]-2),]) # sums count constants, but not additional files
sortedmostcommon = sort(mostcommonaas[1:20], decreasing=TRUE,  index.return = TRUE)

overallfreq = colSums(aadata[1:(dim(aadata)[1]-4),1:20]) / sum(aadata[1:(dim(aadata)[1]-4),1:20])

for (s in 1:(dim(aadata)[1]-4)) {
	simion_ts_sum = 0
	for (i in 1:20) {
		ch_sq_val = (normcounts[s,i] - overallfreq[i])^2 / overallfreq[i]
		simion_ts_sum = simion_ts_sum + ch_sq_val
	}
	print(simion_ts_sum * 20 )
}
aachisq = chisq.test(aadata[1:(dim(aadata)[1]-4),1:20])

chisqsums = rowSums( (aachisq$observed - aachisq$expected)^2 / aachisq[["expected"]] )
plot( aadata[1:(dim(aadata)[1]-4),21], chisqsums, frame.plot=FALSE, xlim=c(0,30000), xlab="Gaps", ylab="Chi-square of AA freq", pch=16, col="#0868ac")
speciesnames = sapply(strsplit(row.names(aadata)[1:(dim(aadata)[1]-4)],"_"), function(x){x[1]})
text( aadata[1:(dim(aadata)[1]-4),21]+10, chisqsums, speciesnames, pos=4)


#
leucinenc = normcounts[1:(dim(aadata)[1]-4),10]
library(nortest)
ad.test(leucinenc)

qqnorm(leucinenc)
qqline(leucinenc, col="red")

sortnormcounts = normcounts[1:(dim(aadata)[1]-4),sortedmostcommon$ix]
par(mfrow=c(4,5))
for (i in 1:20) {
qplot = qqnorm(sortnormcounts[,i])
ymax = max(qplot$y)
ymin = min(qplot$y)
text(-2,ymax-(ymax-ymin)*0.1,names(sortedmostcommon$x)[i], cex=1.5)
text(-2,ymax-(ymax-ymin)*0.25,round(ad.test(sortnormcounts[,i])$p.value,digits=3), cex=1.2)
qqline(sortnormcounts[,i], col='red')
}


par(mfrow=c(1,1))
outputfile = gsub("([\\w/]+)\\....","\\1.pdf",inputfile,perl=TRUE)
pdf(file=outputfile, width=7, height=6.5)
par(mar=c(4,5,3.5,1))

plot(c(1,20),c(0,1),xlim=c(1,20),ylim=c(0,0.15),type='n', axes=FALSE, ylab="Fraction of sites", xlab="Amino acid", main=paste("AA frequency from",datasetname), cex.lab=1.4)
axis(2, cex.axis=1.3)
axis(1, at=1:20, labels=names(aadata)[sortedmostcommon$ix] )
# for loop to only take first 97 entries, meaning all normal taxa
for (i in 1:(dim(aadata)[1]-4)) {
points(1:20, (normcounts[i,1:20])[sortedmostcommon$ix], cex=1.4, pch=16, col="#00000044" )
}
constaadists = allconst[1:20]/mostcommonaas[1:20]*10
constnormby005 = constaadists / 0.05

points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#36a800", pch=15, cex=1.4)
#points(1:20, (holozoaconst[1:20]/sum(holozoaconst[1:20]))[sortedmostcommon$ix], col="#71c872", pch=17, cex=1.3)
#lines(1:20, constnormby005[sortedmostcommon$ix]/20, col="blue", lwd=3)
points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#ea5d1f", pch=15, cex=1.4)
points(1:20, (allhuman[1:20]/sum(allhuman[1:20]))[sortedmostcommon$ix], col="#62698d", pch=17, cex=1.4)

legendcolors = c("#000000", "#36a800", "#ea5d1f" , "#62698d" )
legend(10,0.15, legend=c("All species",paste("Constant sites (",sum(allconst),")",sep=''),"Full-length human prots","All human proteins"), pch=c(16,15,15,17), col=legendcolors, cex=1.3 )

dev.off()
#

png(file=gsub("([\\w/]+)\\....","\\1.png",inputfile,perl=TRUE), width=700, height=600)
par(mar=c(4,5,4,1))
plot(c(1,20),c(0,1),xlim=c(1,20),ylim=c(0,0.15),type='n', axes=FALSE, ylab="Fraction of sites", xlab="Amino acid", main=paste("AA frequency from",datasetname), cex.lab=1.7, cex.main=1.3)
axis(2, cex.axis=1.5)
axis(1, at=1:20, labels=names(aadata)[sortedmostcommon$ix], cex.axis=1.4 )
for (i in 1:97) {
points(1:20, (normcounts[i,1:20])[sortedmostcommon$ix], cex=1.6, pch=16, col="#00000044" )
}
constaadists = allconst[1:20]/mostcommonaas[1:20]*10
constnormby005 = constaadists / 0.05

points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#36a800", pch=15, cex=1.6)
#points(1:20, (holozoaconst[1:20]/sum(holozoaconst[1:20]))[sortedmostcommon$ix], col="#71c872", pch=17, cex=1.3)
#lines(1:20, constnormby005[sortedmostcommon$ix]/20, col="blue", lwd=3)
points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#ea5d1f", pch=15, cex=1.6)
points(1:20, (allhuman[1:20]/sum(allhuman[1:20]))[sortedmostcommon$ix], col="#62698d", pch=17, cex=1.6)

legendcolors = c("#000000", "#36a800", "#ea5d1f" , "#62698d" )
legend(10,0.15, legend=c("All species",paste("Constant sites (",sum(allconst),")",sep=''),"Full-length human prots","All human proteins"), pch=c(16,15,15,17), col=legendcolors, cex=1.6 )

dev.off()