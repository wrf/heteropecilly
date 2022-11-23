# compare amino acid frequency across datasets

inputfile = "~/git/heteropecilly/aa_counts/simion2017_aa_counts.tab"
datasetname = "Simion 2017"

simiondata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )

aadata = simiondata
allconst = aadata[dim(aadata)[1]-3,]
print(allconst)
alfullhumanprots = aadata[dim(aadata)[1]-1,]
allhuman = aadata[dim(aadata)[1],]

normsites = rowSums(aadata[,1:20]) # excludes gaps and X

normcounts = aadata[,1:20]/normsites

mostcommonaas = colSums(aadata[1:(dim(aadata)[1]-2),]) # sums count constants, but not additional files
simsortedmostcommon = sort(mostcommonaas[1:20], decreasing=TRUE,  index.return = TRUE)

pdf(file="~/git/heteropecilly/aa_counts/aa_counts_from_5_sets.pdf", width=7, height=6.5)
par(mar=c(4,5,3.5,1))


plot(c(1,20),c(0,1),xlim=c(1,20),ylim=c(0,0.15),type='n', axes=FALSE, ylab="Fraction of sites", xlab="Amino acid", main=paste("AA frequency from metazoan supermatrices"), cex.lab=1.4)
axis(2, cex.axis=1.3)
axis(1, at=1:20, labels=names(aadata)[sortedmostcommon$ix] )



aapositions = 1:20
aamaxes = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,max)[simsortedmostcommon$ix]
aamins = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,min)[simsortedmostcommon$ix]

polygon(c(aapositions, rev(aapositions)), c(aamaxes,rev(aamins)), border="#2211cc", col="#2211cc33")
points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#2211cc", pch=17, cex=1.4, lwd=2)
#points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#2211cc", pch=15, cex=1.4)

inputfile = "~/git/heteropecilly/aa_counts/cannon2016_aa_counts.tab"
aadata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )
allconst = aadata[dim(aadata)[1]-3,]
print(allconst)
alfullhumanprots = aadata[dim(aadata)[1]-1,]
normsites = rowSums(aadata[,1:20]) # excludes gaps and X
normcounts = aadata[,1:20]/normsites
aamaxes = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,max)[simsortedmostcommon$ix]
aamins = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,min)[simsortedmostcommon$ix]
polygon(c(aapositions, rev(aapositions)), c(aamaxes,rev(aamins)), border="#cccc99", col="#cccc9933")
points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#cccc99", pch=17, cex=1.4, lwd=2)
#points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#cccc99", pch=15, cex=1.4)

inputfile = "~/git/heteropecilly/aa_counts/whelan2017_aa_counts.tab"
aadata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )
allconst = aadata[dim(aadata)[1]-3,]
print(allconst)
alfullhumanprots = aadata[dim(aadata)[1]-1,]
normsites = rowSums(aadata[,1:20]) # excludes gaps and X
normcounts = aadata[,1:20]/normsites
aamaxes = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,max)[simsortedmostcommon$ix]
aamins = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,min)[simsortedmostcommon$ix]
polygon(c(aapositions, rev(aapositions)), c(aamaxes,rev(aamins)), border="#d95f02", col="#d95f0233")
points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#d95f02", pch=17, cex=1.4, lwd=2)
#points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#d95f02", pch=15, cex=1.4)

inputfile = "~/git/heteropecilly/aa_counts/borowiec2015_aa_counts.tab"
aadata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )
allconst = aadata[dim(aadata)[1]-3,]
print(allconst)
alfullhumanprots = aadata[dim(aadata)[1]-1,]
normsites = rowSums(aadata[,1:20]) # excludes gaps and X
normcounts = aadata[,1:20]/normsites
aamaxes = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,max)[simsortedmostcommon$ix]
aamins = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,min)[simsortedmostcommon$ix]
polygon(c(aapositions, rev(aapositions)), c(aamaxes,rev(aamins)), border="#1b9e77", col="#1b9e7733")
points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#1b9e77", pch=17, cex=1.4, lwd=2)
#points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#1b9e77", pch=15, cex=1.4)

inputfile = "~/git/heteropecilly/aa_counts/philippe2009_aa_counts.tab"
aadata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )
allconst = aadata[dim(aadata)[1]-3,]
print(allconst)
alfullhumanprots = aadata[dim(aadata)[1]-1,]
normsites = rowSums(aadata[,1:20]) # excludes gaps and X
normcounts = aadata[,1:20]/normsites
aamaxes = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,max)[simsortedmostcommon$ix]
aamins = apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,min)[simsortedmostcommon$ix]
polygon(c(aapositions, rev(aapositions)), c(aamaxes,rev(aamins)), border="#e7298a", col="#e7298a33")
points(1:20, (allconst[1:20]/sum(allconst[1:20]))[sortedmostcommon$ix], col="#e7298a", pch=17, cex=1.4, lwd=2)
#points(1:20, (alfullhumanprots[1:20]/sum(alfullhumanprots[1:20]))[sortedmostcommon$ix], col="#e7298a", pch=15, cex=1.4)
constvsnorm = sort( allconst[1:20]/sum(allconst[1:20]) /  apply(normcounts[1:(dim(aadata)[1]-4),1:20],2,median) , decreasing=TRUE)
print(constvsnorm)

#points(1:20, (allhuman[1:20]/sum(allhuman[1:20]))[sortedmostcommon$ix], col="#000000", pch=17, cex=1.4)

legend(12.9,0.155,legend=c("Whelan 2017","Cannon 2016", "Simion 2017", "Borowiec 2015", "Philippe 2009"), pch=17, lwd=2, col=c("#d95f02","#cccc99", "#2211cc", "#1b9e77", "#e7298a" ), pt.bg=c("#d95f0233","#cccc9933", "#2211cc33", "#1b9e7733", "#e7298a33" ), cex=1.3 )
#legend(14,0.104,legend=c("Constant sites","Human prots"), pch=c(23,17), pt.lwd=2, col=c("#d95f02","#000000"), cex=1.3 )


#legendcolors = c("#000000", "#36a800", "#ea5d1f" , "#62698d" )
#legend(10,0.15, legend=c("All species",paste("Constant sites (",sum(allconst),")",sep=''),"Full-length human prots","All human proteins"), pch=c(16,15,15,17), col=legendcolors, cex=1.3 )

dev.off()
#