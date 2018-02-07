# plot chi-squared test vs gaps for amino acid counts
# created by wrf 2018-01-29

args = commandArgs(trailingOnly=TRUE)

#inputfile = "~/git/heteropecilly/aa_counts/philippe2009_aa_counts.tab"
inputfile = args[1]

outputfile = gsub("([\\w/]+)\\....","\\1.chisq.pdf",inputfile,perl=TRUE)
aadata = read.table(inputfile, header=TRUE, sep='\t', row.names=c(1) )

numtaxa = match("Constants",row.names(aadata))-1

aadatastrict = aadata[1:numtaxa,1:20]
effectivesites = rowSums(aadatastrict)
aascales = scale( aadatastrict / effectivesites * 1000 )
zmeans = rowMeans(aascales)

aachisq = chisq.test(aadata[1:numtaxa,1:20])
chisqsums = rowSums( (aachisq$observed - aachisq$expected)^2 / aachisq[["expected"]] )

df = as.integer(aachisq$parameter)
chisqsig = qchisq(0.99, df)

numsites = sum(aadata[1,])
#xmax = max(pretty(numsites))
ymax = max(pretty(max(chisqsums)))

xlabel = paste("Gaps % (out of",numsites,"sites)")

colorindex = as.integer(chisqsums>=chisqsig)+1
colorset = c("#0868ac", "#e34a33")

gappercent = aadata[1:(dim(aadata)[1]-4),21]/numsites*100

pdf(file=outputfile, width=12, height=7)
par(mar=c(4.5,4.5,4,2), mfrow=c(1,2))
#plot( aadata[1:(dim(aadata)[1]-4),21], log(chisqsums), frame.plot=FALSE, xlim=c(0,xmax), ylim=c(2.5,7.5), log="y", main=inputfile, xlab=xlabel, ylab="Chi-square of AA freq", pch=16, col="#0868ac", cex.lab=1.2)
plot( gappercent, chisqsums, frame.plot=FALSE, xlim=c(0,100), ylim=c(0,ymax), main=inputfile, xlab=xlabel, ylab="Chi-square of AA freq", pch=16, col=colorset[colorindex], cex.lab=1.2)
speciesnames = sapply(strsplit(row.names(aadata)[1:numtaxa],"_"), function(x){x[1]})
text( gappercent+1, chisqsums, speciesnames, pos=4)
mtext( paste("n taxa =",numtaxa) , side=3, line=1, at=c(-5), cex=1.3 )
#mtext( expression(paste( chi^2, "=", )) , side=3, line=1, at=c(-5), cex=1.3 )
par(mar=c(4.5,4,4,1))
plot( gappercent, zmeans, xlim=c(0,100), pch=16, frame.plot=FALSE, xlab=xlabel, ylab="Average Z-score of AA freq", main=inputfile, col=colorset[colorindex], cex.lab=1.2)
text( gappercent, zmeans, speciesnames, pos=4)

dev.off()