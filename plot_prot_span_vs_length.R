#!/usr/bin/env Rscript
# generate figure of trimmed length versus alignment span
# v1 created 2017-11-20

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/git/heteropecilly/cannon2016_align_pair_stats.tab"
outputfile = gsub("([\\w/]+)\\....","\\1.pdf",inputfile,perl=TRUE)

pairdat = read.table(inputfile, header=TRUE, sep="\t")

trimspandiff = pairdat[["spanPercent"]] - pairdat[["trimmedPercent"]]

colorset = c( '#b10026aa','#e31a1caa','#fc4e2aaa', '#fd8d3caa','#feb24caa','#fed976aa', '#ffeda0aa', '#edf8b1aa','#c7e9b4aa','#7fcdbbaa', '#41b6c4aa','#1d91c0aa','#225ea8aa', '#0c2c84aa') 
alignlength = pairdat[["refProtLength"]]
#bins = unique( c(pretty(c(min(alignlength), median(alignlength)),7), pretty(c(median(alignlength),max(alignlength)),7) ) )
bins = pretty( range(alignlength), 14 )
lengthgroups = cut(alignlength,breaks=bins)

partlists = strsplit(as.character(pairdat[["partition"]]),"_")
partitions = unlist(lapply(partlists, function(x) x[3]))
shortparts = lapply(strsplit(partitions, "-"), function(x) x[1])

morethan35 = trimspandiff > 0.3

xpos = 1:length(trimspandiff)

pdf(outputfile, width=4+length(trimspandiff)%/%16, height=7)
par(mar=c(6,4.5,4,1))
plot(xpos, trimspandiff, type='n', xlab="", ylab="Percent difference between number of sites and span", main=inputfile, axes=FALSE, cex.lab=1.3 )
partxpos = seq(0,length(trimspandiff),10)
abline( v=partxpos, col="lightgray" )
points(xpos, trimspandiff, pch=21, bg=colorset[lengthgroups])
axis(2, cex=1.4)
axis(1, at=partxpos, labels=FALSE )
text(xpos[morethan35], trimspandiff[morethan35]+0.01, shortparts[morethan35], pos=4 )
mtext(partitions,side=1,at=c(xpos),las=2, cex=0.7, line=1)
dev.off()


#