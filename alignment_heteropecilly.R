#

hpdata = read.table("~/git/heteropecilly/hp_by_protein_w_const.tab", header=TRUE, sep="\t")

constants = hpdata[["C"]] + hpdata[["c"]]
gaps = hpdata[["gaps"]]
alignsites = hpdata[,3]

bandcolors = c("#62698d", "#5f6eb9", "#5d72d9", "#838fdb", "#b4b4dc", "#dcb4b4", "#dd7f7f", "#dd5d5d", "#de3737", "#df0000", "#36a800", "#646464")
#bandcolors = c("#62698d", "#5f6eb9", "#5d72d9", "#838fdb", "#b4b4dc", "#dcb4b4", "#dd7f7f", "#dd5d5d", "#de3737", "#df0000", "#71c872", "#36a800", "#646464")

hponly = hpdata[,6:15]
hpwithsites = cbind(hponly, constants, alignsites)
hpwithsites_df = data.frame(hpwithsites)

drawhplines <- function(x){
#expectedcount = x[1:11]/c(x[12])
expectedcount = (x[1:11] - rep(x[12]/10,11) ) / c(x[12])
#lines( 0:10, expectedcount)
points( 0:10, expectedcount, pch=20, col=bandcolors)
}
plot(0,0,type='n', ylim=c(-0.1,0.5), xlim=c(0,11))
segments(0,0.1,11)
segments(0,-0.1,11)
apply(hpwithsites_df, 1, drawhplines)


makehphist <- function(x,y){
bins = seq(-0.1,0.5,0.01)
#expectedcount = (x[1:11] - rep(x[12]/10,11) ) / c(x[12])

expectedcount = (x - y/10) / y
expectedcount[expectedcount>0.2] = 0.2
hph = hist(expectedcount, breaks=bins, plot=FALSE)
hphc = hph$counts
print(hphc)
return(hphc)
}

pdf(file="~/git/heteropecilly/hp_by_protein_w_const.pdf", width=8, height=7)
#hpexpected = apply(hpwithsites[1:11], 2, makehphist, alignsites)
hpexpected = apply(hponly, 2, makehphist, alignsites-constants)

plot(0,0,type='n', ylim=c(0,300), xlim=c(-0.1,0.2), xlab="Ehp", ylab="Number of sites")

for (i in 1:10) {
lines(seq(-0.1,0.2,0.01), hpexpected[,i][1:31], col=bandcolors[i], lwd=3)
}
dev.off()



#