# plot stats from alignment pairs of alignment protein to reference protein

library(VennDiagram)
library(gridExtra)

simiondata = read.table("~/git/heteropecilly/pair_stats/simion2017_align_pair_stats.tab", header=TRUE, sep="\t")

spanpct = simiondata[,7]
ungappedpct = simiondata[,4]
refprotlength = simiondata[,8]

# simiondata[refprotlength>=2000,]

simionsitefraction = sum(simiondata[,3])/sum(simiondata[,8])
simionmeancompleteness = mean(simiondata[,3]/simiondata[,8])

whelandata = read.table("~/git/heteropecilly/pair_stats/whelan2017_align_pair_stats.tab", header=TRUE, sep="\t")

whelansitefraction = sum(whelandata[,3])/sum(whelandata[,8])
whelanmeancompleteness = mean(whelandata[,3]/whelandata[,8])
# whelandata[whelandata[,8]>=1000,]

borowiecdata = read.table("~/git/heteropecilly/pair_stats/borowiec2015_align_pair_stats.tab", header=TRUE, sep="\t")

borowiecsitefraction = sum(borowiecdata[,3])/sum(borowiecdata[,8])

philippedata = read.table("~/git/heteropecilly/pair_stats/philippe2009_align_pair_stats.tab", header=TRUE, sep="\t")

philippesitefraction = sum(philippedata[,3])/sum(philippedata[,8])
philippemeancompleteness = mean(philippedata[,3]/philippedata[,8])


w_in_s_prots = match(whelandata[,2], simiondata[,2])
w_in_s_prots = w_in_s_prots[!is.na(w_in_s_prots)]
sharedprots = refprotlength[w_in_s_prots]

w_in_b_prots = match(whelandata[,2], borowiecdata[,2])
w_in_b_prots = w_in_b_prots[!is.na(w_in_b_prots)]

s_in_b_prots = match(simiondata[,2], borowiecdata[,2])
s_in_b_prots = s_in_b_prots[!is.na(s_in_b_prots)]

p_in_s_prots = match(philippedata[,2], simiondata[,2])
p_in_s_prots = p_in_s_prots[!is.na(p_in_s_prots)]

w_in_s_in_b = match(simiondata[,2][w_in_s_prots], borowiecdata[,2])
w_in_s_in_b = w_in_s_in_b[!is.na(w_in_s_in_b)]

erwinprots = c("sp|P05062|ALDOB_HUMAN","sp|P31153|METK2_HUMAN","sp|P06576|ATPB_HUMAN","sp|P60174|TPIS_HUMAN","sp|P68104|EF1A1_HUMAN","sp|P08237|PFKAM_HUMAN","sp|P08237|PFKAM_HUMAN","sp|P04040|CATA_HUMAN")

pdf(file="~/git/heteropecilly/pair_stats/align_pair_stats_simion-v-whelan.pdf", width=8, height=7)
#png(file="~/git/heteropecilly/pair_stats/align_pair_stats_simion-v-whelan.png", width=800, height=700)
plot(refprotlength, ungappedpct, col="#2211cc33", pch=16, main="Fraction of each human protein kept in final alignments", xlab="Human reference protein length (AAs)", ylab="Fraction coverage in alignment", cex.axis=1.3, cex.lab=1.3)
refprothist = hist(refprotlength,breaks=seq(0,5000,25),plot=FALSE)
#lines(refprothist$breaks[2:length(refprothist$breaks)], refprothist$counts/max(refprothist$counts)/2, lwd=3 )

points(whelandata[,8], whelandata[,4], col="#d95f0233", pch=16)
points(borowiecdata[,8], borowiecdata[,4], col="#1b9e7733", pch=16)
points(philippedata[,8], philippedata[,4], col="#e7298a33", pch=16)

legend(3000,1.01,legend=c("Whelan 2017","Simion 2017", "Borowiec 2015", "Philippe 2009"), pch=15, col=c("#d95f02","#2211cc", "#1b9e77", "#e7298a" ), cex=1.3 )

text(3000,0.70,"212 of 212 prots", col=c("#d95f02"), cex=1.3, pos=4)
text(3500,0.65, paste(round(whelansitefraction,digits=3),"of reference sites"), col=c("#d95f02"), cex=1.3)

text(3000,0.58,"1499 of 1719 prots", col=c("#2211cc"), cex=1.3, pos=4)
text(3500,0.53, paste(round(simionsitefraction,digits=3),"of reference sites"), col=c("#2211cc"), cex=1.3)

text(3000,0.46,"1056 of 1080 prots", col=c("#1b9e77"), cex=1.3, pos=4)
text(3500,0.41, paste(round(borowiecsitefraction,digits=3),"of reference sites"), col=c("#1b9e77"), cex=1.3)

text(3000,0.34,"128 of 128 prots", col=c("#e7298a"), cex=1.3, pos=4)
text(3500,0.29, paste(round(philippesitefraction,digits=3),"of reference sites"), col=c("#e7298a"), cex=1.3)

#text(3600,0.40,paste(length(sharedprots),"proteins\ncommon to both sets"), col=c("#881188"), cex=1.3)
#par( fig = c(grconvertX(c(3000,4200), from="user", to="ndc"), grconvertY(c(0.1,0.4), from="user", to="ndc")), new=TRUE, mar=c(0,0,0,0))
#subvp = viewport(x=0.57, y=0.1, width=0.33, height=0.3)
#grid.arrange( gTree(triplevenn), gTree(triplevenn), gTree( triplevenn ) , ncol=3)
#grid.draw(triplevenn)
dev.off()

png(file="~/git/heteropecilly/pair_stats/align_pair_stats_simion-v-whelan.png", width=600, height=500)
plot(refprotlength, ungappedpct, col="#2211cc66", pch=16, main="Fraction of each human protein kept in final alignments", xlab="Human reference protein length (AAs)", ylab="Fraction coverage in alignment", cex.axis=1.3, cex.lab=1.3)
refprothist = hist(refprotlength,breaks=seq(0,5000,25),plot=FALSE)
#lines(refprothist$breaks[2:length(refprothist$breaks)], refprothist$counts/max(refprothist$counts)/2, lwd=3 )

points(whelandata[,8], whelandata[,4], col="#d95f0266", pch=16)
points(borowiecdata[,8], borowiecdata[,4], col="#1b9e7766", pch=16)

legend(3000,1.0,legend=c("Whelan 2017","Simion 2017", "Borowiec 2015"), pch=15, col=c("#d95f02","#2211cc", "#1b9e77" ), cex=1.3 )

text(3000,0.70,"212 of 212 prots", col=c("#d95f02"), cex=1.3, pos=4)
text(3500,0.65, paste(round(whelansitefraction,digits=3),"of reference sites"), col=c("#d95f02"), cex=1.3)

text(3000,0.58,"1499 of 1719 prots", col=c("#2211cc"), cex=1.3, pos=4)
text(3500,0.53, paste(round(simionsitefraction,digits=3),"of reference sites"), col=c("#2211cc"), cex=1.3)

text(3000,0.46,"1056 of 1080 prots", col=c("#1b9e77"), cex=1.3, pos=4)
text(3500,0.41, paste(round(borowiecsitefraction,digits=3),"of reference sites"), col=c("#1b9e77"), cex=1.3)
dev.off()




pdf(file="~/git/heteropecilly/pair_stats/align_pairs_3-way_venn.pdf", width=7, height=7)
triplevenn = draw.triple.venn( area1=length(whelandata[,8]) , area2=length(simiondata[,8]), area3=length(borowiecdata[,8]), n12=length(w_in_s_prots), n23=length(s_in_b_prots) , n13=length(w_in_b_prots), n123=length(w_in_s_in_b) ,  fill=c( "#cc181833" , "#2211cc33" , "#18cc1833"), category=c("Whelan 2017","Simion 2017","Borowiec 2015") , cex=2, cat.cex=2, cat.pos=c(-20,20,180) )
dev.off()

png(file="~/git/heteropecilly/pair_stats/align_pairs_3-way_venn.png", width=450, height=450)
triplevenn = draw.triple.venn( area1=length(whelandata[,8]) , area2=length(simiondata[,8]), area3=length(borowiecdata[,8]), n12=length(w_in_s_prots), n23=length(s_in_b_prots) , n13=length(w_in_b_prots), n123=length(w_in_s_in_b) ,  fill=c( "#cc181833" , "#2211cc33" , "#18cc1833"), category=c("Whelan 2017","Simion 2017","Borowiec 2015") , cex=2, cat.cex=2, cat.pos=c(-20,20,180) )
dev.off()










#
