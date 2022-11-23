# plot stats from alignment pairs of alignment protein to reference protein

library(VennDiagram)
library(gridExtra)

simiondata = read.table("~/git/heteropecilly/pair_stats/simion2017_align_pair_stats.tab", header=TRUE, sep="\t")

spanpct = simiondata[,7]
ungappedpct = simiondata[,4]
refprotlength = simiondata[,8]

simion_namesplit = sapply(as.character(simiondata[,2]), function(x) unlist(strsplit(x, split='|', fixed=TRUE))[3] )
simion_shortname = sapply(as.character(simion_namesplit), function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1] )
#simion_shortname
simionlongprots = (refprotlength - ungappedpct*1000) >= 1800


simionsitefraction = sum(simiondata[,3])/sum(simiondata[,8])
simionmeancompleteness = mean(simiondata[,3]/simiondata[,8])

whelandata = read.table("~/git/heteropecilly/pair_stats/whelan2017_align_pair_stats.tab", header=TRUE, sep="\t")
whelan_namesplit = sapply(as.character(whelandata[,2]), function(x) unlist(strsplit(x, split='|', fixed=TRUE))[3] )
whelan_shortname = sapply(as.character(whelan_namesplit), function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1] )
whelan_longprots = (whelandata[,8] - whelandata[,4]*2000) >= 1000
whelansitefraction = sum(whelandata[,3])/sum(whelandata[,8])
whelanmeancompleteness = mean(whelandata[,3]/whelandata[,8])
# whelandata[whelandata[,8]>=1000,]

borowiecdata = read.table("~/git/heteropecilly/pair_stats/borowiec2015_align_pair_stats.tab", header=TRUE, sep="\t")
borowiec_namesplit = sapply(as.character(borowiecdata[,2]), function(x) unlist(strsplit(x, split='|', fixed=TRUE))[3] )
borowiec_shortname = sapply(as.character(borowiec_namesplit), function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1] )
borowiec_longprots = (borowiecdata[,8] - borowiecdata[,4]*1500) >= 1000
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

pdf(file="~/git/heteropecilly/pair_stats/align_pair_stats_CS.pdf", width=8, height=7)
par(mar=c(4.5,4.5,4,1) )
plot(refprotlength, ungappedpct, xlim=c(150,4700), col="#a08a2433", pch=16, frame.plot=FALSE, main="Fraction of each human protein kept in final alignments", xlab="Human reference protein length (AAs)", ylab="Fraction coverage in alignment", cex.axis=1.3, cex.lab=1.3)
refprothist = hist(refprotlength,breaks=seq(0,5000,25),plot=FALSE)
#lines(refprothist$breaks[2:length(refprothist$breaks)], refprothist$counts/max(refprothist$counts)/2, lwd=3 )
points(whelandata[,8], whelandata[,4], col="#d95f0233", pch=16)
points(borowiecdata[,8], borowiecdata[,4], col="#1b9e7733", pch=16)
#points(philippedata[,8], philippedata[,4], col="#e7298a33", pch=16)

text(refprotlength[simionlongprots]+10, ungappedpct[simionlongprots], simion_shortname[simionlongprots], pos=4, col="#a08a24")
text(whelandata[,8][whelan_longprots]+10, whelandata[,4][whelan_longprots], whelan_shortname[whelan_longprots], pos=4, col="#d95f02")
text(borowiecdata[,8][borowiec_longprots]+10, borowiecdata[,4][borowiec_longprots], borowiec_shortname[borowiec_longprots], pos=4, col="#1b9e77")

legend(3200,1.01,legend=c("Whelan 2017","Simion 2017", "Borowiec 2015"), pch=15, col=c("#d95f02","#a08a24", "#1b9e77"), cex=1.3 )
#legend(3000,1.01,legend=c("Whelan 2017","Simion 2017", "Borowiec 2015", "Philippe 2009"), pch=15, col=c("#d95f02","#2211cc", "#1b9e77", "#e7298a" ), cex=1.3 )

text(3300,0.70,"212 of 212 prots", col=c("#d95f02"), cex=1.3, pos=4)
text(3900,0.65, paste(round(whelansitefraction,digits=3),"of reference sites"), col=c("#d95f02"), cex=1.3)

text(3300,0.58,"1499 of 1719 prots", col=c("#a08a24"), cex=1.3, pos=4)
text(3900,0.53, paste(round(simionsitefraction,digits=3),"of reference sites"), col=c("#a08a24"), cex=1.3)

text(3300,0.46,"1056 of 1080 prots", col=c("#1b9e77"), cex=1.3, pos=4)
text(3900,0.41, paste(round(borowiecsitefraction,digits=3),"of reference sites"), col=c("#1b9e77"), cex=1.3)

#text(3000,0.34,"128 of 128 prots", col=c("#e7298a"), cex=1.3, pos=4)
#text(3500,0.29, paste(round(philippesitefraction,digits=3),"of reference sites"), col=c("#e7298a"), cex=1.3)

#text(3600,0.40,paste(length(sharedprots),"proteins\ncommon to both sets"), col=c("#881188"), cex=1.3)
#par( fig = c(grconvertX(c(3000,4200), from="user", to="ndc"), grconvertY(c(0.1,0.4), from="user", to="ndc")), new=TRUE, mar=c(0,0,0,0))
#subvp = viewport(x=0.57, y=0.1, width=0.33, height=0.3)
#grid.arrange( gTree(triplevenn), gTree(triplevenn), gTree( triplevenn ) , ncol=3)
#grid.draw(triplevenn)
dev.off()

png(file="~/git/heteropecilly/pair_stats/align_pair_stats_simion-v-whelan.png", width=600, height=500)
plot(refprotlength, ungappedpct, col="#a08a2466", pch=16, main="Fraction of each human protein kept in final alignments", xlab="Human reference protein length (AAs)", ylab="Fraction coverage in alignment", cex.axis=1.3, cex.lab=1.3)
refprothist = hist(refprotlength,breaks=seq(0,5000,25),plot=FALSE)
#lines(refprothist$breaks[2:length(refprothist$breaks)], refprothist$counts/max(refprothist$counts)/2, lwd=3 )

points(whelandata[,8], whelandata[,4], col="#d95f0266", pch=16)
points(borowiecdata[,8], borowiecdata[,4], col="#1b9e7766", pch=16)

legend(3000,1.0,legend=c("Whelan 2017","Simion 2017", "Borowiec 2015"), pch=15, col=c("#d95f02","#a08a24", "#1b9e77" ), cex=1.3 )

text(3000,0.70,"212 of 212 prots", col=c("#d95f02"), cex=1.3, pos=4)
text(3500,0.65, paste(round(whelansitefraction,digits=3),"of reference sites"), col=c("#d95f02"), cex=1.3)

text(3000,0.58,"1499 of 1719 prots", col=c("#a08a24"), cex=1.3, pos=4)
text(3500,0.53, paste(round(simionsitefraction,digits=3),"of reference sites"), col=c("#a08a24"), cex=1.3)

text(3000,0.46,"1056 of 1080 prots", col=c("#1b9e77"), cex=1.3, pos=4)
text(3500,0.41, paste(round(borowiecsitefraction,digits=3),"of reference sites"), col=c("#1b9e77"), cex=1.3)
dev.off()




ryandata = read.table("~/git/heteropecilly/pair_stats/ryan2013-EST_align_pair_stats.tab", header=TRUE, sep="\t")
ryan_namesplit = sapply(as.character(ryandata[,2]), function(x) unlist(strsplit(x, split='|', fixed=TRUE))[3] )
ryan_shortname = sapply(as.character(ryan_namesplit), function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1] )
ryan_longprots = (ryandata[,8] - ryandata[,4]*1500) >= 400
ryansitefraction = sum(ryandata[,3])/sum(ryandata[,8])

whelan_d1_data = read.table("~/git/heteropecilly/pair_stats/whelan2015-d1opi_align_pair_stats.tab", header=TRUE, sep="\t")
whelan_d1_namesplit = sapply(as.character(whelan_d1_data[,2]), function(x) unlist(strsplit(x, split='|', fixed=TRUE))[3] )
whelan_d1_shortname = sapply(as.character(whelan_d1_namesplit), function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1] )
whelan_d1_longprots = (whelan_d1_data[,8] - whelan_d1_data[,4]*2000) >= 800
wheland1_sitefraction = sum(whelan_d1_data[,3])/sum(whelan_d1_data[,8])

cannon_data = read.table("~/git/heteropecilly/pair_stats/cannon2016_align_pair_stats.tab", header=TRUE, sep="\t")
cannon_namesplit = sapply(as.character(cannon_data[,2]), function(x) unlist(strsplit(x, split='|', fixed=TRUE))[3] )
cannon_shortname = sapply(as.character(cannon_namesplit), function(x) unlist(strsplit(x, split='_', fixed=TRUE))[1] )
cannon_longprots = (cannon_data[,8] - cannon_data[,4]*2000) >= 800
cannonsitefraction = sum(cannon_data[,3])/sum(cannon_data[,8])

source("~/R_packages/vioplot.R") # use modified violin plot source code
violincolors =  c("#2211cc" , "#d95f02" , "#1b9e77", "#cc1818", "#a08a24", "#e7298a")
violinnames = c("Ryan\n2013", "Whelan\n2015", "Borowiec\n2015", "Cannon\n2016", "Simion\n2017", "Philippe\n2009")

pdf(file="~/git/heteropecilly/pair_stats/align_pair_stats_RWB.pdf", width=9, height=6)
layout(mat=matrix(c(1,1,2,1,1,3),2,3,byrow=TRUE),heights=c(2,1,1),widths=c(0.8,1,1))
par(mar=c(4.5,4.5,4,0) )
plot(ryandata[,8], ryandata[,4], xlim=c(100,4400), ylim=c(0,1), col="#2211cc33", pch=16, frame.plot=FALSE, main="Fraction of each human protein kept in final alignments", xlab="Human reference protein length (AAs)", ylab="Fraction coverage in alignment", cex.axis=1.4, cex.lab=1.4)
points(whelan_d1_data[,8], whelan_d1_data[,4], col="#d95f0233", pch=16)
points(borowiecdata[,8], borowiecdata[,4], col="#1b9e7733", pch=16)

text(ryandata[,8][ryan_longprots]+10, ryandata[,4][ryan_longprots], ryan_shortname[ryan_longprots], pos=4, col="#2211cc")
text(whelan_d1_data[,8][whelan_d1_longprots]+10, whelan_d1_data[,4][whelan_d1_longprots], whelan_d1_shortname[whelan_d1_longprots], pos=4, col="#d95f02")
text(borowiecdata[,8][borowiec_longprots]+10, borowiecdata[,4][borowiec_longprots], borowiec_shortname[borowiec_longprots], pos=4, col="#1b9e77")

legend(2800,1.01,legend=c("Ryan 2013", "Whelan 2015", "Borowiec 2015"), pch=15, col=c("#2211cc", "#d95f02", "#1b9e77"), cex=1.5 )

#text(3300,0.70,"396 of 406 prots", col=c("#2211cc"), cex=1.3, pos=4)
#text(3900,0.65, paste(round(ryansitefraction,digits=3),"of reference sites"), col=c("#2211cc"), cex=1.3)
#text(3300,0.58,"248 of 251 prots", col=c("#d95f02"), cex=1.3, pos=4)
#text(3900,0.53, paste(round(whelansitefraction,digits=3),"of reference sites"), col=c("#d95f02"), cex=1.3)
#text(3300,0.46,"1056 of 1080 prots", col=c("#1b9e77"), cex=1.3, pos=4)
#text(3900,0.41, paste(round(borowiecsitefraction,digits=3),"of reference sites"), col=c("#1b9e77"), cex=1.3)
mtext("A",3,at=-100, cex=2)

par(mar=c(4,1,4,1))
vioplot(ryandata[,4], whelan_d1_data[,4], borowiecdata[,4], cannon_data[,4], simiondata[,4], philippedata[,4], col=violincolors, names=violinnames, cex.lab=0.84)
mtext("Fraction coverage in alignment",2,line=2.5, cex=0.9)
mtext("B",3,at=-1, cex=2)

plot(0,0,type='n',axes=FALSE, frame.plot=FALSE,xlab="", ylab="")
mtext("C",3,at=-1.5, line=3, cex=2)
dev.off()

#library(vioplot)

pdf(file="~/git/heteropecilly/pair_stats/align_pair_stats_violin_plots.pdf", width=8, height=7)
par(mar=c(4,4,4,1))
vioplot(ryandata[,4], whelan_d1_data[,4], borowiecdata[,4], cannon_data[,4], simiondata[,4], philippedata[,4], col=violincolors, names=violinnames, cex.lab=1)
mtext("Fraction coverage in alignment",2,line=2.5,cex=1.3)
dev.off()

w1_in_r_prots = match(whelan_d1_data[,2], ryandata[,2])
w1_in_r_prots = w1_in_r_prots[!is.na(w1_in_r_prots)]

w1_in_b_prots = match(whelan_d1_data[,2], borowiecdata[,2])
w1_in_b_prots = w1_in_b_prots[!is.na(w1_in_b_prots)]

r_in_b_prots = match(ryandata[,2], borowiecdata[,2])
r_in_b_prots = r_in_b_prots[!is.na(r_in_b_prots)]

w1_in_r_in_b = match(ryandata[,2][w1_in_r_prots], borowiecdata[,2])
w1_in_r_in_b = w1_in_r_in_b[!is.na(w1_in_r_in_b)]

pdf(file="~/git/heteropecilly/pair_stats/align_pair_stats_venn_RWB.pdf", width=6, height=6)
triplevenn = draw.triple.venn( area1=length(ryandata[,8]) , area2=length(whelan_d1_data[,8]), area3=length(borowiecdata[,8]), n12=length(w1_in_r_prots), n23=length(w1_in_b_prots) , n13=length(r_in_b_prots), n123=length(w1_in_r_in_b) ,  fill=c( "#2211cc33" , "#cc181833" , "#18cc1833"), category=c("Ryan 2013","Whelan 2015","Borowiec 2015") , cex=2, cat.cex=2, cat.pos=c(-20,20,180) )
dev.off()




pdf(file="~/git/heteropecilly/pair_stats/align_pair_stats_CSP.pdf", width=8, height=7)
par(mar=c(4.5,4.5,4,1) )
plot(refprotlength, ungappedpct, xlim=c(150,4700), col="#a08a2433", pch=16, frame.plot=FALSE, main="Fraction of each human protein kept in final alignments", xlab="Human reference protein length (AAs)", ylab="Fraction coverage in alignment", cex.axis=1.3, cex.lab=1.3)
points(cannon_data[,8], cannon_data[,4], col="#cc181833", pch=16)
#points(philippedata[,8], philippedata[,4], col="#e7298a33", pch=16)
text(refprotlength[simionlongprots]+10, ungappedpct[simionlongprots], simion_shortname[simionlongprots], pos=4, col="#a08a24")
text(cannon_data[,8][cannon_longprots]+10, cannon_data[,4][cannon_longprots], cannon_shortname[cannon_longprots], pos=4, col="#cc1818")
legend(3200,1.01,legend=c("Cannon 2016","Simion 2017"), pch=15, col=c("#cc1818","#a08a24"), cex=1.3 )
#legend(3200,1.01,legend=c("Cannon 2016","Simion 2017", "Philippe 2009"), pch=15, col=c("#cc1818","#a08a24", "#e7298a"), cex=1.3 )
text(3300,0.70,"212 of 212 prots", col=c("#cc1818"), cex=1.3, pos=4)
text(3900,0.65, paste(round(cannonsitefraction,digits=3),"of reference sites"), col=c("#cc1818"), cex=1.3)
text(3300,0.58,"1499 of 1719 prots", col=c("#a08a24"), cex=1.3, pos=4)
text(3900,0.53, paste(round(simionsitefraction,digits=3),"of reference sites"), col=c("#a08a24"), cex=1.3)
dev.off()



pdf(file="~/git/heteropecilly/pair_stats/align_pairs_WSB_venn.pdf", width=7, height=7)
triplevenn = draw.triple.venn( area1=length(whelandata[,8]) , area2=length(simiondata[,8]), area3=length(borowiecdata[,8]), n12=length(w_in_s_prots), n23=length(s_in_b_prots) , n13=length(w_in_b_prots), n123=length(w_in_s_in_b) ,  fill=c( "#cc181833" , "#2211cc33" , "#18cc1833"), category=c("Whelan 2017","Simion 2017","Borowiec 2015") , cex=2, cat.cex=2, cat.pos=c(-20,20,180) )
dev.off()

png(file="~/git/heteropecilly/pair_stats/align_pairs_WSB_venn.png", width=450, height=450)
triplevenn = draw.triple.venn( area1=length(whelandata[,8]) , area2=length(simiondata[,8]), area3=length(borowiecdata[,8]), n12=length(w_in_s_prots), n23=length(s_in_b_prots) , n13=length(w_in_b_prots), n123=length(w_in_s_in_b) ,  fill=c( "#cc181833" , "#2211cc33" , "#18cc1833"), category=c("Whelan 2017","Simion 2017","Borowiec 2015") , cex=2, cat.cex=2, cat.pos=c(-20,20,180) )
dev.off()










#
