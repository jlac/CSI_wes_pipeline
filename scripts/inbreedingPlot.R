library(ggplot2)

fstat = read.table("inbreeding/picardMetrics.variant_calling_detail_metrics", header=T,skip=6)
png("inbreeding/Heterozygous_to_Homozygous_Ratio_mqc.png",width = 1000)
ggplot(data=fstat, aes(x=SAMPLE_ALIAS, y=HET_HOMVAR_RATIO)) + geom_bar(aes(color=HET_HOMVAR_RATIO<1.5),stat="identity",fill="steelblue") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7)) + ylim(0,3)
dev.off()

png("inbreeding/Mean_Homozygous_Tract_Length_mqc.png",width = 1000)
roh = read.table("inbreeding/ROH.hom.indiv", header=T)
roh2 = as.data.frame(roh)
roh2$barcolor[roh2$KBAVG<1900]<-"Normal"
roh2$barcolor[roh2$KBAVG>1900]<-"Elevated_Homozygosity"
roh2$barcolor[roh2$KBAVG>1975]<-"Consanguinity"
roh2$barcolor<-as.factor(roh2$barcolor)
categories <- c("Normal", "Elevated_Homozygosity", "Consanguinity")
cols <- c("grey","blue","red")
ggplot(data=roh2, aes(x=FID, y=KBAVG, col=barcolor, fill=barcolor)) + geom_bar(stat="identity") + scale_color_manual(values=c("Normal"="black","Elevated_Homozygosity"="black","Consanguinity"="black"), guide=FALSE,) + scale_fill_manual(values=c("Normal"="grey","Elevated_Homozygosity"="blue","Consanguinity"="red")) + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7),legend.position="bottom", legend.box = "horizontal") + labs(fill = "Status") + scale_x_discrete(name ="Phenotip ID") + scale_y_continuous(name ="Avg. Homoyzygous Window Length (KB)", limits=c(0,3500))
dev.off()