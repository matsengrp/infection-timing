{library(ggplot2)
library(gridExtra)
}
#read in data files. Path my need to be altered to reflect current file locations. 
a<-read.csv("~/Desktop/Projects/Molecular Clock/DataAnalysis/Run3March2020/AllRunsAvg.csv")
w<-read.csv("~/Desktop/VLAPD.csv")
s<-read.csv("~/Desktop/VLbyregion.csv")

#sampleplot to check library is in
ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()

###20200325 VL v APD ANALYSIS FOR RUNS 1-3 
#
##all VLs good and failed samples all regions 1x3 scatter. #success indicated by point color
#{a<-ggplot(data=s, aes(x=s$Month, y=s$VLall))+
#  geom_point()+
#  geom_point(aes(y=s$F1VLgood, color='F1 sequenced'))+
#  geom_point(aes(y=s$F1VLbad, color='F1 sequence fail'))+
#  scale_y_continuous(trans = "log10")+
#  labs(x="F1", y="Viral Load")+
#  theme(legend.position = "none")
#
##VL good and failed samples F2
#b<-ggplot(data=s, aes(x=s$Month, y=s$VLall))+
#  geom_point()+
#  geom_point(aes(y=s$F2VLgood, color='F2 sequenced'))+
#  geom_point(aes(y=s$F2VLbad, color='F2 sequence fail'))+
#  scale_y_continuous(trans = "log10")+
#  labs(x="F2",y="")+
#  theme(legend.position = "none", axis.title.y =element_blank#())
#
##VL good and failed samples F3
#c<-ggplot(data=s, aes(x=s$Month, y=s$VLall))+
#  geom_point()+
#  geom_point(aes(y=s$F3VLgood, color='F3 sequenced'))+
#  geom_point(aes(y=s$F3VLbad, color='F3 sequence fail'))+
#  scale_y_continuous(trans = "log10")+
#  labs(x="F3")+
#  theme(axis.title.y =element_blank())
#  
#
#grid.arrange(a,b,c,ncol=3, bottom="Time of Infection (month)")
#}
#
##VLvTOI new sheet loess fit log transformed each patient #separately (11 plots)
#{ggplot(data=w, aes(x=w$Month,y=w$VL, group=IDNUM))+
#  geom_point(aes(color=IDNUM))+
#  geom_smooth(method = 'loess', se=FALSE, aes(color=IDNUM))+
#  labs(x='Time of Infection (Month)', y='Viral Load')+
#  scale_y_continuous(trans = 'log10')+
#  facet_wrap(~w$IDNUM, ncol=4)
#}  
#
#
##VLvTOI new sheet loess fit NOT log transformed each patient #separately (11 plots)
#{ggplot(data=w, aes(x=w$Month,y=w$VL, group=IDNUM))+
#  geom_point(aes(color=IDNUM))+
#  geom_smooth(method = 'loess', se=FALSE, aes(color=IDNUM))+
#  labs(x='Time of Infection (Month)', y='Viral Load')+
#  coord_cartesian(ylim=c(0,50000000))+
#  facet_wrap(~w$IDNUM, ncol=4)
#  
#}
#
##ALL VLvTOI loess fit colored by patient no log
#{ggplot(data=w, aes(x=w$Month,y=w$VL, group=IDNUM))+
#  geom_point(aes(color=IDNUM))+
#  geom_smooth(method = 'loess', se=FALSE, aes(color=IDNUM))+
#  labs(x='Time of Infection (Month)', y='Viral Load')+
#  coord_cartesian(ylim=c(0,50000000))
#}
#
##ALL VLvTOI loess fit colored by patient no log with failed #samples indicated by black 
#{ggplot(data=s, aes(x=s$Month,y=s$VLall, group=IDNUM))+
#    geom_point(aes(y=s$F1VLbad))+
#    geom_point(aes(y=s$F2VLbad))+
#    geom_point(aes(y=s$F3VLbad))+
#    geom_smooth(method = 'loess', se=FALSE, aes(color=IDNUM))+
#    labs(x='Time of Infection (Month)', y='Viral Load')+
#    coord_cartesian(ylim=c(0,50000000))
#}
#
##ALL VLvTOI loess fit colored by patient log transformed
#{ggplot(data=w, aes(x=w$Month,y=w$VL, group=IDNUM))+
#    geom_point(aes(color=IDNUM))+
#    geom_smooth(method = 'loess', se=FALSE, aes(color=IDNUM))+
#    labs(x='Time of Infection (Month)', y='Viral Load')+
#    scale_y_continuous(trans = "log10")
#}
#
##ALL VLvTOI loess for colored by patient log transformed with #failed samples indicated by black
#{ggplot(data=s, aes(x=s$Month,y=s$VLall, group=s$IDNUM))+
#    geom_point(aes(y=s$F1VLbad))+
#    geom_point(aes(y=s$F2VLbad))+
#    geom_point(aes(y=s$F3VLbad))+
#    geom_smooth(method = 'loess', se=FALSE, aes(color=IDNUM))+
#    labs(x='Time of Infection (Month)', y='Viral Load')+
#    scale_y_continuous(trans = "log10")
#}
#
##APD1 all regions split by patient (11 plots)
#{ggplot(data=w, aes(x=w$Month,y=w$F1APD1, group=IDNUM))+
#  geom_point(aes(y=w$F1APD1,color="F1APD1"))+
#  geom_point(aes(y=w$F2APD1,color="F2APD1"))+
#  geom_point(aes(y=w$F3APD1,color="F3APD1"))+
#  geom_smooth(method = 'lm', se=FALSE, aes(y=w$F1APD1, #color="F1APD1"))+
#  geom_smooth(method = 'lm', se=FALSE, aes(y=w$F2APD1, #color="F2APD1"))+
#  geom_smooth(method = 'lm', se=FALSE, aes(y=w$F3APD1, #color="F3APD1"))+
#  labs(x='Time of Infection (Month)', y='APD1')+
#  facet_wrap(~w$IDNUM, ncol=4)+
#  scale_y_continuous(trans="log10")
#}
#
##VL v all APDs colored by region log transformed 
#{ggplot(data=w, aes(x=w$VL,y=w$F1APD1))+
#  geom_point(aes(y=w$F1APD1,color="F1APD1"))+
#  geom_point(aes(y=w$F2APD1,color="F2APD1"))+
#  geom_point(aes(y=w$F3APD1,color="F3APD1"))+
#  labs(x='log Viral Load', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans="log10")+
#  geom_smooth(method = 'lm', se=FALSE)
#}
#  
##VL v all APDs colored by region not log transformed 
#{ggplot(data=w, aes(x=w$VL,y=w$F1APD1))+
#    geom_point(aes(y=w$F1APD1,color="F1APD1"))+
#    geom_point(aes(y=w$F2APD1,color="F2APD1"))+
#    geom_point(aes(y=w$F3APD1,color="F3APD1"))+
#    labs(x='log Viral Load', y='APD 1%')+
#    coord_cartesian(ylim=c(0,0.025), xlim = c(0, 7500000))+
#    geom_smooth(method = 'lm', se=FALSE)
#}
#
##VL v F1 APDs
#{ggplot(data=w, aes(x=w$VL,y=w$F1APD1))+
#  geom_point(aes(y=w$F1APD1,color="F1APD1"))+
#  labs(x='log Viral Load', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans="log10")+
#  geom_smooth(method = 'lm', se=FALSE)
#}
#  
##VL v F2 APDs
#{ggplot(data=w, aes(x=w$VL,y=w$F2APD1))+
#  geom_point(aes(y=w$F1APD1,color="F2APD1"))+
#  labs(x='log Viral Load', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans="log10")+
#  geom_smooth(method = 'lm', se=FALSE)
#}
#  
##VL v F3 APDs
#{ggplot(data=w, aes(x=w$VL,y=w$F3APD1))+
#  geom_point(aes(y=w$F1APD1,color="F3APD1"))+
#  labs(x='log Viral Load', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans="log10")+
#  geom_smooth(method = 'lm', se=FALSE)
#}
#
##VL v F1, F2, F3 together VL log transformed (1x3 plot)
#{
##VL v F1 APDs
#a<-ggplot(data=w, aes(x=w$VL,y=w$F1APD1))+
#  geom_point(aes(y=w$F1APD1))+
#  labs(x='F1', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans="log10")+
#  geom_smooth(method = 'lm', se=FALSE)
#
##VL v F2 APDs
#b<-ggplot(data=w, aes(x=w$VL,y=w$F2APD1))+
#  geom_point(aes(y=w$F1APD1))+
#  labs(x='F2', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans="log10")+
#  geom_smooth(method = 'lm', se=FALSE)
#
##VL v F3 APDs
#c<-ggplot(data=w, aes(x=w$VL,y=w$F3APD1))+
#  geom_point(aes(y=w$F1APD1))+
#  labs(x='F3', y='APD 1%')+
#  coord_cartesian(ylim=c(0,0.025), xlim = c(10000, 10000000))+
#  scale_x_continuous(trans=)+
#  geom_smooth(method = 'lm', se=FALSE)

#grid.arrange(a,b,c,ncol=3, top="VL v APD by Region", #bottom="log Viral Load")
#  }

##APD1s over time not seprated by patient (1x3 plot)
#{
#  #VL v F1 APDs
#  a<-ggplot(data=w, aes(x=w$Month,y=w$F1APD1))+
#    geom_point(aes(y=w$F1APD1))+
#    labs(x='F1', y='APD 1%')+
#    coord_cartesian(ylim=c(0,0.025), xlim = c(0, 24))+
#    geom_smooth(method = 'lm', se=FALSE)
#  
#  #VL v F2 APDs
#  b<-ggplot(data=w, aes(x=w$Month,y=w$F2APD1))+
#    geom_point(aes(y=w$F1APD1))+
#    labs(x='F2', y='APD 1%')+
#    coord_cartesian(ylim=c(0,0.025), xlim = c(0, 24))+
#    geom_smooth(method = 'lm', se=FALSE)
#  
#  #VL v F3 APDs
#  c<-ggplot(data=w, aes(x=w$Month,y=w$F3APD1))+
#    geom_point(aes(y=w$F3APD1))+
#    labs(x='F3', y='APD 1%')+
#    coord_cartesian(ylim=c(0,0.025), xlim = c(0, 24))+
#    geom_smooth(method = 'lm', se=FALSE)
#  
#  grid.arrange(a,b,c,ncol=3, top="APD over Time by Region", #bottom="Month")
#}
#

##20200314 ETI v TI ANALYSIS FOR RUNS 1-3

#separate data by fragment 
{Fragment<-split(a, a$Fragment)
F1<-Fragment$F1
F2<-Fragment$F2
F3<-Fragment$F3

#separate data by run 
Run<-split(a,a$Run)
Run1<-Run$`Run 1`
Run2<-Run$`Run 2`
Run3<-Run$`Run 3`

#separate data by frgament and run
Run1F<-split(Run1, Run1$Fragment)
Run1F1<-Run1F$F1
Run1F2<-Run1F$F2
Run1F3<-Run1F$F3

Run2F<-split(Run2, Run2$Fragment)
Run2F1<-Run2F$F1
Run2F2<-Run2F$F2
Run2F3<-Run2F$F3

Run3F<-split(Run3, Run3$Fragment)
Run3F1<-Run3F$F1
Run3F2<-Run3F$F2
Run3F3<-Run3F$F3



}

#perfect correlation plot 
{ggplot(data=Run3F1,aes(y=Run3F1$ActualTOI..year.,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,7),ylim=c(0,7))
} 


#plot VL v APD1 with equations 1x3 with logtrans for VL, colored by patient 
{
  m <- lm(F1$VL ~ F1$AvgAPD1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=F1,aes(y=AvgAPD1,x=VL))+
    geom_point(aes(color=Sample))+
    stat_smooth(method="lm",se=TRUE)+
    labs(x="F1", y="APD (1% threshold)")+ 
    theme(legend.position = "none")+coord_cartesian(xlim=c(1e+05, 1e+07), ylim=c(0,0.025))+
    annotate(geom="text", x=1e+06, y=0.025, label=texta)+
    scale_x_continuous(trans = 'log10')
  
  m <- lm(F2$VL ~ F2$AvgAPD1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textb)
  
  
  
  b1<-ggplot(data=F2,aes(y=AvgAPD1,x=VL))+
    geom_point(aes(color=Sample))+
    stat_smooth(method="lm",se=TRUE)+labs(x="F2", y="")+
    theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+
    coord_cartesian(xlim=c(1e+05, 1e+07), ylim=c(0,0.025))+
    annotate(geom="text", x=1e+06, y=0.025, label=textb)+
    scale_x_continuous(trans = 'log10')
  
  m <- lm(F3$VL ~ F3$AvgAPD1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  
  c1<-ggplot(data=F3,aes(y=AvgAPD1,x=VL))+
    geom_point(aes(color=Sample))+
    stat_smooth(method="lm",se=TRUE)+
    labs(x="F3", y="")+
    theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+
    coord_cartesian(xlim=c(1e+05, 1e+07), ylim=c(0,0.025))+annotate(geom="text", x=1e+06, y=0.025, label=textc)+
    scale_x_continuous(trans = 'log10')
  
  s<-grid.arrange(a1,b1,c1,ncol=3, bottom='VL')
}

#plot VL over time per patient loess ONLY T1, T2, T3
{a<-ggplot(data=F1,aes(y=VL,x=ActualTOI..Month.,group=Sample))+
    geom_point(aes(color=Sample))+
    stat_smooth(method="loess",se=FALSE, aes(color=Sample))+
    labs(x="", y="Viral Load")+ 
    coord_cartesian(xlim=c (0,24))
  
  grid.arrange(a,ncol=1, top="All Runs VL v Time of Infection", bottom="Time of Infection (month)")
}

#plot VL over time per patient log transformed loess ONLY T1, T2, T3
{a<-ggplot(data=F1,aes(y=VL,x=ActualTOI..Month.,group=Sample))+
    geom_point(aes(color=Sample))+
    stat_smooth(method="loess",se=FALSE, aes(color=Sample))+
    labs(x="", y="Viral Load")+ 
    coord_cartesian(xlim=c (0,24))+
    scale_y_continuous(trans = 'log10')
  
  grid.arrange(a,ncol=1, top="All Runs VL v Time of Infection", bottom="Time of Infection (month)")
}

#plot VL v APD1 with equations 1x3 
{
  m <- lm(F1$VL ~ F1$AvgAPD1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  a1<-ggplot(data=F1,aes(y=AvgAPD1,x=VL))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F1", y="APD (1% threshold)")+ theme(legend.position = "none")+coord_cartesian(ylim=c(0,0.025))+annotate(geom="text", x=3, y=6, label=texta)
  
  m <- lm(F2$VL ~ F2$AvgAPD1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  
  
  b1<-ggplot(data=F2,aes(y=AvgAPD1,x=VL))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(ylim=c(0,0.025))+annotate(geom="text", x=3, y=6, label=textb)
  
  m <- lm(F3$VL ~ F3$AvgAPD1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  
  c1<-ggplot(data=F3,aes(y=AvgAPD1,x=VL))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F3", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(ylim=c(0,0.025))+annotate(geom="text", x=3, y=6, label=textc)
  
  s<-grid.arrange(a1,b1,c1,ncol=3, top='VL v APD1', bottom='Viral Load')
}


 
#plot Run 1, 2, 3, and all Avg ETI v TI for APD 1 with equations 4x3
{
#Run1
m <- lm(Run1F1$ActualTOI..year. ~ Run1F1$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)

a1<-ggplot(data=Run1F1,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("Run 1")+labs(x="", y="Calculated TOI (years)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=texta)

m <- lm(Run1F2$ActualTOI..year. ~ Run1F2$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textb <- paste("y = ",b,"x + ",a, sep="")
print(textb)

b1<-ggplot(data=Run1F2,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("  ")+labs(x="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textb)

m <- lm(Run1F3$ActualTOI..year. ~ Run1F3$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textc <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=Run1F3,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("  ")+labs(x="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textc)

q<-grid.arrange(a1,b1,c1,ncol=3)

#Run2
m <- lm(Run2F1$ActualTOI..year. ~ Run2F1$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)

a1<-ggplot(data=Run2F1,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("Run 2")+labs(x="", y="Calculated TOI (years)")+theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=texta)

m <- lm(Run2F2$ActualTOI..year. ~ Run2F2$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textb <- paste("y = ",b,"x + ",a, sep="")
print(textb)

b1<-ggplot(data=Run2F2,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textb)

m <- lm(Run2F3$ActualTOI..year. ~ Run2F3$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textc <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=Run2F3,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textc)

r<-grid.arrange(a1,b1,c1,ncol=3)

#Run3
m <- lm(Run3F1$ActualTOI..year. ~ Run3F1$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)

a1<-ggplot(data=Run3F1,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("Run 3")+labs(x="", y="Calculated TOI (years)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=texta)

m <- lm(Run3F2$ActualTOI..year. ~ Run3F2$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textb <- paste("y = ",b,"x + ",a, sep="")
print(textb)

b1<-ggplot(data=Run3F2,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textb)

m <- lm(Run3F3$ActualTOI..year. ~ Run3F3$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textc <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=Run3F3,aes(y=AvgTOI1,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textc)

s<-grid.arrange(a1,b1,c1,ncol=3)

#all runs 
m <- lm(F1$ActualTOI..year. ~ F1$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)

a1<-ggplot(data=F1,aes(AvgTOI1,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="red", fill="red")+ggtitle("All Runs")+labs(x="F1", y="Calculated TOI (years)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=texta)

m <- lm(F2$ActualTOI..year. ~ F2$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textb <- paste("y = ",b,"x + ",a, sep="")
print(textb)  

b1<-ggplot(data=F2,aes(AvgTOI1,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orangered1", fill="orangered1")+ggtitle("  ")+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textb)

m <- lm(F3$ActualTOI..year. ~ F3$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textc <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=F3,aes(AvgTOI1,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orange1", fill="orange1")+ggtitle("  ")+labs(x="F3", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textc)

t<-grid.arrange(a1,b1,c1,ncol=3)

#plot all together 4x4

grid.arrange(q,r,s,t,nrow=4, top="All Runs Calculated TOI (APD 1%) v Actual TOI", bottom="Time of Infection (years)")
}

#plot Run 1, 2, 3, and all Avg ETI v TI for APD 10 with equations 4x3
{
  #Run1
  m <- lm(Run1F1$ActualTOI..year. ~ Run1F1$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=Run1F1,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("Run 1")+labs(x="", y="Calculated TOI (years)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=texta)
  
  m <- lm(Run1F2$ActualTOI..year. ~ Run1F2$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textb <- paste("y = ",b,"x + ",a, sep="")
  print(textb)
  
  b1<-ggplot(data=Run1F2,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("  ")+labs(x="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textb)
  
  m <- lm(Run1F3$ActualTOI..year. ~ Run1F3$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textc <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  c1<-ggplot(data=Run1F3,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("  ")+labs(x="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textc)
  
  q<-grid.arrange(a1,b1,c1,ncol=3)
  
  #Run2
  m <- lm(Run2F1$ActualTOI..year. ~ Run2F1$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=Run2F1,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("Run 2")+labs(x="", y="Calculated TOI (years)")+theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=texta)
  
  m <- lm(Run2F2$ActualTOI..year. ~ Run2F2$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textb <- paste("y = ",b,"x + ",a, sep="")
  print(textb)
  
  b1<-ggplot(data=Run2F2,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textb)
  
  m <- lm(Run2F3$ActualTOI..year. ~ Run2F3$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textc <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  c1<-ggplot(data=Run2F3,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textc)
  
  r<-grid.arrange(a1,b1,c1,ncol=3)
  
  #Run3
  m <- lm(Run3F1$ActualTOI..year. ~ Run3F1$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=Run3F1,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+ggtitle("Run 3")+labs(x="", y="Calculated TOI (years)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=texta)
  
  m <- lm(Run3F2$ActualTOI..year. ~ Run3F2$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textb <- paste("y = ",b,"x + ",a, sep="")
  print(textb)
  
  b1<-ggplot(data=Run3F2,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textb)
  
  m <- lm(Run3F3$ActualTOI..year. ~ Run3F3$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textc <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  c1<-ggplot(data=Run3F3,aes(y=AvgTOI10,x=ActualTOI..year.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ggtitle("  ")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,7))+annotate(geom="text", x=1.5, y=6, label=textc)
  
  s<-grid.arrange(a1,b1,c1,ncol=3)
  
  #all runs 
  m <- lm(F1$ActualTOI..year. ~ F1$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=F1,aes(AvgTOI10,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="red", fill="red")+ggtitle("All Runs")+labs(x="F1", y="Calculated TOI (years)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=texta)
  
  m <- lm(F2$ActualTOI..year. ~ F2$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textb <- paste("y = ",b,"x + ",a, sep="")
  print(textb)  
  
  b1<-ggplot(data=F2,aes(AvgTOI10,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orangered1", fill="orangered1")+ggtitle("  ")+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textb)
  
  m <- lm(F3$ActualTOI..year. ~ F3$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  textc <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  c1<-ggplot(data=F3,aes(AvgTOI10,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orange1", fill="orange1")+ggtitle("  ")+labs(x="F3", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textc)
  
  t<-grid.arrange(a1,b1,c1,ncol=3)
  
  #plot all together 4x4
  
  grid.arrange(q,r,s,t,nrow=4, top="All Runs Calculated TOI (APD 10%) v Actual TOI", bottom="Time of Infection (years)")
}


#plot F1, F2, F3 Avg ETI Run 1 APD1 with equations 1x3
{
m <- lm(Run1F1$ActualTOI..Month. ~ Run1F1$AvgTOI10)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(textc)

a1<-ggplot(data=Run1F1,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F1", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=texta)

m <- lm(Run1F2$ActualTOI..Month. ~ Run1F2$AvgTOI10)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(textc)



b1<-ggplot(data=Run1F2,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textb)

m <- lm(Run1F3$ActualTOI..Month. ~ Run1F3$AvgTOI10)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(textc)


c1<-ggplot(data=Run1F3,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F3", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textc)

s<-grid.arrange(a1,b1,c1,ncol=3)
}

#plot APDs over time (months) for Run1 samples all fragments 3x3 funky legend no equations 
{a<-ggplot(data=Run1F1,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
b<-ggplot(data=Run1F2,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
c<-ggplot(data=Run1F3,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))

q<-grid.arrange(a,b,c,ncol=3)


#plot APD5 over time (months) for Run1 samples all fragments  
a<-ggplot(data=Run1F1,aes(y=Run1F1$APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run1F2,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run1F3,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  r<-grid.arrange(a,b,c,ncol=3)


#plot APD10 over time (months) for Run1 samples all fragments  
a<-ggplot(data=Run1F1,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run1F2,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run1F3,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  s<-grid.arrange(a,b,c,ncol=3)



grid.arrange(q,r,s,nrow=3, top="Run 1 APDs v Time of Infection", bottom="Time of Infection (months)")
}

#plot APD1 over time (months) for Run2 samples all fragments  3x3 funky legend no equations 
{a<-ggplot(data=Run2F1,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run2F2,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run2F3,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  q<-grid.arrange(a,b,c,ncol=3)


#plot APD5 over time (months) for Run2 samples all fragments  
a<-ggplot(data=Run2F1,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run2F2,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run2F3,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  r<-grid.arrange(a,b,c,ncol=3)


#plot APD10 over time (months) for Run2 samples all fragments  
a<-ggplot(data=Run2F1,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run2F2,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run2F3,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  s<-grid.arrange(a,b,c,ncol=3)


grid.arrange(q,r,s,nrow=3, top="Run 2 APDs v Time of Infection", bottom="Time of Infection (months)")

}

#plot APDs over time (months) for Run3 samples all fragments 3x3 funky legend no equations  
{a<-ggplot(data=Run3F1,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run3F2,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run3F3,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  q<-grid.arrange(a,b,c,ncol=3)


#plot APD5 over time (months) for Run3 samples all fragments  
a<-ggplot(data=Run3F1,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run3F2,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run3F3,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  r<-grid.arrange(a,b,c,ncol=3)


#plot APD10 over time (months) for Run3 samples all fragments  
a<-ggplot(data=Run3F1,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=Run3F2,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=Run3F3,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  s<-grid.arrange(a,b,c,ncol=3)
  
  grid.arrange(q,r,s,nrow=3, top="Run 3 APDs v Time of Infection", bottom="Time of Infection (months)")
  
}




#plot APD1 over time (months) for all samples all fragments 3x3 funky legend no equations 
{a<-ggplot(data=F1,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=F2,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=F3,aes(y=APD_1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.position="none",axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  q<-grid.arrange(a,b,c,ncol=3)

#plot APD5 over time (months) for all samples all fragments  
a<-ggplot(data=F1,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=F2,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=F3,aes(y=APD_5,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  r<-grid.arrange(a,b,c,ncol=3)


#plot APD10 over time (months) for all samples all fragments  
a<-ggplot(data=F1,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=F2,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=F3,aes(y=APD_10,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  s<-grid.arrange(a,b,c,ncol=3)
  
  grid.arrange(q,r,s,nrow=3, top="All Runs APDs v Time of Infection", bottom="Time of Infection (months)")
  
}




#plot Avg APDs over time (months) for all samples all fragments 3x3 funky legend no equations 
{a<-ggplot(data=F1,aes(y=AvgAPD1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=F2,aes(y=AvgAPD1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=F3,aes(y=AvgAPD1,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  q<-grid.arrange(a,b,c,ncol=3)


#plot Avg APD5 over time (months) for all samples all fragments #saved in excel file as APD2 instead of 5
a<-ggplot(data=F1,aes(y=AvgAPD2,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=F2,aes(y=AvgAPD2,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=F3,aes(y=AvgAPD2,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  r<-grid.arrange(a,b,c,ncol=3)


#plot Avg APD10 over time (months) for all samples all fragments #saved in excel file as APD2 instead of 5
a<-ggplot(data=F1,aes(y=AvgAPD3,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  b<-ggplot(data=F2,aes(y=AvgAPD3,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  c<-ggplot(data=F3,aes(y=AvgAPD3,x=ActualTOI..Month.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,0.025))
  
  s<-grid.arrange(a,b,c,ncol=3)
  grid.arrange(q,r,s,nrow=3, top="All Runs APDs v Time of Infection", bottom="Time of Infection (months)")
  
}



#plot Avg ETI APD1 over time (years) for Run 3 1x3 no equations 
{a<-ggplot(data=Run3F1,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  b<-ggplot(data=Run3F2,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  c<-ggplot(data=Run3F3,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  
  grid.arrange(a,b,c,ncol=3, top="Run 3 ETI (APD1) v Time of Infection", bottom="Time of Infection (months)")
}

#plot Avg ETI APD5 over time (years) for Run 3 1x3 no equations 
{a<-ggplot(data=Run3F1,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  b<-ggplot(data=Run3F2,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  c<-ggplot(data=Run3F3,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  
  grid.arrange(a,b,c,ncol=3, top="Run 3 ETI (APD5) v Time of Infection", bottom="Time of Infection (months)")
}

#plot Avg ETI APD10 over time (years) for Run 3 1x3 no equations 
{a<-ggplot(data=Run3F1,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  b<-ggplot(data=Run3F2,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  c<-ggplot(data=Run3F3,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  
  grid.arrange(a,b,c,ncol=3, top="Run 3 ETI (APD10) v Time of Infection", bottom="Time of Infection (months)")
}




#plot Avg ETI over time (years) for Run 2 3x3 funky legend without equation 
{a<-ggplot(data=Run2F1,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="Calculated TOI (APD 1)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  b<-ggplot(data=Run2F2,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  c<-ggplot(data=Run2F3,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  
  q<-grid.arrange(a,b,c,ncol=3)


#plot Avg ETI APD5 over time (years) for Run 2
a<-ggplot(data=Run2F1,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="Calculated TOI (APD 5)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  b<-ggplot(data=Run2F2,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  c<-ggplot(data=Run2F3,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  
  r<-grid.arrange(a,b,c,ncol=3)


#plot Avg ETI APD10 over time (years) for Run 2
a<-ggplot(data=Run2F1,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F1", y="Calculated TOI (APD 10)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  b<-ggplot(data=Run2F2,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  c<-ggplot(data=Run2F3,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F3", y="")+theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))
  
  s<-grid.arrange(a,b,c,ncol=3)
  grid.arrange(q,r,s,nrow=3, top="Run 2 Calculated TOI v Time of Infection", bottom="Time of Infection (months)")
}

#plot Avg ETI APD1/5/10 over time (years) for Run 1 with equations 3x3 with equations 
{
  m <- lm(Run1F1$ActualTOI..Month. ~ Run1F1$AvgTOI1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=Run1F1,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="APD 1")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=texta)
  
  m <- lm(Run1F2$ActualTOI..Month. ~ Run1F2$AvgTOI1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textb)
  
  b1<-ggplot(data=Run1F2,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textb)
  
  m <- lm(Run1F3$ActualTOI..Month. ~ Run1F3$AvgTOI1)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  c1<-ggplot(data=Run1F3,aes(y=AvgTOI1,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textc)
  
  q<-grid.arrange(a1,b1,c1,ncol=3)


#plot Avg ETI APD5 over time (years) for Run 1

  m <- lm(Run1F1$ActualTOI..Month. ~ Run1F1$AvgTOI5)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  a1<-ggplot(data=Run1F1,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="APD 5")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=texta)
 
  m <- lm(Run1F2$ActualTOI..year. ~ Run1F2$AvgTOI5)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(texta)
  
  b1<-ggplot(data=Run1F2,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textb)
  
  m <- lm(Run1F3$ActualTOI..Month. ~ Run1F3$AvgTOI5)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  c1<-ggplot(data=Run1F3,aes(y=AvgTOI5,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textc)
  
  r<-grid.arrange(a1,b1,c1,ncol=3)


#plot Avg ETI APD10 over time (years) for Run 1

  m <- lm(Run1F1$ActualTOI..Month. ~ Run1F1$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  a1<-ggplot(data=Run1F1,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F1", y="APD 10")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=texta)
 
  m <- lm(Run1F2$ActualTOI..Month. ~ Run1F2$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  
  
  b1<-ggplot(data=Run1F2,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textb)
  
  m <- lm(Run1F3$ActualTOI..Month. ~ Run1F3$AvgTOI10)
  a <- signif(coef(m)[1], digits = 2)
  b <- signif(coef(m)[2], digits = 2)
  texta <- paste("y = ",b,"x + ",a, sep="")
  print(textc)
  
  
  c1<-ggplot(data=Run1F3,aes(y=AvgTOI10,x=ActualTOI..Month.))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=TRUE)+labs(x="F3", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,24),ylim=c(0,7))+annotate(geom="text", x=3, y=6, label=textc)
  
  s<-grid.arrange(a1,b1,c1,ncol=3)

grid.arrange(q,r,s,nrow=3, top="Run 1 Calculated TOI  v Actual TOI", bottom="Time of Infection (months)")
}



#plot ETI APD1 over time (years) for all samples all fragments by PATIENT 1x3  
{a<-ggplot(data=F1,aes(CalculatedTOI1,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="Calculated Time of Infection (years) (APD 1)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  b<-ggplot(data=F2,aes(CalculatedTOI1,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  c<-ggplot(data=F3,aes(CalculatedTOI1,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  
  grid.arrange(a,b,c,ncol=3, top="All Runs Estimated Time of Infection v Time of Infection", bottom="Time of Infection (years)")
}

#plot ETI APD5 over time (years) for all samples all fragments by PATIENT  1x3  
{a<-ggplot(data=F1,aes(CalculatedTOI5,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="Calculated Time of Infection (years) (APD 5)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  b<-ggplot(data=F2,aes(CalculatedTOI5,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  c<-ggplot(data=F3,aes(CalculatedTOI5,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  
  grid.arrange(a,b,c,ncol=3, top="All Runs Estimated Time of Infection v Time of Infection", bottom="Time of Infection (years)")
}

#plot ETI APD10 over time (years) for all samples all fragments by PATIENT 1x3  
{a<-ggplot(data=F1,aes(CalculatedTOI10,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="Calculated Time of Infection (years) (APD 10)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  b<-ggplot(data=F2,aes(CalculatedTOI10,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  c<-ggplot(data=F3,aes(CalculatedTOI10,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  
  grid.arrange(a,b,c,ncol=3, top="All Runs Estimated Time of Infection v Time of Infection", bottom="Time of Infection (years)")
}





#plot Avg ETI APD1 over time (years) for all samples all fragments by PATIENT 1x3 
{a<-ggplot(data=F1,aes(AvgTOI1,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="Calculated Time of Infection (years) (APD 1)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  b<-ggplot(data=F2,aes(AvgTOI1,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  c<-ggplot(data=F3,aes(AvgTOI1,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  
  grid.arrange(a,b,c,ncol=3, top="All Runs Estimated Time of Infection v Time of Infection", bottom="Time of Infection (years)")
}

#plot Avg ETI APD5 over time (years) for all samples all fragments by PATIENT  1x3  
{a<-ggplot(data=F1,aes(AvgTOI5,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="Calculated Time of Infection (years) (APD 5)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  b<-ggplot(data=F2,aes(AvgTOI5,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  c<-ggplot(data=F3,aes(AvgTOI5,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  
  grid.arrange(a,b,c,ncol=3, top="All Runs Estimated Time of Infection v Time of Infection", bottom="Time of Infection (years)")
}

#plot Avg ETI APD10 over time (years) for all samples all fragments by PATIENT   1x3  
{a<-ggplot(data=F1,aes(AvgTOI10,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="Calculated Time of Infection (years) (APD 10)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  b<-ggplot(data=F2,aes(AvgTOI10,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  c<-ggplot(data=F3,aes(AvgTOI10,x=ActualTOI..year.,group=Sample))+geom_point(aes(color=Sample))+stat_smooth(method="lm",se=FALSE, aes(color=Sample))+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))
  
  grid.arrange(a,b,c,ncol=3, top="All Runs Estimated Time of Infection v Time of Infection", bottom="Time of Infection (years)")
}


#plot Avg ETI all APDs over time (years) for all samples all fragments with equations :) 3x3
{
#linear regression equation   
m <- lm(F1$ActualTOI..year. ~ F1$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)
 
a1<-ggplot(data=F1,aes(AvgTOI1,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="red", fill="red")+labs(x="", y="Calc TOI (years) (APD 1)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=texta)

m <- lm(F2$ActualTOI..year. ~ F2$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textb <- paste("y = ",b,"x + ",a, sep="")
print(textb)  
  
b1<-ggplot(data=F2,aes(AvgTOI1,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orangered1", fill="orangered1")+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textb)
  
m <- lm(F3$ActualTOI..year. ~ F3$AvgTOI1)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textc <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=F3,aes(AvgTOI1,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orange1", fill="orange1")+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textc)
  
  q<-grid.arrange(a1,b1,c1,ncol=3)


#plot Avg ETI APD5 over time (years) for all samples all fragments with equations  
m <- lm(F1$ActualTOI..year. ~ F1$AvgTOI5)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)
a1<-ggplot(data=F1,aes(AvgTOI5,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="red", fill="red")+labs(x="", y="Calc TOI (years) (APD 5)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=texta)
 
m <- lm(F2$ActualTOI..year. ~ F2$AvgTOI5)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(textb)
b1<-ggplot(data=F2,aes(AvgTOI5,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orangered1", fill="orangered1")+labs(x="", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textb)

m <- lm(F3$ActualTOI..year. ~ F3$AvgTOI5)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=F3,aes(AvgTOI5,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orange1", fill="orange1")+labs(x="", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textc)
                                                                                                                                                                                                                                                                                                                       
  
  r<-grid.arrange(a1,b1,c1,ncol=3)


#plot Avg ETI APD10 over time (years) for all samples all fragments with equations :) 

m <- lm(F1$ActualTOI..year. ~ F1$AvgTOI10)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
texta <- paste("y = ",b,"x + ",a, sep="")
print(texta)

a1<-ggplot(data=F1,aes(AvgTOI10,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="red", fill="red")+labs(x="F1", y="Calc TOI (years) (APD 10)")+ theme(legend.position = "none")+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=texta)
  
m <- lm(F2$ActualTOI..year. ~ F2$AvgTOI10)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textb <- paste("y = ",b,"x + ",a, sep="")
print(textb)

b1<-ggplot(data=F2,aes(AvgTOI10,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orangered1", fill="orangered1")+labs(x="F2", y="")+ theme(legend.position = "none",axis.title.y=element_blank() , axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textb)
  
m <- lm(F3$ActualTOI..year. ~ F3$AvgTOI10)
a <- signif(coef(m)[1], digits = 2)
b <- signif(coef(m)[2], digits = 2)
textc <- paste("y = ",b,"x + ",a, sep="")
print(textc)

c1<-ggplot(data=F3,aes(AvgTOI10,x=ActualTOI..year.))+geom_point(aes(shape=Run))+stat_smooth(method="lm",se=TRUE,color="orange1", fill="orange1")+labs(x="F3", y="")+theme(legend.title=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_cartesian(xlim=c (0,2),ylim=c(0,6))+annotate(geom="text", x=1.5, y=6, label=textc)
  
  s<-grid.arrange(a1,b1,c1,ncol=3)
 
  
  grid.arrange(q,r,s,nrow=3, top="All Runs Calculated TOI v Actual TOI", bottom="Time of Infection (years)")

}



#Random Maggie code: 

# Create a plot for ETI vs. actual TI for APD for all runs facetted by fragment
plot_ETI_TI_regression(all_runs, fragment = "all", apd = 5, run = "all", facet = "Fragment", points = "Run")

# Create a plot for ETI vs. actual TI for APD1 for all runs facetted by fragment
plot_ETI_TI_regression(all_runs, fragment = "all", apd = 1, run = "all", facet = "Fragment", points = "Run")


```
theme(legend.position ="none") 
    text = c()
    for (frame in c("F1", "F2", "F3")){
        
        if (frame == "F1"){
            frame_data = F1
            index = 1
        } else if (frame == "F2"){
            frame_data = F2
            index = 2
        } else if (frame == "F3"){
            frame_data = F3
            index = 3
        }
        #linear regression equation   
        slope <- lm(frame_data[["ActualTOI..year."]] ~ frame_data[[avgETI]])
        a <- signif(coef(slope)[1], digits = 2)
        b <- signif(coef(slope)[2], digits = 2)
        text[index] = paste("y = ",b,"x + ",a, sep="")
    }
    dat_text <- data.frame(label = text, Fragment = c("F1", "F2", "F3"))
    title = paste0("ETI vs. Actual TI by Fragment for APD",apd)




