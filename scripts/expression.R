library("agricolae", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("scales", lib.loc="/usr/local/lib/R/site-library")
library("ggplot2", lib.loc="/usr/local/lib/R/site-library")
library("plyr", lib.loc="/usr/local/lib/R/site-library")
setwd("~/kairuku/Nanopore_seq/expression")
save.image("~/kairuku/Nanopore_seq/expression/expression.RData")

####import data(qPCR)
rawdata<-read.delim(file = "/home/yichun/kairuku/Nanopore_seq/expression/expression_raw", header = TRUE)
rawdata$Gene<-gsub("CC1G_01973", "CC2G_008991", rawdata$Gene)
rawdata$Gene<-gsub("CC1G_02204", "CC2G_002608", rawdata$Gene)
rawdata$Gene<-gsub("GSK3", "CC2G_005966", rawdata$Gene)
rawdata<-subset(rawdata, Stage == "M" | Stage == "K")
summary.expression<-ddply(rawdata, c("Stage", "Gene", "Strain"), 
                          summarise, 
                          mean = mean(Expression),
                          std = sd(Expression),
                          cv = sd(Expression)/mean(Expression)
)
summary.expression


mean.CC2G_008991<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_008991")
mean.CC2G_002608<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_002608")
mean.CC2G_005966<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_005966")
mean.CC2G_009974<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_009974")
mean.CC2G_013873<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_013873")
mean.CC2G_007150<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_007150")
mean.CC2G_007168<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_007168")
mean.CC2G_008645<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_008645")
mean.CC2G_009607<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_009607")
mean.CC2G_012681<-subset(summary.expression, Stage == "M" & Strain == "326wt" & Gene == "CC2G_012681")

for (i in (1:nrow(rawdata))){
  if (rawdata$Gene[i]== "CC2G_008991"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_008991$mean
  } else if (rawdata$Gene[i]== "CC2G_002608"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_002608$mean
  } else if (rawdata$Gene[i]== "CC2G_005966"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_005966$mean
  } else if (rawdata$Gene[i]== "CC2G_009974"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_009974$mean
  } else if (rawdata$Gene[i]== "CC2G_013873"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_013873$mean
  } else if (rawdata$Gene[i]== "CC2G_007150"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_007150$mean
  } else if (rawdata$Gene[i]== "CC2G_007168"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_007168$mean
  } else if (rawdata$Gene[i]== "CC2G_008645"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_008645$mean
  } else if (rawdata$Gene[i]== "CC2G_009607"){
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_009607$mean
  } else {
    rawdata$Norm[i]<- rawdata$Expression[i]/mean.CC2G_012681$mean
  } 
}

summary.expression.norm<-ddply(rawdata, c("Stage", "Gene", "Strain"), 
                               summarise, 
                               mean = mean(Norm),
                               std = sd(Norm),
                               se = sd(Norm)/(3^0.5),
                               cv = sd(Norm)/mean(Norm)
)
View(summary.expression.norm)
summary.expression.norm$ysig<-summary.expression.norm$mean+summary.expression.norm$se+0.2
rawdata$group<-paste(rawdata$Stage, rawdata$Strain, sep = "")

####CC2G_008991
raw.CC2G_008991<-subset(rawdata, Gene == "CC2G_008991")

bartlett.CC2G_008991<-bartlett.test(Norm ~ group, data = raw.CC2G_008991)
bartlett.CC2G_008991
aov.CC2G_008991<-aov(Norm ~ group, data = raw.CC2G_008991)
summary(aov.CC2G_008991)
TukeyHSD.CC2G_008991<-TukeyHSD(aov.CC2G_008991)
HSD.CC2G_008991<-HSD.test(aov.CC2G_008991, "group")
HSD.CC2G_008991

####CC2G_002608
raw.CC2G_002608<-subset(rawdata, Gene == "CC2G_002608")

bartlett.CC2G_002608<-bartlett.test(Norm ~ group, data = raw.CC2G_002608)
bartlett.CC2G_002608
aov.CC2G_002608<-aov(Norm ~ group, data = raw.CC2G_002608)
summary(aov.CC2G_002608)
TukeyHSD.CC2G_002608<-TukeyHSD(aov.CC2G_002608)
HSD.CC2G_002608<-HSD.test(aov.CC2G_002608, "group")
HSD.CC2G_002608

####CC2G_005966
raw.CC2G_005966<-subset(rawdata, Gene == "CC2G_005966")

bartlett.CC2G_005966<-bartlett.test(Norm ~ group, data = raw.CC2G_005966)
bartlett.CC2G_005966
aov.CC2G_005966<-aov(Norm ~ group, data = raw.CC2G_005966)
summary(aov.CC2G_005966)
TukeyHSD.CC2G_005966<-TukeyHSD(aov.CC2G_005966)
HSD.CC2G_005966<-HSD.test(aov.CC2G_005966, "group")
HSD.CC2G_005966

####CC2G_009974
raw.CC2G_009974<-subset(rawdata, Gene == "CC2G_009974")

bartlett.CC2G_009974<-bartlett.test(Norm ~ group, data = raw.CC2G_009974)
bartlett.CC2G_009974
aov.CC2G_009974<-aov(Norm ~ group, data = raw.CC2G_009974)
summary(aov.CC2G_009974)
TukeyHSD.CC2G_009974<-TukeyHSD(aov.CC2G_009974)
HSD.CC2G_009974<-HSD.test(aov.CC2G_009974, "group")
HSD.CC2G_009974

####CC2G_013873
raw.CC2G_013873<-subset(rawdata, Gene == "CC2G_013873")

bartlett.CC2G_013873<-bartlett.test(Norm ~ group, data = raw.CC2G_013873)
bartlett.CC2G_013873
aov.CC2G_013873<-aov(Norm ~ group, data = raw.CC2G_013873)
summary(aov.CC2G_013873)
TukeyHSD.CC2G_013873<-TukeyHSD(aov.CC2G_013873)
HSD.CC2G_013873<-HSD.test(aov.CC2G_013873, "group")
HSD.CC2G_013873

####CC2G_007150
raw.CC2G_007150<-subset(rawdata, Gene == "CC2G_007150")

bartlett.CC2G_007150<-bartlett.test(Norm ~ group, data = raw.CC2G_007150)
bartlett.CC2G_007150
aov.CC2G_007150<-aov(Norm ~ group, data = raw.CC2G_007150)
summary(aov.CC2G_007150)
TukeyHSD.CC2G_007150<-TukeyHSD(aov.CC2G_007150)
HSD.CC2G_007150<-HSD.test(aov.CC2G_007150, "group")
HSD.CC2G_007150

####CC2G_007168
raw.CC2G_007168<-subset(rawdata, Gene == "CC2G_007168")

bartlett.CC2G_007168<-bartlett.test(Norm ~ group, data = raw.CC2G_007168)
bartlett.CC2G_007168
aov.CC2G_007168<-aov(Norm ~ group, data = raw.CC2G_007168)
summary(aov.CC2G_007168)
TukeyHSD.CC2G_007168<-TukeyHSD(aov.CC2G_007168)
HSD.CC2G_007168<-HSD.test(aov.CC2G_007168, "group")
HSD.CC2G_007168

####CC2G_008645
raw.CC2G_008645<-subset(rawdata, Gene == "CC2G_008645")

bartlett.CC2G_008645<-bartlett.test(Norm ~ group, data = raw.CC2G_008645)
bartlett.CC2G_008645
aov.CC2G_008645<-aov(Norm ~ group, data = raw.CC2G_008645)
summary(aov.CC2G_008645)
TukeyHSD.CC2G_008645<-TukeyHSD(aov.CC2G_008645)
HSD.CC2G_008645<-HSD.test(aov.CC2G_008645, "group")
HSD.CC2G_008645

####CC2G_009607
raw.CC2G_009607<-subset(rawdata, Gene == "CC2G_009607")

bartlett.CC2G_009607<-bartlett.test(Norm ~ group, data = raw.CC2G_009607)
bartlett.CC2G_009607
aov.CC2G_009607<-aov(Norm ~ group, data = raw.CC2G_009607)
summary(aov.CC2G_009607)
TukeyHSD.CC2G_009607<-TukeyHSD(aov.CC2G_009607)
HSD.CC2G_009607<-HSD.test(aov.CC2G_009607, "group")
HSD.CC2G_009607

####CC2G_012681
raw.CC2G_012681<-subset(rawdata, Gene == "CC2G_012681")

bartlett.CC2G_012681<-bartlett.test(Norm ~ group, data = raw.CC2G_012681)
bartlett.CC2G_012681
aov.CC2G_012681<-aov(Norm ~ group, data = raw.CC2G_012681)
summary(aov.CC2G_012681)
TukeyHSD.CC2G_012681<-TukeyHSD(aov.CC2G_012681)
HSD.CC2G_012681<-HSD.test(aov.CC2G_012681, "group")
HSD.CC2G_012681
