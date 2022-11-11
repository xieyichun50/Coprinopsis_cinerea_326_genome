###################################
####venn diagram of DEGs
###################################
library("VennDiagram", lib.loc="/usr/local/lib/R/site-library")

cc326hm<-read.delim("/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/cc326_hm_genelist", header = FALSE)
cc130g<-read.delim("/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/cc130to326h.common.txt", header = FALSE,na.strings=c("","NA"))
cc326m<-read.delim("/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/cc326m_genelist", header = FALSE, na.strings=c("","NA"))
Mid<-read.delim("/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/muraguchi_PID_CC1G", header = FALSE)

names(cc130g)[1]<-"CC1G"
names(cc130g)[2]<-"CC2G"
names(cc326m)[1]<-"mid"
names(cc326m)[2]<-"CC2G"
names(Mid)[1]<-"mid"
names(Mid)[2]<-"CC1G"

for (i in 1:nrow(cc130g)){
  if (is.na(cc130g$CC2G[i])==FALSE){
    cc130g$combine[i] <- cc130g$CC2G[i]
  } else {
    cc130g$combine[i]=cc130g$CC1G[i]
  }
}

cc326m<-merge(cc326m, Mid, all.x=TRUE, by=c("mid"))
for (i in 1:nrow(cc326m)) {
  if (is.na(cc326m$CC2G[i])==FALSE){
    cc326m$combine[i] <- cc326m$CC2G[i]
  } else if (is.na(cc326m$CC1G[i])==FALSE){
    cc326m$combine[i] <- cc326m$CC1G[i]
  } else {
    cc326m$combine[i] <- cc326m$mid[i]
  }
}


venn.diagram(x = list(Cc130_WGS = cc130g$combine, Cc326_Illumina_only = cc326m$combine, Cc326_Hybrid_assembly = cc326hm$V1),
             filename = "common_gene.png", imagetype = "png" , units = "in", height = 4, width = 4, resolution = 300,
             fill = c("burlywood2", "darkgray", "plum2"), margin = 0.1, cat.dist = 0.1)
venn.diagram(x = list(Cc130_WGS = cc130g$combine, Cc326_Illumina_only = cc326m$combine, Cc326_Hybrid_assembly = cc326hm$V1),
             filename = "common_gene", imagetype = "tiff" , units = "in", height = 4, width = 4, resolution = 300,
             fill = c("burlywood2", "darkgray", "plum2"), margin = 0.1, cat.dist = 0.1)


####ID match
names(cc326hm)[1]="CC2G"
IDmatch3 <- merge(cc326hm, cc326m, by = c("CC2G"), all.x = TRUE)
IDmatch3 <- merge(IDmatch3, cc130g, by = c("CC2G"), all.x = TRUE)
write.table(IDmatch3, file = "/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/IDmatch3.txt")
