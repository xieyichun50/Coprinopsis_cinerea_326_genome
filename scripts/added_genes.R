####Look for genes in the added regions
addgenes<-NA
for (i in 1:nrow(paf.gapped)) {
  alist<-subset(list.gene, region == paf.gapped$Tseq[i] & start >= paf.gapped$gapL[i] & stop <= paf.gapped$gapR[i])
  addgenes<-rbind(addgenes, alist)
}
addgenes<-subset(addgenes, is.na(region) == FALSE)
sum(addgenes$length)

####karyotype
library(RIdeogram, lib.loc = "/usr/local/lib/R/site-library")
karyotype<-read.delim("/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/chrlist.txt", header = TRUE)
gapdense<-subset(paf.gapped, select = c("Tseq", "gapL", "gapR"))
names(gapdense)[1] = "Chr"
names(gapdense)[2] = "Start"
names(gapdense)[3] = "End"
gapdense$Value <- -1

addgenes.dense<-subset(addgenes, select = c("region", "start", "stop"))
names(addgenes.dense)[1] = "Chr"
names(addgenes.dense)[2] = "Start"
names(addgenes.dense)[3] = "End"
addgenes.dense$Value <- 1

gapdense<-rbind(gapdense, addgenes.dense)

ideogram(karyotype = karyotype, overlaid = gapdense)
convertSVG("chromosome.svg", device = "png", width = 32, height = 32, dpi =300)
convertSVG("chromosome.svg", device = "tiff", width = 32, height = 32, dpi =300)
