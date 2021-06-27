library("tidyr", lib.loc="/usr/local/lib/R/site-library")
library("ggplot2", lib.loc="/usr/local/lib/R/site-library")

####Import bcftools result
mutant.SNP<-read.delim(file = "/home/yichun/kairuku/Nanopore_seq/SNP/filelist", header = FALSE)
mutant.SNP$V1<-as.character(mutant.SNP$V1)

bcftools_result_merge<-NA
for (i in 1:nrow(mutant.SNP)) {
  bcftools_result<-read.delim(mutant.SNP$V1[i], header = FALSE)
  names(bcftools_result)[1]="CHROM"
  names(bcftools_result)[2]="POS"
  names(bcftools_result)[3]="ID"
  names(bcftools_result)[4]="REF"
  names(bcftools_result)[5]="ALT"
  names(bcftools_result)[6]="QUAL"
  names(bcftools_result)[7]="FILTER"
  names(bcftools_result)[8]="INFO"
  names(bcftools_result)[9]="FORMAT"
  names(bcftools_result)[10]="OTHER"
  bcftools_result_clean<-subset(bcftools_result, ALT != "<*>")
  bcftools_result_clean$sample=mutant.SNP$V1[i]
  #bcftools_result_clean$sample<- gsub("/home/yichun/Nanopore_seq/SNP/","",bcftools_result_clean$sample)
  #bcftools_result_clean$sample<- gsub(".snp.filter.vcf","",bcftools_result_clean$sample)
  
  bcftools_result_clean <- separate(bcftools_result_clean, "INFO", c("DP", "I16","QS"),
                                    sep = ";", remove = FALSE, convert = TRUE)
  bcftools_result_clean$DP <- gsub("DP=", "", bcftools_result_clean$DP)
  bcftools_result_clean$I16 <- gsub("I16=", "", bcftools_result_clean$I16)
  bcftools_result_clean$QS <- gsub("QS=", "", bcftools_result_clean$QS)
  bcftools_result_clean <- separate(bcftools_result_clean, "I16", 
                                    c("refpassF", "refpassR","altpassF","altpassR",
                                      "sumrefQ","sumrefQ^2","sumaltQ","sumaltQ^2",
                                      "sumrefmapQ","sumrefmapQ^2","sumaltmapQ","sumaltmapQ^2",
                                      "sumrefdistQ","sumrefdistQ^2","sumaltdistQ","sumaltdistQ^2"),
                                    sep = ",", remove = FALSE, convert = TRUE)
  bcftools_result_clean$coverage.ref<-bcftools_result_clean$refpassF+bcftools_result_clean$refpassR
  bcftools_result_clean$coverage.alt<-bcftools_result_clean$altpassF+bcftools_result_clean$altpassR
  bcftools_result_clean$ratio<-bcftools_result_clean$coverage.alt/(bcftools_result_clean$coverage.alt+bcftools_result_clean$coverage.ref)
  bcftools_result_clean <- subset(bcftools_result_clean, 
                                  select = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                                             "ratio","coverage.ref","coverage.alt",
                                             "refpassF", "refpassR","altpassF","altpassR",
                                             "sumrefQ","sumrefQ^2","sumaltQ","sumaltQ^2",
                                             "sumrefmapQ","sumrefmapQ^2","sumaltmapQ","sumaltmapQ^2",
                                             "sumrefdistQ","sumrefdistQ^2","sumaltdistQ","sumaltdistQ^2",
                                             "FORMAT", "OTHER", "sample"))
  bcftools_result_merge <- rbind(bcftools_result_merge, bcftools_result_clean)
}

rm(bcftools_result, bcftools_result_clean)
bcftools_result_merge <- subset(bcftools_result_merge, is.na(CHROM)==FALSE)
bcftools_result_merge <- bcftools_result_merge[order(bcftools_result_merge$CHROM, bcftools_result_merge$POS, bcftools_result_merge$ALT),]
bcftools_result_merge$site <- paste(bcftools_result_merge$CHROM,":",bcftools_result_merge$POS)
write.csv(bcftools_result_merge, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.csv")

bcftools_result_merge.coverage3 <- subset(bcftools_result_merge, coverage.alt >=3)
write.csv(bcftools_result_merge.coverage3, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.coverage3.csv")
bcftools_result_merge.coverage5 <- subset(bcftools_result_merge, coverage.alt >=5)
bcftools_result_merge.r50 <- subset(bcftools_result_merge, coverage.alt >=3 & ratio >=0.5)
SNP.summary<-as.data.frame(xtabs(~sample, bcftools_result_merge.r50))

bcftools_result_merge.c5r70 <- subset(bcftools_result_merge, coverage.alt >=5 & ratio >=0.70)
SNP.summary<-as.data.frame(xtabs(~sample, bcftools_result_merge.c5r70))
SNP.summary
bcftools_result_merge.c5r70.vcf <- subset(bcftools_result_merge.c5r70, 
                                      select = c("CHROM", "POS", "ID", "REF", "ALT","QUAL","FILTER"))
write.csv(bcftools_result_merge.c5r70.vcf, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r70.vcf")
write.csv(bcftools_result_merge.c5r70, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r70.csv")

SNP.summary<-as.data.frame(xtabs(~site, bcftools_result_merge.c5r70))
SNP.summary
xtabs(~Freq, SNP.summary)


bcftools_result_merge.c5r90 <- subset(bcftools_result_merge.c5r70, coverage.alt >=5 & ratio >=0.90)
SNP.summary<-as.data.frame(xtabs(~sample, bcftools_result_merge.c5r90))
SNP.summary
bcftools_result_merge.c5r90.vcf <- subset(bcftools_result_merge.c5r90, 
                                          select = c("CHROM", "POS", "ID", "REF", "ALT","QUAL","FILTER"))
write.csv(bcftools_result_merge.c5r90, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r90.csv")
write.csv(bcftools_result_merge.c5r90.vcf, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r90.vcf")
SNP.summary<-as.data.frame(xtabs(~site, bcftools_result_merge.c5r90))
SNP.summary
xtabs(~Freq, SNP.summary)



