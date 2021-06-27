##########################################
####Summarizing snpEff results
##########################################
Anno.CC2G <- read.delim(file = "/home/yichun/kairuku/cc326_Hybrid_assembly/funannotate/fun_202004/annotate_results/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.annotations.txt",
                        header = TRUE)

Anno.CC2G <- subset(Anno.CC2G, select = c("GeneID", "Feature", "Contig", "Start","Stop", "Strand",
                                          "Name", "Product", "PFAM","EggNog", "COG", "GO.Terms", 
                                          "Secreted", "Membrane", "Protease", "CAZyme"))
Anno.CC2G <-Anno.CC2G[!duplicated(Anno.CC2G$GeneID),]

####4 samples
Anno.t3c5r70s1 <- read.delim(file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.t3.c5r70s1/bcftools_result_merge.t3.c5r70s1.txt", 
                             header = TRUE)

Anno.t3c5r70s1 <- merge(Anno.t3c5r70s1, Anno.CC2G, by = c("GeneID"), all.x = TRUE)

Summary.type <- as.data.frame(xtabs(~Type, Anno.t3c5r70s1))
Summary.type

test<-subset(bcftools_result_merge.t3.c5r70, select = c("CHROM", "POS", "sample"))
names(test)[1] = "Region"
names(test)[2] = "Position"
Anno.t3c5r70s1 <- merge(Anno.t3c5r70s1, test, by = c("Region", "Position"), all.x = TRUE)
write.csv(Anno.t3c5r70s1, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.t3.c5r70s1/Anno.t3c5r70s1.annosample.csv")

Anno.t3c5r70s1.mean <- subset(Anno.t3c5r70s1, Type != "synonymous_variant" & Type !="upstream_gene_variant" & Type != "downstream_gene_variant" & Type != "intergenic_region")
write.csv(Anno.t3c5r70s1.mean, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.t3.c5r70s1/Anno.t3c5r70s1.mean.csv")

Summary.GeneID <- as.data.frame(xtabs(~GeneID + sample, Anno.t3c5r70s1.mean, drop.unused.levels = TRUE))
Summary.GeneID <- subset(Summary.GeneID, Freq > 0)

Summary.GeneIDsamples <- as.data.frame(xtabs(~GeneID, Summary.GeneID))
xtabs(~Freq, Summary.GeneIDsamples)


####13 samples
Anno.t3c5r70s1 <- read.delim(file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r70s1/bcftools_result_merge.c5r70s1.anno.txt", 
                             header = TRUE)

Anno.t3c5r70s1 <- merge(Anno.t3c5r70s1, Anno.CC2G, by = c("GeneID"), all.x = TRUE)

Summary.type <- as.data.frame(xtabs(~Type, Anno.t3c5r70s1))
Summary.type

test<-subset(bcftools_result_merge.c5r70, select = c("CHROM", "POS", "sample"))
names(test)[1] = "Region"
names(test)[2] = "Position"
Anno.t3c5r70s1 <- merge(Anno.t3c5r70s1, test, by = c("Region", "Position"), all.x = TRUE)
write.csv(Anno.t3c5r70s1, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r70s1/Anno.t3c5r70s1.annosample.csv")

Anno.t3c5r70s1.mean <- subset(Anno.t3c5r70s1, Type != "synonymous_variant" & Type !="upstream_gene_variant" & Type != "downstream_gene_variant" & Type != "intergenic_region")
write.csv(Anno.t3c5r70s1.mean, file = "/home/yichun/kairuku/Nanopore_seq/SNP/bcftools_result_merge.c5r70s1/Anno.t3c5r70s1.mean.csv")

Summary.GeneID <- as.data.frame(xtabs(~GeneID + sample, Anno.t3c5r70s1.mean, drop.unused.levels = TRUE))
Summary.GeneID <- subset(Summary.GeneID, Freq > 0)

Summary.GeneIDsamples <- as.data.frame(xtabs(~GeneID, Summary.GeneID))
xtabs(~Freq, Summary.GeneIDsamples)
