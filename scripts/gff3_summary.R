gff3<-read.delim("/home/yichun/cc326_Hybrid_assembly/funannotate/fun_202004/annotate_results/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.clean.gff3",
                 header = FALSE)
names(gff3)[1]="region"
names(gff3)[2]="source"
names(gff3)[3]="type"
names(gff3)[4]="start"
names(gff3)[5]="stop"
names(gff3)[6]="score"
names(gff3)[7]="strand"
names(gff3)[8]="phase"
names(gff3)[9]="attributes"

gff3$length <- gff3$stop-gff3$start+1

gff3.summary<-matrix(data = NA, nrow = 6, ncol = 4)
gff3.summary<-as.data.frame(gff3.summary)
names(gff3.summary)[1]="Number"
names(gff3.summary)[2]="N.Average"
names(gff3.summary)[3]="Length"
names(gff3.summary)[4]="L.Average"
row.names(gff3.summary)[1]="gene"
row.names(gff3.summary)[2]="CDS"
row.names(gff3.summary)[3]="Exon"
row.names(gff3.summary)[4]="Intron"
row.names(gff3.summary)[5]="5’-UTR"
row.names(gff3.summary)[6]="3’-UTR"

list.gene<-subset(gff3, type == "gene")
gff3.summary$Number[1]<-nrow(list.gene)
gff3.summary$Length[1]<-sum(list.gene$length)
gff3.summary$L.Average[1]<-gff3.summary$Length[1]/gff3.summary$Number[1]

list.CDS<-subset(gff3, type == "CDS")
gff3.summary$Number[2]<-nrow(list.CDS)
gff3.summary$Length[2]<-sum(list.CDS$length)
gff3.summary$N.Average[2]<-gff3.summary$Number[2]/gff3.summary$Number[1]
gff3.summary$L.Average[2]<-gff3.summary$Length[2]/gff3.summary$Number[1]

list.exon<-subset(gff3, type == "exon")
gff3.summary$Number[3]<-nrow(list.exon)
gff3.summary$Length[3]<-sum(list.exon$length)
gff3.summary$N.Average[3]<-gff3.summary$Number[3]/gff3.summary$Number[1]
gff3.summary$L.Average[3]<-gff3.summary$Length[3]/gff3.summary$Number[1]

list.CDS<-subset(gff3, type == "five_prime_UTR")
gff3.summary$Number[5]<-nrow(list.CDS)
gff3.summary$Length[5]<-sum(list.CDS$length)
gff3.summary$N.Average[5]<-gff3.summary$Number[5]/gff3.summary$Number[1]
gff3.summary$L.Average[5]<-gff3.summary$Length[5]/gff3.summary$Number[1]

list.3UTR<-subset(gff3, type == "three_prime_UTR")
gff3.summary$Number[6]<-nrow(list.3UTR)
gff3.summary$Length[6]<-sum(list.3UTR$length)
gff3.summary$N.Average[6]<-gff3.summary$Number[6]/gff3.summary$Number[1]
gff3.summary$L.Average[6]<-gff3.summary$Length[6]/gff3.summary$Number[1]

gff3.summary$Number[4]<-gff3.summary$Number[3]-gff3.summary$Number[1]
gff3.summary$Length[4]<-gff3.summary$Length[1]-gff3.summary$Length[3]
gff3.summary$N.Average[4]<-gff3.summary$Number[4]/gff3.summary$Number[1]
gff3.summary$L.Average[4]<-gff3.summary$Length[4]/gff3.summary$Number[1]

write.csv(gff3.summary,
          "/home/yichun/cc326_Hybrid_assembly/funannotate/fun_202004/annotate_results/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.clean.gff3stat.csv")

write.table(list.gene,"/home/yichun/cc326_Hybrid_assembly/funannotate/fun_202004/annotate_results/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gene.txt")
