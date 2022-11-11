paf<-read.delim("/home/yichun/cc326_Hybrid_assembly/genome_align/cc326_hm/cc326m2h.paf", header = FALSE)
names(paf)[1]="Qseq"
names(paf)[2]="QseqL"
names(paf)[3]="Qstart"
names(paf)[4]="Qend"
names(paf)[5]="strand"
names(paf)[6]="Tseq"
names(paf)[7]="TseqL"
names(paf)[8]="Tstart"
names(paf)[9]="Tend"
names(paf)[10]="Nmatch"
names(paf)[11]="Nbase"
names(paf)[12]="MapQ"
names(paf)[13]="tp"
names(paf)[14]="cm"
names(paf)[15]="s1"
names(paf)[16]="s2"
names(paf)[17]="dv"
names(paf)[18]="rl"
paf.sorted<-paf[order(paf$Tseq, paf$Tstart, paf$Tend),]
rownames(paf.sorted)<-1:nrow(paf.sorted)
paf.sorted$gapL<-NA
paf.sorted$gapR<-NA
add.bind<-NA
for (i in 2:nrow(paf.sorted)) {
  if (paf.sorted$Tseq[i]==paf.sorted$Tseq[i-1]) {
    paf.sorted$gapL[i]<-paf.sorted$Tend[i-1]+1
    paf.sorted$gapR[i]<-paf.sorted$Tstart[i]-1
  } else {
    add.sorted<-paf.sorted[i-1,]
    add.bind<-rbind(add.bind, add.sorted)
  }
}

add.sorted<-paf.sorted[nrow(paf.sorted),]
add.bind<-rbind(add.bind, add.sorted)
add.bind$gapL<-add.bind$Tend+1
add.bind$gapR<-add.bind$TseqL
add.bind<-subset(add.bind, is.na(Tseq) == FALSE)

for (i in 1:nrow(paf.sorted)) {
  if (is.na(paf.sorted$gapL[i]==TRUE)) {
    paf.sorted$gapL[i]<-1
    paf.sorted$gapR[i]<-paf.sorted$Tstart[i]-1
  } else {

  }
}

paf.sorted<-rbind(paf.sorted, add.bind)

paf.sorted$gap<-paste(paf.sorted$Tseq,":",paf.sorted$gapL,"-",paf.sorted$gapR)
paf.gapped<-subset(paf.sorted, paf.sorted$gapL-paf.sorted$gapR < 0)
paf.gapped<-paf[order(paf.gapped$Tseq, paf.gapped$Tstart, paf.gapped$Tend),]
rownames(paf.gapped)<-1:nrow(paf.gapped)
write.csv(paf.gapped, file = "paf.gapped.csv")
##In shell
# sed 's/,/\t/g' paf.gapped.csv | cut -f22 | sed 's/"//g' | sed 's/\s//g' > gaplist


paf.term<-subset(paf.gapped, gapL == 1 | gapR == TseqL)

sum(paf.gapped$gapR-paf.gapped$gapL)
sum(paf.term$gapR-paf.term$gapL)

