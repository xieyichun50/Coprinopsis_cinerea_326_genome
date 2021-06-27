library("tidyr", lib.loc="/usr/local/lib/R/site-library")
library("ggplot2", lib.loc="/usr/local/lib/R/site-library")
library("scales", lib.loc="/usr/local/lib/R/site-library")
####Import sniffle result
mutant.SV<-read.delim(file = "/home/yichun/kairuku/Nanopore_seq/SV/filelist", header = FALSE)
mutant.SV$V1<-as.character(mutant.SV$V1)
header<-read.delim(file = "/home/yichun/kairuku/Nanopore_seq/SV/header", header = FALSE)
sniffle_result_merge<-NA
for (i in 1:nrow(mutant.SV)) {
  sniffle_result<-read.delim(mutant.SV$V1[i], header = FALSE)
  names(sniffle_result)[1]="CHROM"
  names(sniffle_result)[2]="POS"
  names(sniffle_result)[3]="ID"
  names(sniffle_result)[4]="REF"
  names(sniffle_result)[5]="ALT"
  names(sniffle_result)[6]="QUAL"
  names(sniffle_result)[7]="FILTER"
  names(sniffle_result)[8]="INFO"
  names(sniffle_result)[9]="FORMAT"
  names(sniffle_result)[10]="OTHER"
  sniffle_result$sample=mutant.SV$V1[i]
  sniffle_result_clean <- sniffle_result
  sniffle_result_clean$INFO <- gsub("Snifflesv1.0.11;STD_quant_start", "Snifflesv1.0.11;CHR2=NA;END=NA;STD_quant_start",
                                    sniffle_result_clean$INFO)
  sniffle_result_clean <- separate(sniffle_result_clean, "INFO", c("Precision", "SVMETHOD", "CHR2", "END", 
                                                                   "STD_quant_start", "STD_quant_stop", 
                                                                   "Kurtosis_quant_start", "Kurtosis_quant_stop",
                                                                   "SVTYPE", "SUPTYPE", "SVLEN", "STRANDS", "RE"),
                                   sep = ";", remove = FALSE, convert = TRUE)
  sniffle_result_clean$SVTYPE <- gsub("SVTYPE=", "", sniffle_result_clean$SVTYPE)
  sniffle_result_clean$SUPTYPE <- gsub("SUPTYPE=", "", sniffle_result_clean$SUPTYPE)
  sniffle_result_clean$SVLEN <- gsub("SVLEN=", "", sniffle_result_clean$SVLEN)
  sniffle_result_clean$STRANDS <- gsub("STRANDS=", "", sniffle_result_clean$STRANDS)
  sniffle_result_clean$RE <- gsub("RE=", "", sniffle_result_clean$RE)
  sniffle_result_clean <- subset(sniffle_result_clean, 
                                select = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                                           "Precision", "CHR2", "END",
                                           "STD_quant_start", "STD_quant_stop", 
                                           "Kurtosis_quant_start", "Kurtosis_quant_stop",
                                           "SVTYPE", "SUPTYPE", "SVLEN", "STRANDS", "RE",
                                           "FORMAT", "OTHER", "sample"))
  sniffle_result_merge <- rbind(sniffle_result_merge, sniffle_result_clean)
}

sniffle_result_filter <- subset(sniffle_result_merge, is.na(CHROM)==FALSE)
sniffle_result_filter <- subset(sniffle_result_filter, is.na(CHR2)==FALSE)
sniffle_result_filter <- sniffle_result_filter[grep(pattern="]scaffold", sniffle_result_filter[["ALT"]], invert = TRUE),]
sniffle_result_filter <- sniffle_result_filter[grep(pattern="\\[scaffold", sniffle_result_filter[["ALT"]], invert = TRUE),]

sniffle_result_filter <- subset(sniffle_result_filter, CHROM != "scaffold_19")

sniffle_result_filter <-sniffle_result_filter[order(sniffle_result_filter$CHROM, sniffle_result_filter$POS, sniffle_result_filter$ALT),]
sniffle_result_filter$site<-paste(sniffle_result_filter$CHROM,":",sniffle_result_filter$POS)

write.csv(sniffle_result_filter, file = "/home/yichun/kairuku/Nanopore_seq/SV/SV_mutant/sniffle_result_merge1.csv")

sniffle_result_filter<- read.csv(file = "/home/yichun/kairuku/Nanopore_seq/SV/SV_mutant/sniffle_result_merge1.filter.csv", 
                                 header = TRUE)

SV.summary<-as.data.frame(xtabs(~sample+SVTYPE, sniffle_result_filter))
SV.summary<-as.data.frame(xtabs(~site+SVTYPE, sniffle_result_filter))
SV.summary<-as.data.frame(xtabs(~RE+sample, sniffle_result_filter))

write.csv(SV.summary, file = "/home/yichun/kairuku/Nanopore_seq/SV/SV_mutant/SV.freq.csv")
SV.summary<-read.delim(file = "/home/yichun/kairuku/Nanopore_seq/SV/SV_mutant/RE_SV_count.txt", header = TRUE)

RE.SV.count.plot<-ggplot(SV.summary, mapping = aes(x = RE, y = Sum, color = factor(sample, 
                                                                                  levels = c("UV554","UV506","UV636","UV423",
                                                                                             "UV663","UV433","UV632","UV512",
                                                                                             "UV562","UV616","UV634","UV657",
                                                                                             "UV645"))))+
  geom_line(size = 0.5)+
  #geom_vline(aes(xintercept = 10), color = "gray40", linetype = "dashed", size = 0.3)+
  geom_point(shape = 20, size = 1)+
  labs(x="Read Coverage", y="Number of Structural Variants", colour = "Isolate (Average Coverage)")+
  scale_x_continuous(breaks = seq(0,30,5))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=6),
                labels = trans_format("log10", math_format(10^.x)),
                limits= c(1, 30000))+
  scale_color_manual(values = c("gray0", "gray10", "gray20",
                               "gray30", "gray40", "gray50", "gray60", "gray70", "gray80",
                               "darkorange", "goldenrod", "hotpink", "indianred4"),
                     labels = c("UV554 (25.11 X)","UV506 (23.76 X)","UV636 (22.30 X)","UV423 (18.36 X)",
                                "UV663 (14.76 X)","UV433 (11.97 X)","UV632 (11.65 X)","UV512 (11.63 X)",
                                "UV562 (10.80 X)","UV616   (9.47 X)","UV634   (9.29 X)","UV657   (3.89 X)",
                                "UV645   (2.89 X)"))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.line = element_line(linetype = "solid", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
RE.SV.count.plot
ggsave(filename = "RE.SV.count.plot.png", width = 8, height = 6, units = "in", dpi = 300)
ggsave(filename = "RE.SV.count.plot.tiff", width = 8, height = 6, units = "in", dpi = 300)


SV.summary.plot<-ggplot(data = subset(SV.summary, Freq >=2),
                        mapping = aes(x = paste(site,SVTYPE), y = Freq))+
  geom_bar(stat = "identity", position = "dodge")+
  coord_flip()+
  labs(y = "Number of events", x = "")+
  #scale_y_continuous(limits = c(1.5,13))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        axis.line = element_line(size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3))
SV.summary.plot
ggsave("SV.summary.tiff", width = 6, height = 9, units = "in", dpi = 300)
ggsave("SV.summary.png", width = 6, height = 9, units = "in", dpi = 300)
