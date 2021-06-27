################
##Summarizing SV calling from sniffles
################
library("tidyr", lib.loc="/usr/local/lib/R/site-library")
library("ggplot2", lib.loc="/usr/local/lib/R/site-library")

####Import sniffle result
sniffle_result<-read.delim(file = "cc326_ont_cc130.simple.vcf", header = TRUE)

####split INFO column
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
sniffle_result_clean$CHROM <- gsub("\\.1", "", sniffle_result_clean$CHROM)
sniffle_result_clean$CHROM <- gsub("AACS020000", "Chr_", sniffle_result_clean$CHROM)

sniffle_result_clean<- subset(sniffle_result_clean, CHROM != "AACS02000068.1", 
                              select = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                                         "Precision", "SVMETHOD", "CHR2", "END",
                                         "STD_quant_start", "STD_quant_stop", 
                                         "Kurtosis_quant_start", "Kurtosis_quant_stop",
                                         "SVTYPE", "SUPTYPE", "SVLEN", "STRANDS", "RE",
                                         "FORMAT", "cc326_ont_minimap2_cc130.sorted.bam"))
SV.summary<-as.data.frame(xtabs(~CHROM+SVTYPE, sniffle_result_clean))

SV.summary.plot<-ggplot(data = SV.summary,
          mapping = aes(x = CHROM, y = Freq, fill = factor(SVTYPE)))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c("Chr_01", "Chr_02", "Chr_03",
                              "Chr_04", "Chr_05", "Chr_06",
                              "Chr_07", "Chr_08", "Chr_09",
                              "Chr_10", "Chr_11", "Chr_12", "Chr_13"))+
  labs(y = "Number of events", x = "Chromosome", fill = "SVTYPE")+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        axis.line = element_line(size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3))
SV.summary.plot
ggsave("SV.summary.tiff", width = 6, height = 4, units = "in", dpi = 300)
ggsave("SV.summary.png", width = 6, height = 4, units = "in", dpi = 300)
