####CC2G_008991
plotdata.CC2G_008991<-subset(summary.expression.norm, Gene == "CC2G_008991")
plot.CC2G_008991<-ggplot(data = plotdata.CC2G_008991, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "Glycogen synthase\n(CC2G_008991)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 6, 1.5))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
#  annotate(geom = "text", label = "c", x = 0.88, y = 1.27, size = 5, colour = "black")+
#  annotate(geom = "text", label = "c", x = 1.12, y = 0.92, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 1.88, y = 2.57, size = 5, colour = "black")+
#  annotate(geom = "text", label = "c", x = 2.12, y = 1.61, size = 5, colour = "black")+
#  annotate(geom = "text", label = "ab", x = 2.88, y = 5.54, size = 5, colour = "black")+
#  annotate(geom = "text", label = "c", x = 3.12, y = 1.47, size = 5, colour = "black")+
#  annotate(geom = "text", label = "a", x = 3.88, y = 6.22, size = 5, colour = "black")+
#  annotate(geom = "text", label = "c", x = 4.12, y = 1.37, size = 5, colour = "black")+
#  annotate(geom = "text", label = "c", x = 4.88, y = 1.70, size = 5, colour = "black")+
#  annotate(geom = "text", label = "c", x = 5.12, y = 1.55, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 5.88, y = 1.93, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 6.12, y = 2.11, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 6.88, y = 2.20, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 7.12, y = 2.49, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_008991
ggsave("CC2G_008991.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_008991.tiff", width = 6, height = 4, units = "in", dpi = 300)

####CC2G_002608
plotdata.CC2G_002608<-subset(summary.expression.norm, Gene == "CC2G_002608")
plot.CC2G_002608<-ggplot(data = plotdata.CC2G_002608, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "Glycogen phosphorylase\n(CC2G_002608)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 4, 1))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
#  annotate(geom = "text", label = "d", x = 0.88, y = 1.28, size = 5, colour = "black")+
#  annotate(geom = "text", label = "d", x = 1.12, y = 1.23, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abcd", x = 1.88, y = 2.28, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 2.12, y = 2.01, size = 5, colour = "black")+
#  annotate(geom = "text", label = "a", x = 2.88, y = 4.03, size = 5, colour = "black")+
#  annotate(geom = "text", label = "a", x = 3.12, y = 4.02, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abc", x = 3.88, y = 3.17, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abcd", x = 4.12, y = 2.32, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 4.88, y = 1.79, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abcd", x = 5.12, y = 2.61, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abcd", x = 5.88, y = 2.60, size = 5, colour = "black")+
#  annotate(geom = "text", label = "ab", x = 6.12, y = 3.30, size = 5, colour = "black")+
#  annotate(geom = "text", label = "cd", x = 6.88, y = 1.70, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abc", x = 7.12, y = 3.28, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_002608

ggsave("CC2G_002608.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_002608.tiff", width = 6, height = 4, units = "in", dpi = 300)


####CC2G_005966
plotdata.CC2G_005966<-subset(summary.expression.norm, Gene == "CC2G_005966")
plot.CC2G_005966<-ggplot(data = plotdata.CC2G_005966, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "CMGC/MAPK protein kinase\n(CC2G_005966)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 4, 1))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
#  annotate(geom = "text", label = "bcd", x = 0.88, y = 1.22, size = 5, colour = "black")+
#  annotate(geom = "text", label = "d", x = 1.12, y = 0.70, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 1.88, y = 1.48, size = 5, colour = "black")+
#  annotate(geom = "text", label = "cd", x = 2.12, y = 1.09, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 2.88, y = 2.16, size = 5, colour = "black")+
#  annotate(geom = "text", label = "abc", x = 3.12, y = 2.81, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 3.88, y = 1.76, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 4.12, y = 1.98, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 4.88, y = 1.38, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 5.12, y = 2.47, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bc", x = 5.88, y = 2.30, size = 5, colour = "black")+
#  annotate(geom = "text", label = "a", x = 6.12, y = 4.06, size = 5, colour = "black")+
#  annotate(geom = "text", label = "bcd", x = 6.88, y = 1.32, size = 5, colour = "black")+
#  annotate(geom = "text", label = "ab", x = 7.12, y = 2.65, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_005966
ggsave("CC2G_005966.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_005966.tiff", width = 6, height = 4, units = "in", dpi = 300)


####CC2G_009974
plotdata.CC2G_009974<-subset(summary.expression.norm, Gene == "CC2G_009974")
plot.CC2G_009974<-ggplot(data = plotdata.CC2G_009974, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "RasGAP\n(CC2G_009974)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 2, 0.5))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_009974
ggsave("CC2G_009974.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_009974.tiff", width = 6, height = 4, units = "in", dpi = 300)


####CC2G_007168
plotdata.CC2G_007168<-subset(summary.expression.norm, Gene == "CC2G_007168")
plot.CC2G_007168<-ggplot(data = plotdata.CC2G_007168, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "alpha,alpha-trehalose-phosphate synthase\n(CC2G_007168)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 2.5, 0.5))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_007168
ggsave("CC2G_007168.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_007168.tiff", width = 6, height = 4, units = "in", dpi = 300)


####HSD.CC2G_007150
plotdata.CC2G_007150<-subset(summary.expression.norm, Gene == "CC2G_007150")
plot.CC2G_007150<-ggplot(data = plotdata.CC2G_007150, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "RasGAP\n(CC2G_007150)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 1.25, 0.25))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  #geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_007150
ggsave("CC2G_007150.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_007150.tiff", width = 6, height = 4, units = "in", dpi = 300)


####CC2G_008645
plotdata.CC2G_008645<-subset(summary.expression.norm, Gene == "CC2G_008645")
plot.CC2G_008645<-ggplot(data = plotdata.CC2G_008645, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "AGC/PKA protein kinase\n(CC2G_008645)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  #geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  scale_y_continuous(breaks = seq(0, 1.5, 0.5))+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_008645
ggsave("CC2G_008645.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_008645.tiff", width = 6, height = 4, units = "in", dpi = 300)


####CC2G_009607
plotdata.CC2G_009607<-subset(summary.expression.norm, Gene == "CC2G_009607")
plot.CC2G_009607<-ggplot(data = plotdata.CC2G_009607, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "Adenylate cyclase\n(CC2G_009607)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 2, 0.5))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_009607
ggsave("CC2G_009607.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_009607.tiff", width = 6, height = 4, units = "in", dpi = 300)


####CC2G_012681
plotdata.CC2G_012681<-subset(summary.expression.norm, Gene == "CC2G_012681")
plot.CC2G_012681<-ggplot(data = plotdata.CC2G_012681, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"))+
  labs(title = "Trahalase\n(CC2G_012681)", y = "", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 4, 1))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 0), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = "none")
plot.CC2G_012681
ggsave("CC2G_012681.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_012681.tiff", width = 6, height = 4, units = "in", dpi = 300)


plotdata.CC2G_009974<-subset(summary.expression.norm, Gene == "CC2G_009974")
plot.CC2G_009974<-ggplot(data = plotdata.CC2G_009974, 
                         aes(x = Strain, y = mean, 
                             fill = factor(Stage, levels = c("M", "K"))))+
  geom_bar(stat = "identity",  position = "dodge", colour = "black", width = 0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+se, width = 0.15), position = position_dodge(0.5))+
  scale_fill_manual(limits = c("M", "K"), values = c("grey50", "white"), 
                    labels = c("Vegetative mycelia (4 days in darkness)", "Mycelia with hyphal knot (2 light-dark cycles)"))+
  labs(title = NULL, y = "Relative expression", x = NULL, colour = NULL)+
  scale_x_discrete(limits = c("326wt","UV423","UV506","UV512", "UV433", "UV554", "UV663"))+
  scale_y_continuous(breaks = seq(0, 2, 0.5))+
  #geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
  #geom_hline(yintercept = 2, color = "black", linetype = "dashed")+
  annotate(geom = "text", label = "", x = 0.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 1.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 2.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 3.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 4.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 5.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.12, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 6.88, y = 1, size = 5, colour = "black")+
  annotate(geom = "text", label = "", x = 7.12, y = 1, size = 5, colour = "black")+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks.y = element_line(colour = "black", size = 0.5), 
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 0), 
        panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = c(0.36, 0.88))
plot.CC2G_009974
ggsave("CC2G_009974s.png", width = 6, height = 4, units = "in", dpi = 300)
ggsave("CC2G_009974s.tiff", width = 6, height = 4, units = "in", dpi = 300)
