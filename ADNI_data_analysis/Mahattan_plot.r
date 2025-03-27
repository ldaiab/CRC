cd output
module load r/4.4.0;R
###################################################################################################
Coxph=read.csv('LOG_CPH.csv')[,-1]
IPCW=read.csv('LOG_IPCW.csv')[,-1]
KLR=read.csv('LOG_KLR.csv')[,-1]
N01=read.csv('LOG_N01.csv')[,-1]
U01=read.csv('LOG_U1.csv')[,-1]
Coxph_NoLOG=read.csv('NoLOG_CPH.csv')[,-1]
IPCW_NoLOG=read.csv('NoLOG_IPCW.csv')[,-1]
KLR_NoLOG=read.csv('NoLOG_KLR.csv')[,-1]
M=as.data.frame(cbind(Coxph[,2],IPCW[,2],KLR[,2],N01[,2],U01[,2],Coxph_NoLOG[,2],IPCW_NoLOG[,2],KLR_NoLOG[,2]))
colnames(M)=c('Coxph','IPCW','KLR','N01','U01','Coxph_NoLOG','IPCW_NoLOG','KLR_NoLOG')
M1=read.csv('../../ADNI_Gene/ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv');
temp=data.frame(list(gene=M1[,6],count=rep(0,length(M1[,6]))))
UG=unique(M1[,6]);L=length(UG)
temp3=duplicated(M1[,6])
for(ii in 1:L)
{
temp2=which(M1[,6]==UG[ii])
temp[temp2,2]=cumsum(temp3[temp2])
}
temp1=paste0(M1[,6],'_',temp[,2])
rownames(M)=temp1
FDR=M
p=dim(FDR)[2]
for(ii in 1:p)
{
FDR[,ii]=p.adjust(M[,ii],'fdr')
}
#
colnames(FDR)=c("CPH", "IPCW", "KLR", "CRC[N(0,1)]", "CRC[U(-1,1)]", "CPH (No-Log Trans)","IPCW (No-Log Trans)","KLR (No-Log Trans)")
FDR=FDR[,c(1,6,2,7,3,8,4,5)]
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
data_plot <- FDR %>%
  mutate(Position = row_number(),
         Protein = rownames(FDR)) %>%  # Capture protein names
  pivot_longer(
    cols = -c(Position, Protein),
    names_to = "Method",
    values_to = "FDR"
  ) %>%
  mutate(
    Chromosome = as.numeric(factor(Method,levels =  unique(Method))),
    neg_log_FDR = -log10(FDR)
  ) %>%
  arrange(Chromosome, Position)
data_plot[,2]=gsub('_0','',as.matrix(data_plot[,2]))
data_cum <- data_plot %>%
  group_by(Chromosome) %>%
  summarise(max_pos = max(Position)) %>%
  mutate(add_pos = lag(cumsum(max_pos), default = 0)  ) %>% #
  select(Chromosome, add_pos) %>%
  left_join(data_plot, ., by = "Chromosome") %>%
  arrange(Chromosome, Position) %>%
  mutate(cum_pos = Position + add_pos)
# Threshold data
thresholds <- data.frame(
  value = c(0.2, 0.1, 0.05),
  color = c("orange", "purple", "red"),
  neg_log_value = -log10(c(0.2, 0.1, 0.05))
)
significant_points <- data_cum %>%
  filter(FDR <= 0.1) %>%
  group_by(Protein) %>%
  slice_max(neg_log_FDR, n = 5) %>%  # Keep highest points per protein
  ungroup()
  
method_centers <- data_cum %>%
  group_by(Chromosome) %>%
  summarise(center = mean(cum_pos)) %>%
  pull(center) 

method_labels <- c(
  "CPH",
  "CPH (No-Log Trans)",
  "IPCW",
  "IPCW (No-Log Trans)",
  "KLR",
  "KLR (No-Log Trans)",
  expression(CRC[N(0,1)]),
  expression(CRC[U(-1,1)])
)

pdf('compare1.pdf', width = 17, height = 8) # Narrower width
ggplot(data_cum, aes(x = cum_pos, y = neg_log_FDR)) +
  geom_point(aes(color = factor(Chromosome)), alpha = 0.7, size = 1) +
  geom_hline(data = thresholds,
             aes(yintercept = neg_log_value, linetype = factor(value)),
             color = "#D3D3D3", linewidth = 0.8) +
  geom_text_repel(
    data = significant_points,
    aes(label = Protein),
    size = 7,
    box.padding = 0.3,
    max.overlaps = 30,
    segment.color = "black",
    min.segment.length = 0.1
  ) +
  scale_color_manual(
    name = "Methods",
    values = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928","#000000","#36013F","#006D5B"), # Colorblind-friendly
    #labels = levels(factor(data_plot$Method))
	labels = method_labels
  ) +
  scale_linetype_manual(
    name = "FDR Thresholds",
    values = c("dotted", "dashed", "longdash"),
    labels = c("0.05", "0.10", "0.20")
  ) +
  scale_x_continuous(
    #label = levels(factor(data_plot$Method)),
	labels = method_labels,
    breaks = method_centers
  ) +
  labs(x = "Statistical Methods", y = "-log10(FDR)") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    legend.key.size = unit(1.5, "lines"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.margin = margin(1, 3, 1, 1, "cm")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)), # Larger legend dots
    linetype = guide_legend(override.aes = list(linewidth = 1))
  )
  dev.off()