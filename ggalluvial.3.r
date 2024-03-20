library(sankeyD3)
library(patchwork)
library(webshot)

library(magrittr)

countCellCount = readRDS("countCellCount.RData")

colorlist = hue_pal()(20)
colorlist = rep(brewer.pal(12,"Paired")[1:10],2)
colorlist = c(brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"))
pdf("mb.raw.spatial.fig.ggalluvial.pdf")
ggplot(countCellCount, aes(y = Freq,axis1 = factor(annotation,levels = c("Meninges","Cortex_L1","Cortex-L2/3","Cortex_L4","Cortex_L5","Cortex_L6","Lateral-ventral cortex","Cortical amygdalar area","Posterior amygdalar nucleus","Fiber tract","Subiculum","Stratum oriens of CA1", "CA1","Stratum lacunosum/raditum of CA1","Molecular layer of dentate gyrus","Dentate gyrus","Cavity","Thalamus","Midbrain","Substandia nigra/Ventral tegmental area")), axis2 = factor(celltype,levels = c("VLMC2","TEGLU7","TEGLU4","TEGLU8","TEGLU10","TEGLU2","TEGLU3","TEGLU12","TEGLU15","TEGLU22","TEGLU6","TEGLU23","TEGLU24","ACTE2","DGGRC2","DEGLU1","MOL2","ACNT1","TEINH18")), fill = annotation, alpha=Freq))+
  geom_flow()+
  geom_alluvium(aes(fill = annotation), curve_type = "sine") +
  # scale_fill_manual(values = c(Brown = "#70493D", Hazel = "#E2AC76",Green = "#3F752B", Blue  = "#81B0E4")) +
  scale_fill_manual(values = colorlist) +
  guides(fill = FALSE) +
  geom_stratum(alpha = .2) +
  geom_text(stat = "stratum", size=3, aes(label = after_stat(stratum)), reverse = T) +
  scale_x_continuous(breaks = 1:2, expand = c(0,0), labels = c("annotation", "celltype")) +
  ggtitle("") +
  xlab("") + ylab("") +
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) 

ggplot(countCellCount, aes(y = celltyperatio,axis1 = factor(annotation,levels = c("Meninges","Cortex_L1","Cortex-L2/3","Cortex_L4","Cortex_L5","Cortex_L6","Lateral-ventral cortex","Cortical amygdalar area","Posterior amygdalar nucleus","Fiber tract","Subiculum","Stratum oriens of CA1", "CA1","Stratum lacunosum/raditum of CA1","Molecular layer of dentate gyrus","Dentate gyrus","Cavity","Thalamus","Midbrain","Substandia nigra/Ventral tegmental area")), axis2 = factor(celltype,levels = c("VLMC2","TEGLU7","TEGLU4","TEGLU8","TEGLU10","TEGLU2","TEGLU3","TEGLU12","TEGLU15","TEGLU22","TEGLU6","TEGLU23","TEGLU24","ACTE2","DGGRC2","DEGLU1","MOL2","ACNT1","TEINH18")), fill = annotation, alpha=celltyperatio))+
  geom_flow()+
  geom_alluvium(aes(fill = annotation), curve_type = "sine") +
  # scale_fill_manual(values = c(Brown = "#70493D", Hazel = "#E2AC76",Green = "#3F752B", Blue  = "#81B0E4")) +
  scale_fill_manual(values = colorlist) +
  guides(fill = FALSE) +
  geom_stratum(alpha = .2) +
  geom_text(stat = "stratum", size=3, aes(label = after_stat(stratum)), reverse = T) +
  scale_x_continuous(breaks = 1:2, expand = c(0,0), labels = c("annotation", "celltype")) +
  ggtitle("") +
  xlab("") + ylab("") +
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) 

ggplot(countCellCount, aes(y = annotationratio,axis1 = factor(annotation,levels = c("Meninges","Cortex_L1","Cortex-L2/3","Cortex_L4","Cortex_L5","Cortex_L6","Lateral-ventral cortex","Cortical amygdalar area","Posterior amygdalar nucleus","Fiber tract","Subiculum","Stratum oriens of CA1", "CA1","Stratum lacunosum/raditum of CA1","Molecular layer of dentate gyrus","Dentate gyrus","Cavity","Thalamus","Midbrain","Substandia nigra/Ventral tegmental area")), axis2 = factor(celltype,levels = c("VLMC2","TEGLU7","TEGLU4","TEGLU8","TEGLU10","TEGLU2","TEGLU3","TEGLU12","TEGLU15","TEGLU22","TEGLU6","TEGLU23","TEGLU24","ACTE2","DGGRC2","DEGLU1","MOL2","ACNT1","TEINH18")), fill = annotation, alpha=annotationratio))+
  geom_flow()+
  geom_alluvium(aes(fill = annotation), curve_type = "sine") +
  # scale_fill_manual(values = c(Brown = "#70493D", Hazel = "#E2AC76",Green = "#3F752B", Blue  = "#81B0E4")) +
  scale_fill_manual(values = colorlist) +
  guides(fill = FALSE) +
  geom_stratum(alpha = .2) +
  geom_text(stat = "stratum", size=3, aes(label = after_stat(stratum)), reverse = T) +
  scale_x_continuous(breaks = 1:2, expand = c(0,0), labels = c("annotation", "celltype")) +
  ggtitle("") +
  xlab("") + ylab("") +
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) 
dev.off()
