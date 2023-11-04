library(ggplot2)
library(RColorBrewer)
library(patchwork)
# https://zhuanlan.zhihu.com/p/630484390

```
3753	Astroependymal	P95
416	Excitatory_Neuroblasts	P95
30976	Excitatory_Neurons	P95
631	Inhibitory_Neuroblasts	P95
5885	Inhibitory_Neurons	P95
1326	Microglia	P95
90	OEC	P95
6891	Oligodendrocytes	P95
719	OPC	P95
1163	Vascular	P95
4461	Astroependymal	P11
5507	Excitatory_Neuroblasts	P11
16498	Excitatory_Neurons	P11
7262	Inhibitory_Neuroblasts	P11
7762	Inhibitory_Neurons	P11
1002	Microglia	P11
88	OEC	P11
1102	Oligodendrocytes	P11
2609	OPC	P11
1332	Vascular	P11		
4317	Astroependymal	P2
954	Excitatory_Neuroblasts	P2
16371	Excitatory_Neurons	P2
8818	Inhibitory_Neuroblasts	P2
15072	Inhibitory_Neurons	P2
236	Microglia	P2
74	OEC	P2
1072	Oligodendrocytes	P2
266	OPC	P2
1096	Vascular	P2
```
mydat = read.table("C:\\Users\\dell\\Desktop\\map.stat.txt", header=F)
names(mydat) = c("count","celltype","age")
mydat$age = factor(mydat$age,levels = c("P2","P11","P95"))

p1 = ggplot(data = mydat, mapping = aes(x = celltype, y = count, fill = age)) +
  geom_bar(stat = 'identity', position = 'fill',colour = 'black') +
  scale_fill_manual(values = brewer.pal(length(unique(mydat$age)),"Spectral")) +
  # scale_fill_brewer(palette="Dark2") +
  ylab('Fraction') + 
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) 

p2 = ggplot(data = mydat, mapping = aes(x = age, y = count, fill = celltype)) +
  geom_bar(stat = 'identity', position = 'fill',colour = 'black') +
  scale_fill_manual(values = brewer.pal(length(unique(mydat$celltype)),"Spectral")) +
  # scale_fill_brewer(palette="Dark2") +
  # scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(mydat$celltype)))) +
  ylab('Fraction') + 
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() + 
  theme(panel.grid=element_blank()) 


p1 / p2


