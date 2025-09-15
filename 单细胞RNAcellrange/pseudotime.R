library(magrittr)
library(dplyr)
library(reshape2)
library(ggrepel)
library(ggplot2)
library(patchwork)

setwd("E:/BaiduSyncdisk/Bioinformatic/06.scCHIATAC/01.3Dmodeling/pseudotime")

meta_k = read.table("k562.final.stat.txt", header = T, sep = "\t")
meta_p = read.table("patski.final.stat.txt", header = T, sep = "\t")

count_k = read.table("k562.PET.bedpe.gz.stat", header = T, sep = "\t")
count_p = read.table("patski.PET.bedpe.gz.stat", header = T, sep = "\t")

for (i in c("k", "p")) {
  count = get(paste0("count_", i))
  meta = get(paste0("meta_", i))
  
  # 计算各个group中心点
  data = merge(meta, count, by = "row.names", all.x = T) %>% 
    select(., -intra, -inter, -intragt1000, -intrale1000)
  colnames(data)[c(1, ncol(data)-2, ncol(data)-1, ncol(data))] = c("cellid", "PC_1", "PC_2", "contact_number")
  center = data.frame(S = c(mean(data[data$Phase == "S", ]$PC_1), mean(data[data$Phase == "S", ]$PC_2)), 
                      G1 = c(mean(data[data$Phase == "G1", ]$PC_1), mean(data[data$Phase == "G1", ]$PC_2)), 
                      G2M = c(mean(data[data$Phase == "G2M", ]$PC_1), mean(data[data$Phase == "G2M", ]$PC_2)), 
                      row.names = c("PC_1", "PC_2"))
  
  # 计算点到上一个时期中心点的距离
  G1 = data[data$Phase == "G1", ] %>% 
    mutate(distance = ((.$PC_1 - center["PC_1", "S"])^2 + (.$PC_2 - center["PC_2", "S"])^2) - ((.$PC_1 - center["PC_1", "G2M"])^2 + (.$PC_2 - center["PC_2", "G2M"])^2)) %>% 
    .[order(-.$distance), ]
  S = data[data$Phase == "S", ] %>% 
    mutate(distance = ((.$PC_1 - center["PC_1", "G2M"])^2 + (.$PC_2 - center["PC_2", "G2M"])^2) - ((.$PC_1 - center["PC_1", "G1"])^2 + (.$PC_2 - center["PC_2", "G1"])^2)) %>% 
    .[order(-.$distance), ]
  G2M = data[data$Phase == "G2M", ] %>% 
    mutate(distance = ((.$PC_1 - center["PC_1", "G1"])^2 + (.$PC_2 - center["PC_2", "G1"])^2) - ((.$PC_1 - center["PC_1", "S"])^2 + (.$PC_2 - center["PC_2", "S"])^2)) %>% 
    .[order(-.$distance), ]
  
  data = rbind(G1, S, G2M)
  
  assign(paste0("data_", i), data)
}

rm(i, count, meta, data, G1, S, G2M, center)
write.table(data_k, file = "k562.dis.stat.txt", quote = F, sep = "\t", row.names = F)
write.table(data_p, file = "patski.dis.stat.txt", quote = F, sep = "\t", row.names = F)

#awk -F "\t" -v phase="" -v num=0 -v group=0 'NR > 1 {if(phase != $(NF-4)){phase=$(NF-4);num=$(NF-1);group=0}else{if(num < 1000000){num+=$(NF-1)}else{num=$(NF-1);group+=1}};print $0"\t"phase"_"group}'

for (type in c("k562", "patski")) {
  data = read.table(paste0(type, ".dis.stat.txt"), header = T, sep = "\t")
  colnames(data)[c(ncol(data)-2, ncol(data))] = c("ncontact", "group")
  stat = aggregate(data$ncontact, by = list(data$group), sum)
  data = data[data$group %in% stat[stat$x >= 1000000, ]$Group.1, ]
  meta = aggregate(data[, c("PC_1", "PC_2", "distance")], by = list(data$group), mean) %>% 
    merge(., stat, by = "Group.1", all.x = T)
  colnames(meta) = c("group", "PC_1", "PC_2", "distance", "contact_num")
  # write.table(meta, file = paste0(type, ".meta.stat.txt"), quote = F, sep = "\t", row.names = F)
  # write.table(data, file = paste0(type, ".pseudotime.stat.txt"), quote = F, sep = "\t", row.names = F)
  data = aggregate(data[, c("PC_1", "PC_2")], by = list(data$Phase, data$group), mean) %>% 
    mutate(Group.1 = .$Group.1 %>% factor(., levels = c("G1", "S", "G2M"))) %>% 
    mutate(color = gsub(".*_", "", .$Group.2) %>% as.numeric) %>% 
    .[order(.$Group.1, .$color), ] %>% 
    mutate(color = 1:nrow(.))
  colnames(data) = c("Phase", "group", "PC_1", "PC_2", "color")
  data$group = gsub("^.*_", "", data$group)
  #data$color = data$color / max(data$color)
  
  p = ggplot(data, aes(x = PC_2, y = PC_1)) +
    geom_point(aes(color = color, shape = Phase)) +
    geom_text_repel(aes(label = group)) +
    scale_color_distiller(name = paste0(type, " pseudotime"), palette = "RdYlBu") +
    labs(title = type, x = "PC 2", y = "PC 1") +
    theme_bw() +
    theme(title = element_text(size = 14), 
          legend.text = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 12), 
          axis.text = element_text(color = "black", size =8))
  
  assign(type, p)
}

p = k562 + patski + 
  plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = "A")

pdf(file = "pca_separate.label.pdf", width = 22, height = 10, pointsize = 8)
print(p)
dev.off()

#绘制dis计算示意图
data = read.table("k562.dis.stat.txt", header = T, sep = "\t")
data$group = data$metacell
center = data.frame(S = c(mean(data[data$Phase == "S", ]$PC_1), mean(data[data$Phase == "S", ]$PC_2)), 
                    G1 = c(mean(data[data$Phase == "G1", ]$PC_1), mean(data[data$Phase == "G1", ]$PC_2)), 
                    G2M = c(mean(data[data$Phase == "G2M", ]$PC_1), mean(data[data$Phase == "G2M", ]$PC_2)), 
                    row.names = c("PC_1", "PC_2")) %>% t %>% as.data.frame %>% 
  mutate(Phase = rownames(.)) %>% 
  mutate(group = paste0(.$Phase, " center "))
plotdata = rbind(center, data[1, c("PC_1", "PC_2", "Phase", "group")])
plotdata$Phase = factor(plotdata$Phase, levels = c("G1", "S", "G2M"))
p = ggplot(plotdata, aes(x = PC_2, y = PC_1)) +
  geom_segment(aes(x = plotdata[plotdata$group == "G1_0", ]$PC_2, 
                   xend = plotdata[plotdata$Phase == "S", ]$PC_2, 
                   y = plotdata[plotdata$group == "G1_0", ]$PC_1, 
                   yend = plotdata[plotdata$Phase == "S", ]$PC_1), 
               linetype = 2, linewidth = 1, alpha = 0.1, color = "#00BA38") +
  geom_segment(aes(x = plotdata[plotdata$group == "G1_0", ]$PC_2, 
                   xend = plotdata[plotdata$Phase == "G2M", ]$PC_2, 
                   y = plotdata[plotdata$group == "G1_0", ]$PC_1, 
                   yend = plotdata[plotdata$Phase == "G2M", ]$PC_1), 
               linetype = 2, linewidth = 1, alpha = 0.1, color = "#619CFF") +
  geom_point(aes(color = Phase, shape = Phase), size = 3) +
  geom_text_repel(aes(label = group)) +
  labs(x = "PC 2", y = "PC 1") +
  theme_bw() +
  theme(panel.grid = element_blank())
  
  
pdf(file = "dis.pdf", width = 5, height = 4, pointsize = 8)
print(p)
dev.off()


percentage <- function(x) { return (x / sum(x)) }

##### k562 order by contact #####
finaldat = read.table("k562.finaldat.txt", header = T, sep = "\t")
finaldat$cellid = rownames(finaldat)

data = read.table("k562.dis.stat.txt", header = T, sep = "\t", row.names = 1)
data = data[, c("S.Score", "G2M.Score")] %>%
  mutate(Gene_diff = .$S.Score - .$G2M.Score) %>%
  mutate(Pesudotime_dis = 1:nrow(.) / nrow(.)) %>%
  .[finaldat$cellid, ] %>%
  .[grep("^NA", rownames(.), invert = T), ] %>% 
  mutate(order = 1:nrow(.))

data1 = merge(data, finaldat[, c("cellid", "type")], by.x = "row.names", by.y = "cellid", all.x = T)

# data2 = finaldat[finaldat$cellid %in% rownames(data), ]
# cellid = data2$cellid
# data2 = data2[, grep("^bin", colnames(data2))]  %>%
#   sapply(., percentage) %>%
#   as.data.frame %>% t %>% as.data.frame
# colnames(data2) = cellid
# data2[is.na(data2)] = 0
# data2 = scale(data2, center = F, scale = colSums(data2)) %>%
#   as.data.frame %>%
#   mutate(region = factor(c(rep("20k-2M", 52), rep("2M-12M", 20), rep(">12M", nrow(.)-72)),
#                          levels = c("20k-2M", "2M-12M", ">12M")))
# data2 = aggregate(select(data2, -region), by = list(data2$region), sum) %>%
#   select(., -Group.1) %>% t %>% as.data.frame
# colnames(data2) = c("20k-2M", "2M-12M", ">12M")
# data2$Contact_diff = data2$`20k-2M` - data2$`2M-12M`
# write.table(data2, file = "k562.data2.txt", quote = F, sep = "\t")
data2 = read.table("k562.data2.txt", header = T, sep = "\t", check.names = F)

data = merge(data1, data2, by.x = "Row.names", by.y = "row.names", all.x = T) %>% .[order(.$order), ]

# #拟合Gene_diff，并将其归一化
# tmp = mgcv::gam(Gene_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Gene_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Gene_diff = tmp$Gene_diff
# #拟合Contact_diff，并将其归一化
# tmp = mgcv::gam(Contact_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Contact_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Contact_diff = tmp$Contact_diff

data = data %>% 
  melt(id = c("Row.names", "order")) %>% 
  mutate(order = as.numeric(.$order)) %>% 
  mutate(variable = as.character(.$variable))
phase = data[data$variable == "type", ]$value
data[data$variable == "type", ]$variable = data[data$variable == "type", ]$value

#添加分组信息，周期信息为一组，Gene_diff和Contact_diff为一组
data$group = data$variable
data[data$variable %in% unique(phase), ]$group = "Cell phase"
data[data$variable %in% c("Contact_diff", "Gene_diff"), ]$group = "Diff"
data[data$group == "Cell phase", ]$value = 1
data$group = factor(data$group, levels = c("20k-2M", "2M-12M", ">12M", "S.Score", "G2M.Score", "Diff", "Pesudotime_dis", "Cell phase"))
data$value = as.numeric(data$value)

mytheme = theme(panel.border = element_blank(), 
                panel.grid = element_blank(), 
                axis.line = element_line(colour = "black"), 
                axis.text.x = element_blank(), 
                axis.text.y = element_text(colour = "black"), 
                axis.ticks.x = element_blank())

p1 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == ">12M", ], mapping = aes(x = order, y = value), color = "#E64B35") +
  labs(x = "", y = "", title = "The propotion of contact (>12M)") +
  theme_bw() +
  mytheme

p2 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "2M-12M", ], mapping = aes(x = order, y = value), color = "#4DBBD5") +
  labs(x = "", y = "", title = "The propotion of contact (2M-12M)") +
  theme_bw() +
  mytheme

p3 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "20k-2M", ], mapping = aes(x = order, y = value), color = "#00A087") +
  labs(x = "", y = "", title = "The propotion of contact (20k-2M)") +
  theme_bw() +
  mytheme

p4 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "S.Score", ], mapping = aes(x = order, y = value), color = "#3C5488") +
  labs(x = "", y = "", title = "S phase mark gene score") +
  theme_bw() +
  mytheme

p5 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "G2M.Score", ], mapping = aes(x = order, y = value), color = "#F39B7F") +
  labs(x = "", y = "", title = "G2M phase marker gene score") +
  theme_bw() +
  mytheme

# p6 = ggplot(data[data$group == "Diff", ], aes(x = order, y = value, color = variable)) +
#   geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
#   geom_line(linewidth = 1) +
#   labs(x = "", y = "", title = "Diff", color = "Diff") +
#   theme_bw() +
#   mytheme

p6 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "Diff", ], mapping = aes(x = order, y = value, color = variable), linewidth = 1) +
  scale_y_continuous(name = "Phase marker gene score Diff", sec.axis = sec_axis(~ (. + 0.3), name = "Contact (2k-2M & 2M-12M) Diff")) +
  labs(x = "", y = "", title = "Marker gene and contact Diff", color = "Diff") +
  theme_bw() +
  mytheme +
  theme(legend.position = "none")

p7 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 182, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 182, xmax = 541, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 541, xmax = 5709, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 5709, xmax = 34980, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34980, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "Pesudotime_dis", ], mapping = aes(x = order, y = value), color = "#8491B4") +
  labs(x = "", y = "", title = "Pseudotime") +
  theme_bw() +
  mytheme

# p8 = ggplot(data[data$group == "Cell phase", ] %>% mutate(variable = factor(.$variable, levels = c("Post-M", "G1", "Early_S", "Late_S-G2", "Pre-M"))), aes(x = order, y = value, color = variable)) +
#   geom_line(linewidth = 2) +
#   labs(x = "", y = "", title = "", color = "Cell Phase") +
#   scale_color_manual(values = c(`Post-M` = "red", `G1` = "blue", `Early_S` = "green", `Late_S-G2`="purple",`Pre-M`="orange")) +
#   theme_bw() +
#   theme(panel.border = element_blank(), 
#         panel.grid = element_blank(), 
#         axis.text = element_blank(), 
#         axis.ticks = element_blank())

layout = "
  AABBCC
  DDEEFF
  GG####
  "

p = p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
  plot_layout(design = layout, guides = 'collect') + 
  plot_annotation(tag_levels = "A")

pdf(file = "pseudotime.k562.pdf", width = 14, height = 10, pointsize = 8)
print(p)
dev.off()


##### patski order by contact #####
finaldat = read.table("patski.finaldat.txt", header = T, sep = "\t")
finaldat$cellid = rownames(finaldat)

data = read.table("patski.dis.stat.txt", header = T, sep = "\t", row.names = 1)
data = data[, c("S.Score", "G2M.Score")] %>%
  mutate(Gene_diff = .$S.Score - .$G2M.Score) %>%
  mutate(Pesudotime_dis = 1:nrow(.) / nrow(.)) %>%
  .[finaldat$cellid, ] %>%
  .[grep("^NA", rownames(.), invert = T), ] %>% 
  mutate(order = 1:nrow(.))

data1 = merge(data, finaldat[, c("cellid", "type")], by.x = "row.names", by.y = "cellid", all.x = T)

# data2 = finaldat[finaldat$cellid %in% rownames(data), ]
# cellid = data2$cellid
# data2 = data2[, grep("^bin", colnames(data2))]  %>%
#   sapply(., percentage) %>%
#   as.data.frame %>% t %>% as.data.frame
# colnames(data2) = cellid
# data2[is.na(data2)] = 0
# data2 = scale(data2, center = F, scale = colSums(data2)) %>%
#   as.data.frame %>%
#   mutate(region = factor(c(rep("20k-2M", 52), rep("2M-12M", 20), rep(">12M", nrow(.)-72)),
#                          levels = c("20k-2M", "2M-12M", ">12M")))
# data2 = aggregate(select(data2, -region), by = list(data2$region), sum) %>%
#   select(., -Group.1) %>% t %>% as.data.frame
# colnames(data2) = c("20k-2M", "2M-12M", ">12M")
# data2$Contact_diff = data2$`20k-2M` - data2$`2M-12M`
# write.table(data2, file = "patski.data2.txt", quote = F, sep = "\t")
data2 = read.table("patski.data2.txt", header = T, sep = "\t", check.names = F)

data = merge(data1, data2, by.x = "Row.names", by.y = "row.names", all.x = T) %>% .[order(.$order), ]

# #拟合Gene_diff，并将其归一化
# tmp = mgcv::gam(Gene_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Gene_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Gene_diff = tmp$Gene_diff
# #拟合Contact_diff，并将其归一化
# tmp = mgcv::gam(Contact_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Contact_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Contact_diff = tmp$Contact_diff

data = data %>% 
  melt(id = c("Row.names", "order")) %>% 
  mutate(order = as.numeric(.$order)) %>% 
  mutate(variable = as.character(.$variable))
phase = data[data$variable == "type", ]$value
data[data$variable == "type", ]$variable = data[data$variable == "type", ]$value

#添加分组信息，周期信息为一组，Gene_diff和Contact_diff为一组
data$group = data$variable
data[data$variable %in% unique(phase), ]$group = "Cell phase"
data[data$variable %in% c("Contact_diff", "Gene_diff"), ]$group = "Diff"
data[data$group == "Cell phase", ]$value = 1
data$group = factor(data$group, levels = c("20k-2M", "2M-12M", ">12M", "S.Score", "G2M.Score", "Diff", "Pesudotime_dis", "Cell phase"))
data$value = as.numeric(data$value)

mytheme = theme(panel.border = element_blank(), 
                panel.grid = element_blank(), 
                axis.line = element_line(colour = "black"), 
                axis.text.x = element_blank(), 
                axis.text.y = element_text(colour = "black"), 
                axis.ticks.x = element_blank())

p1 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == ">12M", ], mapping = aes(x = order, y = value), color = "#E64B35") +
  labs(x = "", y = "", title = "The propotion of contact (>12M)") +
  theme_bw() +
  mytheme

p2 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "2M-12M", ], mapping = aes(x = order, y = value), color = "#4DBBD5") +
  labs(x = "", y = "", title = "The propotion of contact (2M-12M)") +
  theme_bw() +
  mytheme

p3 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "20k-2M", ], mapping = aes(x = order, y = value), color = "#00A087") +
  labs(x = "", y = "", title = "The propotion of contact (20k-2M)") +
  theme_bw() +
  mytheme

p4 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "S.Score", ], mapping = aes(x = order, y = value), color = "#3C5488") +
  labs(x = "", y = "", title = "S phase mark gene score") +
  theme_bw() +
  mytheme

p5 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "G2M.Score", ], mapping = aes(x = order, y = value), color = "#F39B7F") +
  labs(x = "", y = "", title = "G2M phase marker gene score") +
  theme_bw() +
  mytheme

# p6 = ggplot(data[data$group == "Diff", ], aes(x = order, y = value, color = variable)) +
#   geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
#   geom_line(linewidth = 1) +
#   labs(x = "", y = "", title = "Diff", color = "Diff") +
#   theme_bw() +
#   mytheme

p6 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "Diff", ], mapping = aes(x = order, y = value, color = variable), linewidth = 1) +
  scale_y_continuous(name = "Phase marker gene score Diff", sec.axis = sec_axis(~ (. + 0.3), name = "Contact (2k-2M & 2M-12M) Diff")) +
  labs(x = "", y = "", title = "Marker gene and contact Diff", color = "Diff") +
  theme_bw() +
  mytheme +
  theme(legend.position = "none")

p7 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 70, ymin=-Inf, ymax = Inf), fill = "red", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 70, xmax = 1396, ymin=-Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 1396, xmax = 28842, ymin=-Inf, ymax = Inf), fill = "green", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 28842, xmax = 34875, ymin=-Inf, ymax = Inf), fill = "purple", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 34875, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "orange", alpha = 0.1) +
  geom_smooth(data[data$group == "Pesudotime_dis", ], mapping = aes(x = order, y = value), color = "#8491B4") +
  labs(x = "", y = "", title = "Pseudotime") +
  theme_bw() +
  mytheme

# p8 = ggplot(data[data$group == "Cell phase", ] %>% mutate(variable = factor(.$variable, levels = c("Post-M", "G1", "Early_S", "Late_S-G2", "Pre-M"))), aes(x = order, y = value, color = variable)) +
#   geom_line(linewidth = 2) +
#   labs(x = "", y = "", title = "", color = "Cell Phase") +
#   scale_color_manual(values = c(`Post-M` = "red", `G1` = "blue", `Early_S` = "green", `Late_S-G2`="purple",`Pre-M`="orange")) +
#   theme_bw() +
#   theme(panel.border = element_blank(), 
#         panel.grid = element_blank(), 
#         axis.text = element_blank(), 
#         axis.ticks = element_blank())

layout = "
  AABBCC
  DDEEFF
  GG####
  "

p = p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
  plot_layout(design = layout, guides = 'collect') + 
  plot_annotation(tag_levels = "A")

pdf(file = "pseudotime.patski.pdf", width = 14, height = 10, pointsize = 8)
print(p)
dev.off()

##### k562 order by gene #####
finaldat = read.table("k562.finaldat.txt", header = T, sep = "\t")
finaldat$cellid = rownames(finaldat)

data = read.table("k562.dis.stat.txt", header = T, sep = "\t", row.names = 1)
data = data[, c("S.Score", "G2M.Score", "Phase")] %>%
  mutate(Gene_diff = .$S.Score - .$G2M.Score) %>%
  mutate(Pesudotime_dis = 1:nrow(.)) %>%
  .[finaldat$cellid, ] %>%
  .[grep("^NA", rownames(.), invert = T), ] %>% 
  mutate(order = 1:nrow(.))

data2 = read.table("k562.data2.txt", header = T, sep = "\t", check.names = F)

data = merge(data, data2, by = "row.names", all.x = T) %>% .[order(.$Pesudotime_dis), ]

# #拟合Gene_diff，并将其归一化
# tmp = mgcv::gam(Gene_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Gene_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Gene_diff = tmp$Gene_diff
# #拟合Contact_diff，并将其归一化
# tmp = mgcv::gam(Contact_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Contact_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Contact_diff = tmp$Contact_diff

data = data %>% 
  melt(id = c("Row.names", "Pesudotime_dis")) %>% 
  mutate(Pesudotime_dis = as.numeric(.$Pesudotime_dis)) %>% 
  mutate(variable = as.character(.$variable))
phase = data[data$variable == "Phase", ]$value
data[data$variable == "Phase", ]$variable = data[data$variable == "Phase", ]$value

#添加分组信息，周期信息为一组，Gene_diff和Contact_diff为一组
data$group = data$variable
data[data$variable %in% unique(phase), ]$group = "Cell phase"
data[data$variable %in% c("Contact_diff", "Gene_diff"), ]$group = "Diff"
data[data$group == "Cell phase", ]$value = 1
data$group = factor(data$group, levels = c("20k-2M", "2M-12M", ">12M", "S.Score", "G2M.Score", "Diff", "order", "Cell phase"))
data$value = as.numeric(data$value)

mytheme = theme(panel.border = element_blank(), 
                panel.grid = element_blank(), 
                axis.line = element_line(colour = "black"), 
                axis.text.x = element_blank(), 
                axis.text.y = element_text(colour = "black"), 
                axis.ticks.x = element_blank())

p1 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == ">12M", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#E64B35") +
  labs(x = "", y = "", title = "The propotion of contact (>12M)") +
  theme_bw() +
  mytheme

p2 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "2M-12M", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#4DBBD5") +
  labs(x = "", y = "", title = "The propotion of contact (2M-12M)") +
  theme_bw() +
  mytheme

p3 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "20k-2M", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#00A087") +
  labs(x = "", y = "", title = "The propotion of contact (20k-2M)") +
  theme_bw() +
  mytheme

p4 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "S.Score", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#3C5488") +
  labs(x = "", y = "", title = "S phase mark gene score") +
  theme_bw() +
  mytheme

p5 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "G2M.Score", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#F39B7F") +
  labs(x = "", y = "", title = "G2M phase marker gene score") +
  theme_bw() +
  mytheme

# p6 = ggplot(data[data$group == "Diff", ], aes(x = order, y = value, color = variable)) +
#   geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
#   geom_line(linewidth = 1) +
#   labs(x = "", y = "", title = "Diff", color = "Diff") +
#   theme_bw() +
#   mytheme

p6 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "Diff", ], mapping = aes(x = Pesudotime_dis, y = value, color = variable), linewidth = 1) +
  scale_y_continuous(name = "Phase marker gene score Diff", sec.axis = sec_axis(~ (. + 0.3), name = "Contact (2k-2M & 2M-12M) Diff")) +
  labs(x = "", y = "", title = "Marker gene and contact Diff", color = "Diff") +
  theme_bw() +
  mytheme +
  theme(legend.position = "none")

p7 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "order", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#8491B4") +
  labs(x = "", y = "", title = "Contact order") +
  theme_bw() +
  mytheme

# p8 = ggplot(data[data$group == "Cell phase", ] %>% mutate(variable = factor(.$variable, levels = c("G1", "S", "G2M"))), aes(x = Pesudotime_dis, y = value, color = variable)) +
#   geom_line(linewidth = 2) +
#   labs(x = "", y = "", title = "", color = "Cell Phase") +
#   scale_color_manual(values = c(`G1` = "#F8766D", `S` = "#00BA38", `G2M`="#619CFF")) +
#   theme_bw() +
#   theme(panel.border = element_blank(), 
#         panel.grid = element_blank(), 
#         axis.text = element_blank(), 
#         axis.ticks = element_blank())


layout = "
  AABBCC
  DDEEFF
  GG####
  "

p = p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
  plot_layout(design = layout, guides = 'collect') + 
  plot_annotation(tag_levels = "A")

pdf(file = "pseudotime.gene.k562.pdf", width = 14, height = 10, pointsize = 8)
print(p)
dev.off()

##### patski order by gene #####
finaldat = read.table("patski.finaldat.txt", header = T, sep = "\t")
finaldat$cellid = rownames(finaldat)

data = read.table("patski.dis.stat.txt", header = T, sep = "\t", row.names = 1)
data = data[, c("S.Score", "G2M.Score", "Phase")] %>%
  mutate(Gene_diff = .$S.Score - .$G2M.Score) %>%
  mutate(Pesudotime_dis = 1:nrow(.)) %>%
  .[finaldat$cellid, ] %>%
  .[grep("^NA", rownames(.), invert = T), ] %>% 
  mutate(order = 1:nrow(.))

data2 = read.table("patski.data2.txt", header = T, sep = "\t", check.names = F)

data = merge(data, data2, by = "row.names", all.x = T) %>% .[order(.$Pesudotime_dis), ]

# #拟合Gene_diff，并将其归一化
# tmp = mgcv::gam(Gene_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Gene_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Gene_diff = tmp$Gene_diff
# #拟合Contact_diff，并将其归一化
# tmp = mgcv::gam(Contact_diff ~ s(order, bs = "cs"), data = data)
# tmp = data.frame(order = tmp$model$order, Contact_diff = (tmp$fitted.values - min(tmp$fitted.values)) / (max(tmp$fitted.values) - min(tmp$fitted.values))) %>%
#   .[order(.$order), ]
# data$Contact_diff = tmp$Contact_diff

data = data %>% 
  melt(id = c("Row.names", "Pesudotime_dis")) %>% 
  mutate(Pesudotime_dis = as.numeric(.$Pesudotime_dis)) %>% 
  mutate(variable = as.character(.$variable))
phase = data[data$variable == "Phase", ]$value
data[data$variable == "Phase", ]$variable = data[data$variable == "Phase", ]$value

#添加分组信息，周期信息为一组，Gene_diff和Contact_diff为一组
data$group = data$variable
data[data$variable %in% unique(phase), ]$group = "Cell phase"
data[data$variable %in% c("Contact_diff", "Gene_diff"), ]$group = "Diff"
data[data$group == "Cell phase", ]$value = 1
data$group = factor(data$group, levels = c("20k-2M", "2M-12M", ">12M", "S.Score", "G2M.Score", "Diff", "order", "Cell phase"))
data$value = as.numeric(data$value)

mytheme = theme(panel.border = element_blank(), 
                panel.grid = element_blank(), 
                axis.line = element_line(colour = "black"), 
                axis.text.x = element_blank(), 
                axis.text.y = element_text(colour = "black"), 
                axis.ticks.x = element_blank())

p1 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == ">12M", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#E64B35") +
  labs(x = "", y = "", title = "The propotion of contact (>12M)") +
  theme_bw() +
  mytheme

p2 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "2M-12M", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#4DBBD5") +
  labs(x = "", y = "", title = "The propotion of contact (2M-12M)") +
  theme_bw() +
  mytheme

p3 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "20k-2M", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#00A087") +
  labs(x = "", y = "", title = "The propotion of contact (20k-2M)") +
  theme_bw() +
  mytheme

p4 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "S.Score", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#3C5488") +
  labs(x = "", y = "", title = "S phase mark gene score") +
  theme_bw() +
  mytheme

p5 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "G2M.Score", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#F39B7F") +
  labs(x = "", y = "", title = "G2M phase marker gene score") +
  theme_bw() +
  mytheme

# p6 = ggplot(data[data$group == "Diff", ], aes(x = order, y = value, color = variable)) +
#   geom_rect(mapping = aes(xmin = 0, xmax = 11227, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 11227, xmax = 27536, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
#   geom_rect(mapping = aes(xmin = 27536, xmax = 40595, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
#   geom_line(linewidth = 1) +
#   labs(x = "", y = "", title = "Diff", color = "Diff") +
#   theme_bw() +
#   mytheme

p6 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "Diff", ], mapping = aes(x = Pesudotime_dis, y = value, color = variable), linewidth = 1) +
  scale_y_continuous(name = "Phase marker gene score Diff", sec.axis = sec_axis(~ (. + 0.3), name = "Contact (2k-2M & 2M-12M) Diff")) +
  labs(x = "", y = "", title = "Marker gene and contact Diff", color = "Diff") +
  theme_bw() +
  mytheme +
  theme(legend.position = "none")

p7 = ggplot() +
  geom_rect(mapping = aes(xmin = 0, xmax = 18595, ymin=-Inf, ymax = Inf), fill = "#F8766D", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 18595, xmax = 27529, ymin=-Inf, ymax = Inf), fill = "#00BA38", alpha = 0.1) +
  geom_rect(mapping = aes(xmin = 27529, xmax = 35515, ymin=-Inf, ymax = Inf), fill = "#619CFF", alpha = 0.1) +
  geom_smooth(data[data$group == "order", ], mapping = aes(x = Pesudotime_dis, y = value), color = "#8491B4") +
  labs(x = "", y = "", title = "Contact order") +
  theme_bw() +
  mytheme

# p8 = ggplot(data[data$group == "Cell phase", ] %>% mutate(variable = factor(.$variable, levels = c("G1", "S", "G2M"))), aes(x = Pesudotime_dis, y = value, color = variable)) +
#   geom_line(linewidth = 2) +
#   labs(x = "", y = "", title = "", color = "Cell Phase") +
#   scale_color_manual(values = c(`G1` = "#F8766D", `S` = "#00BA38", `G2M`="#619CFF")) +
#   theme_bw() +
#   theme(panel.border = element_blank(), 
#         panel.grid = element_blank(), 
#         axis.text = element_blank(), 
#         axis.ticks = element_blank())


layout = "
  AABBCC
  DDEEFF
  GG####
  "

p = p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
  plot_layout(design = layout, guides = 'collect') + 
  plot_annotation(tag_levels = "A")

pdf(file = "pseudotime.gene.patski.pdf", width = 14, height = 10, pointsize = 8)
print(p)
dev.off()