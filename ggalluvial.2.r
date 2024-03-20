library(sankeyD3)
library(patchwork)
library(webshot)

nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1
sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
  Value = "weight", NodeID = "name",nodeWidth =10,units = 'TWh',
  height=300,width=300,colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
  numberFormat=".0f",fontSize = 8)


library(magrittr)
URL <- paste0("https://cdn.rawgit.com/christophergandrud/networkD3/",
              "master/JSONdata/energy.json")
Energy <- jsonlite::fromJSON("C:\\Users\\dell\\Desktop\\energy.json")
p1 = sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
  Target = "target", Value = "value", NodeID = "name",
  units = "TWh", fontSize = 12, nodeWidth = 30)
p2 = sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
  Target = "target", Value = "value", NodeID = "name",
  units = "TWh", fontSize = 12, nodeWidth = 30) %>% saveNetwork(file = 'Net1.html')


colorlist = c(brewer.pal(8,"Dark2"),brewer.pal(12,"Paired"))
colorlist2 = c(colorlist[1],colorlist[3],colorlist[4],colorlist[4],colorlist[5],colorlist[6],colorlist[6],colorlist[7],colorlist[7],colorlist[8],colorlist[11],colorlist[12],colorlist[12],colorlist[14],colorlist[16],colorlist[18],colorlist[19],colorlist[19],colorlist[20])
countCellCount = readRDS("countCellCount.RData")
nodes = c(c("Meninges","Cortex_L1","Cortex-L2/3","Cortex_L4","Cortex_L5","Cortex_L6","Lateral-ventral cortex","Cortical amygdalar area","Posterior amygdalar nucleus","Fiber tract","Subiculum","Stratum oriens of CA1", "CA1","Stratum lacunosum/raditum of CA1","Molecular layer of dentate gyrus","Dentate gyrus","Cavity","Thalamus","Midbrain","Substandia nigra/Ventral tegmental area"),c("VLMC2","TEGLU7","TEGLU4","TEGLU8","TEGLU10","TEGLU2","TEGLU3","TEGLU12","TEGLU15","TEGLU22","TEGLU6","TEGLU23","TEGLU24","ACTE2","DGGRC2","DEGLU1","MOL2","ACNT1","TEINH18"))
# nodesdat = data.frame(name=nodes,yOrder=1:39)
nodesdat = data.frame(name=nodes,yOrder=1:39,nodecolor = c(colorlist,colorlist2))
nodesdat$name = factor(nodesdat$name, levels = nodes)
nodesdat1 = data.frame(celltype=nodes, target  = 0:(length(nodes)-1))
nodesdat2 = data.frame(annotation=nodes, source = 0:(length(nodes)-1))
countCellCountdat = left_join(left_join(countCellCount,nodesdat1, by="celltype"),nodesdat2,by="annotation")
countCellCountdat$linkcolor = "grey"
countCellCountdat[countCellCountdat$annotation=="Meninges" & countCellCountdat$celltype=="VLMC2",]$linkcolor = colorlist[1]
# countCellCountdat[countCellCountdat$annotation=="Cortex_L1" & countCellCountdat$celltype=="ACTE2",]$linkcolor = colorlist[2]
countCellCountdat[countCellCountdat$annotation=="Cortex-L2/3" & countCellCountdat$celltype=="TEGLU7",]$linkcolor = colorlist[3]
countCellCountdat[countCellCountdat$annotation=="Cortex_L4" & countCellCountdat$celltype=="TEGLU4",]$linkcolor = colorlist[4]
countCellCountdat[countCellCountdat$annotation=="Cortex_L4" & countCellCountdat$celltype=="TEGLU8",]$linkcolor = colorlist[4]
countCellCountdat[countCellCountdat$annotation=="Cortex_L5" & countCellCountdat$celltype=="TEGLU10",]$linkcolor = colorlist[5]
countCellCountdat[countCellCountdat$annotation=="Cortex_L6" & countCellCountdat$celltype=="TEGLU2",]$linkcolor = colorlist[6]
countCellCountdat[countCellCountdat$annotation=="Cortex_L6" & countCellCountdat$celltype=="TEGLU3",]$linkcolor = colorlist[6]
countCellCountdat[countCellCountdat$annotation=="Lateral-ventral cortex" & countCellCountdat$celltype=="TEGLU12",]$linkcolor = colorlist[7]
countCellCountdat[countCellCountdat$annotation=="Lateral-ventral cortex" & countCellCountdat$celltype=="TEGLU15",]$linkcolor = colorlist[7]
countCellCountdat[countCellCountdat$annotation=="Cortical amygdalar area" & countCellCountdat$celltype=="TEGLU22",]$linkcolor = colorlist[8]
countCellCountdat[countCellCountdat$annotation=="Subiculum" & countCellCountdat$celltype=="TEGLU6",]$linkcolor = colorlist[11]
countCellCountdat[countCellCountdat$annotation=="Stratum oriens of CA1" & countCellCountdat$celltype=="TEGLU23",]$linkcolor = colorlist[12]
countCellCountdat[countCellCountdat$annotation=="Stratum oriens of CA1" & countCellCountdat$celltype=="TEGLU24",]$linkcolor = colorlist[12]
countCellCountdat[countCellCountdat$annotation=="Stratum lacunosum/raditum of CA1" & countCellCountdat$celltype=="ACTE2",]$linkcolor = colorlist[14]
countCellCountdat[countCellCountdat$annotation=="Dentate gyrus" & countCellCountdat$celltype=="DGGRC2",]$linkcolor = colorlist[16]
countCellCountdat[countCellCountdat$annotation=="Thalamus" & countCellCountdat$celltype=="DEGLU1",]$linkcolor = colorlist[18]
countCellCountdat[countCellCountdat$annotation=="Midbrain" & countCellCountdat$celltype=="MOL2",]$linkcolor = colorlist[19]
countCellCountdat[countCellCountdat$annotation=="Midbrain" & countCellCountdat$celltype=="ACNT1",]$linkcolor = colorlist[19]

countCellCountdat[countCellCountdat$annotation=="Meninges" & countCellCountdat$celltype=="VLMC2",]$linkcolor = "color1"
countCellCountdat[countCellCountdat$annotation=="Cortex-L2/3" & countCellCountdat$celltype=="TEGLU7",]$linkcolor = "color3"
countCellCountdat[countCellCountdat$annotation=="Cortex_L4" & countCellCountdat$celltype=="TEGLU4",]$linkcolor = "color4"
countCellCountdat[countCellCountdat$annotation=="Cortex_L4" & countCellCountdat$celltype=="TEGLU8",]$linkcolor = "color4"
countCellCountdat[countCellCountdat$annotation=="Cortex_L5" & countCellCountdat$celltype=="TEGLU10",]$linkcolor = "color5"
countCellCountdat[countCellCountdat$annotation=="Cortex_L6" & countCellCountdat$celltype=="TEGLU2",]$linkcolor = "color6"
countCellCountdat[countCellCountdat$annotation=="Cortex_L6" & countCellCountdat$celltype=="TEGLU3",]$linkcolor = "color6"
countCellCountdat[countCellCountdat$annotation=="Lateral-ventral cortex" & countCellCountdat$celltype=="TEGLU12",]$linkcolor = "color7"
countCellCountdat[countCellCountdat$annotation=="Lateral-ventral cortex" & countCellCountdat$celltype=="TEGLU15",]$linkcolor = "color7"
countCellCountdat[countCellCountdat$annotation=="Cortical amygdalar area" & countCellCountdat$celltype=="TEGLU22",]$linkcolor = "color8"
countCellCountdat[countCellCountdat$annotation=="Subiculum" & countCellCountdat$celltype=="TEGLU6",]$linkcolor = "color11"
countCellCountdat[countCellCountdat$annotation=="Stratum oriens of CA1" & countCellCountdat$celltype=="TEGLU23",]$linkcolor = "color12"
countCellCountdat[countCellCountdat$annotation=="Stratum oriens of CA1" & countCellCountdat$celltype=="TEGLU24",]$linkcolor = "color12"
countCellCountdat[countCellCountdat$annotation=="Stratum lacunosum/raditum of CA1" & countCellCountdat$celltype=="ACTE2",]$linkcolor = "color14"
countCellCountdat[countCellCountdat$annotation=="Dentate gyrus" & countCellCountdat$celltype=="DGGRC2",]$linkcolor = "color16"
countCellCountdat[countCellCountdat$annotation=="Thalamus" & countCellCountdat$celltype=="DEGLU1",]$linkcolor = "color18"
countCellCountdat[countCellCountdat$annotation=="Midbrain" & countCellCountdat$celltype=="MOL2",]$linkcolor = "color19"
countCellCountdat[countCellCountdat$annotation=="Midbrain" & countCellCountdat$celltype=="ACNT1",]$linkcolor = "color19"

# colrr = c(colorlist[1],colorlist[3],colorlist[4],colorlist[5],colorlist[6],colorlist[7],colorlist[8],colorlist[11],colorlist[12],colorlist[14],colorlist[16],colorlist[18],colorlist[19])
# my_color <- 'd3.scaleOrdinal() .domain(["grey","#1B9E77","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#B2DF8A","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#FFFF99"]) .range(["grey","#1B9E77","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#B2DF8A","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#FFFF99"])'
my_color <- 'd3.scaleOrdinal() .domain(["grey","color1","color3","color4","color5","color6","color7","color8","color11","color12","color14","color16","color18","color19"]) .range(["#80808039","#1B9E77","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#B2DF8A","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#FFFF99"])'
pdf("mb.raw.spatial.fig.ggalluvial.2.pdf")
sankeyNetwork(Links = countCellCountdat[,c("source","target","Freq","linkcolor")], Nodes = nodesdat, Source = "source",
              Target = "target", Value = "Freq", NodeID = "name", yOrderComparator = htmlwidgets::JS("function(a, b) { return a.yOrder - b.yOrder; }"),
              LinkGroup = "linkcolor", colourScale = my_color, NodeColor = "nodecolor", units = "TWh", fontSize = 12, nodeWidth = 30,showNodeValues=F)
dev.off()

sankeyNetwork(Links = countCellCountdat[,c("source","target","Freq","linkcolor")], Nodes = nodesdat, Source = "source",
              Target = "target", Value = "Freq", NodeID = "name", yOrderComparator = htmlwidgets::JS("function(a, b) { return a.yOrder - b.yOrder; }"),
              LinkGroup = "linkcolor", colourScale = my_color, NodeColor = "nodecolor", units = "TWh", fontSize = 12, nodeWidth = 30,showNodeValues=F) %>% saveNetwork(file = 'Net2.html')


if(!is_phantomjs_installed()){
  install_phantomjs()
}
p = sankeyNetwork(Links = countCellCountdat[,c("source","target","Freq","linkcolor")], Nodes = nodesdat, Source = "source",
              Target = "target", Value = "Freq", NodeID = "name", yOrderComparator = htmlwidgets::JS("function(a, b) { return a.yOrder - b.yOrder; }"),
              LinkGroup = "linkcolor", colourScale = my_color, NodeColor = "nodecolor", units = "TWh", fontSize = 12, nodeWidth = 30,showNodeValues=F)
saveNetwork(p,"sankey.html")
webshot("sankey.html" , "sankey.pdf")
