library(canvasXpress)

y <- read.table("http://www.canvasxpress.org/data/cX-irist-dat.txt", header=TRUE, sep="\t", 
                quote="", row.names=1, fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
# 读取分组信息
z <- read.table("http://www.canvasxpress.org/data/cX-irist-var.txt", header=TRUE, sep= "\t", 
                quote="", row.names=1, fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
# 绘制三维散点图，主要参数为数据、分组、分组列、置信椭圆列、图表类型以及相关标签               
canvasXpress(data      = y,
             varAnnot  = z,
             colorBy   = "Species",
             ellipseBy = "Species",
             graphType = "Scatter3D",
             title     = "Iris Data Set",
             xAxis     = list("Sepal.Length"),
             yAxis     = list("Petal.Width"),
             zAxis     = list("Petal.Length"))

cXsaveAsPDF(cx, "example.pdf")
webshot("test.2.html" , "test.2.html.pdf")
