library(gridExtra)
library(grid)
d - head(iris, 3)
g - tableGrob(d)
grid.newpage()
grid.draw(g)