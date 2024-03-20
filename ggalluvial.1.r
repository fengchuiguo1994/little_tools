library(ggalluvial)

data(mtcars)

ggplot(mtcars,
  aes(axis1 = factor(gear, levels = c(5,4,3)),  # 调整 gear 的顺序为 3, 4, 5
      axis2 = factor(cyl, levels = c(8,6,4)),   # 调整 cyl 的顺序为 4, 6, 8
      y = mpg)) +
  geom_alluvium(aes(fill = factor(vs)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()