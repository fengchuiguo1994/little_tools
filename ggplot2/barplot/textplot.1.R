library(ggplot2)

# 创建一个空的 ggplot 图形
p <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Hello, World!", size = 5, color = "blue") +
  xlim(-1, 1) + ylim(-1, 1) +
  theme_void()
print(p)
# 保存为图片
# ggsave("text_image.png", plot = p, width = 5, height = 5)