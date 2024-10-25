data <- c(1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100)
x = data
# 计算四分位数
Q1 <- quantile(data, 0.25)
Q3 <- quantile(data, 0.75)
IQR <- Q3 - Q1

# 定义离群点的界限
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# 判定离群点
outliers <- data[data < lower_bound | data > upper_bound]

# 打印结果
print("离群点:")
print(outliers)


x = sort(x, decreasing = T)
cutoff = x[1]
for (i in seq(1, length(x)-1)) {
  if (x[i] > x[i]*10) {
    cutoff = x[i]
  }
}
