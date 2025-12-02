library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 定义扇形分割函数
split_into_sectors <- function(data, points_per_sector = 1000) {
  # 计算数据的中心点
  center_x <- mean(data$X)
  center_y <- mean(data$Y)
  
  # 计算每个点相对于中心的极坐标
  delta_x <- data$X - center_x
  delta_y <- data$Y - center_y
  
  # 计算半径和角度（弧度）
  r <- sqrt(delta_x^2 + delta_y^2)
  theta_rad <- atan2(delta_y, delta_x)
  
  # 将角度转换为度数（0°到360°范围）
  theta_deg <- (theta_rad * 180 / pi) %% 360
  
  # 创建包含极坐标的数据框
  data_polar <- data %>%
    mutate(
      r = r,
      theta_rad = theta_rad,
      theta_deg = theta_deg
    )
  
  # 按角度排序
  data_sorted <- data_polar %>% arrange(theta_deg)
  
  # 划分扇形（每1000个连续点一个扇形）
  n <- nrow(data_sorted)
  sectors <- list()
  
  for (i in seq(1, n, by = points_per_sector)) {
    end_idx <- min(i + points_per_sector - 1, n)
    sector_data <- data_sorted[i:end_idx, ]
    if (nrow(sector_data) == 0) next
    
    # 计算该扇形的角度范围
    min_angle <- min(sector_data$theta_deg)
    max_angle <- max(sector_data$theta_deg)
    
    # 存储扇形信息
    sector_info <- list(
      sector_id = length(sectors) + 1,
      points = sector_data %>% select(X, Y),
      min_angle = min_angle,
      max_angle = max_angle,
      center_angle = (min_angle + max_angle) / 2,
      point_count = nrow(sector_data)
    )
    
    sectors[[length(sectors) + 1]] <- sector_info
  }
  
  return(list(sectors = sectors, center_x = center_x, center_y = center_y))
}

# 生成示例数据
set.seed(42)
n_points <- 10000
data <- data.frame(
  X = rnorm(n_points, 0, 10),
  Y = rnorm(n_points, 0, 10)
)

# 执行扇形分割
result <- split_into_sectors(data, 1000)
sectors <- result$sectors
center_x <- result$center_x
center_y <- result$center_y

# 打印结果
cat("总数据点数:", n_points, "\n")
cat("划分的扇形数量:", length(sectors), "\n\n")

for (i in seq_along(sectors)) {
  sector <- sectors[[i]]
  cat(sprintf("扇形 %d: %d 个点, 角度范围: %.1f° - %.1f°\n", 
              sector$sector_id, sector$point_count, 
              sector$min_angle, sector$max_angle))
}


# 可视化函数
plot_sectors <- function(sectors, center_x, center_y) {
  # 准备绘图数据
  plot_data <- data.frame()
  sector_info <- data.frame()
  
  for (sector in sectors) {
    sector_points <- sector$points %>%
      mutate(sector_id = factor(sector$sector_id))
    plot_data <- rbind(plot_data, sector_points)
    
    sector_info <- rbind(sector_info, data.frame(
      sector_id = factor(sector$sector_id),
      min_angle = sector$min_angle,
      max_angle = sector$max_angle
    ))
  }
  
  # 计算最大半径用于绘制边界线
  max_r <- max(sqrt((plot_data$X - center_x)^2 + (plot_data$Y - center_y)^2)) * 1.1
  
  # 创建基础图形
  p <- ggplot() +
    # 绘制数据点
    geom_point(data = plot_data, 
               aes(x = X, y = Y, color = sector_id), 
               size = 1, alpha = 0.6) +
    
    # 绘制中心点
    geom_point(data = data.frame(x = center_x, y = center_y), 
               aes(x = x, y = y), 
               color = "red", size = 4, shape = 4, stroke = 2) +
    
    # 绘制扇形边界线
    geom_segment(data = sector_info,
                 aes(x = center_x, y = center_y,
                     xend = center_x + max_r * cos(min_angle * pi / 180),
                     yend = center_y + max_r * sin(min_angle * pi / 180),
                     color = sector_id),
                 linetype = "dashed", alpha = 0.7) +
    
    geom_segment(data = sector_info,
                 aes(x = center_x, y = center_y,
                     xend = center_x + max_r * cos(max_angle * pi / 180),
                     yend = center_y + max_r * sin(max_angle * pi / 180),
                     color = sector_id),
                 linetype = "dashed", alpha = 0.7) +
    
    # 图形美化
    theme_minimal() +
    labs(title = paste("扇形分割 -", length(sectors), "个扇形，每个1000个点"),
         x = "X", y = "Y", color = "扇形ID") +
    theme(legend.position = "right") +
    coord_fixed() +
    scale_color_viridis_d()
  
  return(p)
}

# 绘制结果
p1 = plot_sectors(sectors, center_x, center_y)

# 处理0°和360°边界的高级版本
split_into_sectors_advanced <- function(data, points_per_sector = 1000) {
  center_x <- mean(data$X)
  center_y <- mean(data$Y)
  
  delta_x <- data$X - center_x
  delta_y <- data$Y - center_y
  
  r <- sqrt(delta_x^2 + delta_y^2)
  theta_rad <- atan2(delta_y, delta_x)
  theta_deg <- (theta_rad * 180 / pi) %% 360
  
  data_polar <- data %>%
    mutate(
      r = r,
      theta_rad = theta_rad,
      theta_deg = theta_deg
    )
  
  # 按角度排序
  data_sorted <- data_polar %>% arrange(theta_deg)
  
  # 检查是否需要处理角度环绕
  first_angle <- data_sorted$theta_deg[1]
  last_angle <- data_sorted$theta_deg[nrow(data_sorted)]
  
  if (first_angle > 0 && last_angle < 360 && (360 - last_angle + first_angle) < 180) {
    cat("检测到角度环绕，进行调整...\n")
    
    # 找到环绕点（角度跳跃最大的地方）
    angle_diffs <- diff(data_sorted$theta_deg)
    wrap_idx <- which.max(angle_diffs) + 1
    
    if (wrap_idx > 1 && wrap_idx <= nrow(data_sorted)) {
      # 重新排列数据
      data_sorted <- rbind(data_sorted[wrap_idx:nrow(data_sorted), ], 
                           data_sorted[1:(wrap_idx-1), ])
    }
  }
  
  # 划分扇形
  n <- nrow(data_sorted)
  sectors <- list()
  
  for (i in seq(1, n, by = points_per_sector)) {
    end_idx <- min(i + points_per_sector - 1, n)
    sector_data <- data_sorted[i:end_idx, ]
    
    if (nrow(sector_data) == 0) next
    
    min_angle <- min(sector_data$theta_deg)
    max_angle <- max(sector_data$theta_deg)
    
    sector_info <- list(
      sector_id = length(sectors) + 1,
      points = sector_data %>% select(X, Y),
      min_angle = min_angle,
      max_angle = max_angle,
      center_angle = (min_angle + max_angle) / 2,
      point_count = nrow(sector_data)
    )
    
    sectors[[length(sectors) + 1]] <- sector_info
  }
  
  return(list(sectors = sectors, center_x = center_x, center_y = center_y))
}

# 使用高级版本
result_advanced <- split_into_sectors_advanced(data, 1000)
sectors_advanced <- result_advanced$sectors

# 绘制高级版本结果
p2 = plot_sectors(sectors_advanced, result_advanced$center_x, result_advanced$center_y)

# 扇形统计信息
sector_stats <- do.call(rbind, lapply(sectors, function(s) {
  data.frame(
    sector_id = s$sector_id,
    point_count = s$point_count,
    min_angle = s$min_angle,
    max_angle = s$max_angle,
    angle_span = s$max_angle - s$min_angle,
    avg_distance = mean(sqrt((s$points$X - center_x)^2 + (s$points$Y - center_y)^2))
  )
}))

print(sector_stats)

# 绘制扇形角度跨度分布
p3 = ggplot(sector_stats, aes(x = factor(sector_id), y = angle_span)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  labs(title = "每个扇形的角度跨度", x = "扇形ID", y = "角度跨度(度)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 绘制平均距离分布
p4 = ggplot(sector_stats, aes(x = factor(sector_id), y = avg_distance)) +
  geom_bar(stat = "identity", fill = "darkorange", alpha = 0.7) +
  labs(title = "每个扇形的平均距离", x = "扇形ID", y = "平均距离") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 创建非均匀分布的数据来测试
set.seed(123)
n_points <- 5000
angle_bias <- runif(n_points, 0, 2*pi)  # 创建角度偏差

data_nonuniform <- data.frame(
  X = 15 * cos(angle_bias) + rnorm(n_points, 0, 3),
  Y = 15 * sin(angle_bias) + rnorm(n_points, 0, 3)
)

# 对非均匀数据应用扇形分割
result_nonuniform <- split_into_sectors_advanced(data_nonuniform, 500)
p5 = plot_sectors(result_nonuniform$sectors, 
             result_nonuniform$center_x, 
             result_nonuniform$center_y)

print(p1+p2+p3+p4+p5)

