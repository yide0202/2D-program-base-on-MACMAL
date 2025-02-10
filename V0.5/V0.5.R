#### 科研1.0



####运行程序


library(Rcpp)
sourceCpp("new.cpp")
# 示例数据
deaths_r <- c(5, 3, 4, 2, 6, 7, 3, 4, 5, 2)
pop_r <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100)

# 运行模型（使用 AIC 和 Binomial 分布）
result <- runClusterSurvival(deaths_r, pop_r, criterion_type=1, dist_type_code=0)

# 查看结果
print(result)

# 可视化聚类概率
library(ggplot2)

N <- length(deaths_r)
time_points <- 1:N

data_plot <- data.frame(
  Time = time_points,
  Death_Rate = deaths_r / pop_r,
  p_i = result$p_i
)

ggplot(data_plot, aes(x = Time)) +
  geom_point(aes(y = Death_Rate), color = "grey", size = 2, alpha=0.6) +
  geom_line(aes(y = p_i), color = "blue", size = 1) +
  labs(title = "Cluster Survival MLE - 每个位置的聚类概率",
       x = "时间区间",
       y = "死亡率 / 聚类概率 p(i)") +
  theme_minimal()






# 加载必要的包
library(Rcpp)
library(ggplot2)
library(parallel)

# 加载 Cluster Survival MLE 的 C++ 代码
# 假设 C++ 代码已保存为 "ClusterSurvivalMLE.cpp"
# 请确保该文件与当前工作目录一致
sourceCpp("new.cpp")  # 请确保 "new.cpp" 中定义了 runClusterSurvival

# 示例数据生成
set.seed(123)
N <- 100  # 时间点数量
true_p <- c(rep(0.05, 30), rep(0.2, 40), rep(0.05, 30))  # 真正的死亡率
pop_r <- rep(100, N)  # 每个时间点的总风险人数
deaths_r <- rbinom(N, pop_r, true_p)  # 每个时间点的事件发生次数

# 运行 Cluster Survival MLE 模型（使用 AIC 和 Binomial 分布）
result <- runClusterSurvival(deaths_r, pop_r, criterion_type = 1, dist_type_code = 0)

# 查看模型结果
print(result)

# 引导法（Bootstrap）计算置信区间
# 定义引导函数
bootstrap_run <- function(deaths, pop, criterion_type, dist_type_code, n_boot = 1000) {
  # 创建集群
  cl <- makeCluster(detectCores() - 1)
  
  # 将必要的变量和函数传递给集群
  clusterExport(cl, c("deaths", "pop", "criterion_type", "dist_type_code"), envir = environment())
  
  # 在每个子进程中加载 C++ 函数
  clusterEvalQ(cl, {
    library(Rcpp)
    sourceCpp("new.cpp")  # 确保每个子进程都加载了 `runClusterSurvival`
  })
  
  # 并行运行引导
  results <- parLapply(cl, 1:n_boot, function(i) {
    sample_indices <- sample(1:length(deaths), replace = TRUE)
    sampled_deaths <- deaths[sample_indices]
    sampled_pop <- pop[sample_indices]
    runClusterSurvival(sampled_deaths, sampled_pop, criterion_type, dist_type_code)
  })
  
  # 关闭集群
  stopCluster(cl)
  
  return(results)
}

# 运行引导
set.seed(123)  # 设置随机种子以获得可重复结果
bootstrap_results <- bootstrap_run(deaths_r, pop_r, criterion_type = 1, dist_type_code = 0, n_boot = 1000)

# 提取每个位置的 p_i
p_i_matrix <- matrix(0, nrow = 1000, ncol = N)  # n_boot=1000

for (j in 1:length(bootstrap_results)) {
  p_i_matrix[j, ] <- bootstrap_results[[j]]$p_i
}

# 计算置信区间（2.5% 和 97.5% 分位数）
p_i_lower <- apply(p_i_matrix, 2, quantile, probs = 0.025)
p_i_upper <- apply(p_i_matrix, 2, quantile, probs = 0.975)

# 添加置信区间到原始结果中
result$p_i_lower <- p_i_lower
result$p_i_upper <- p_i_upper

# 提取聚类区间信息
clusters <- data.frame(
  cs = result$cs,
  ce = result$ce,
  p0 = result$p0,
  pc = result$pc
)

# 计算风险比（HR = pc / p0）
clusters$HR <- clusters$pc / clusters$p0

# 打印聚类区间及其风险比
print(clusters)

# 准备绘图数据
data_plot <- data.frame(
  Time = 1:N,
  Death_Rate = deaths_r / pop_r,
  p_i = result$p_i,
  p_i_lower = result$p_i_lower,
  p_i_upper = result$p_i_upper
)

# 创建综合性生存预测图
ggplot(data_plot, aes(x = Time)) +
  # 绘制实际死亡率点
  geom_point(aes(y = Death_Rate), color = "grey40", size = 2, alpha = 0.6) +
  # 绘制估计的聚类概率线
  geom_line(aes(y = p_i), color = "blue", size = 1) +
  # 添加置信区间带
  geom_ribbon(aes(ymin = p_i_lower, ymax = p_i_upper), alpha = 0.2, fill = "blue") +
  # 高亮显示聚类区间
  geom_rect(data = clusters, aes(xmin = cs, xmax = ce, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +
  # 添加聚类区间的 HR 信息
  geom_segment(data = clusters, aes(x = cs, xend = ce, y = p0, yend = p0),
               color = "darkgreen", size = 1) +
  geom_text(data = clusters, aes(x = (cs + ce) / 2, y = max(data_plot$p_i) + 0.02,
                                 label = sprintf("HR=%.2f", HR)),
            color = "black", size = 3, vjust = 0) +
  # 设置图形标签和主题
  labs(title = "Cluster Survival MLE - 生存预测与聚类概率",
       subtitle = "灰色点：实际死亡率；蓝色线：估计的聚类概率；红色区域：聚类区间",
       x = "时间点",
       y = "死亡率 / 聚类概率 p(i)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )









