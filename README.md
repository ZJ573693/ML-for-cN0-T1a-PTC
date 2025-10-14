# ML-for-cN0-PTMC
PMRA and Shapley-Based Machine Learning for Predicting Lymph Node Metastasis in Central Subregions of Clinically Node-Negative Papillary Thyroid Microcarcinoma: A Prospective Multicenter Validation and Development of a Web Calculator
#0数据处理及装包 

```{r}
install.packages("survival")
install.packages('rrtable')
install.packages('magrittr') 
install.packages("ggplot")
install.packages("dplyr")
install.packages("AER")
library("dplyr")
library("AER")
library(openxlsx) 
library(survival) 
library(rrtable)
library(ggplot2)

# 安装或更新必要的包
required_packages <- c("foreign", "rms", "car", "ggplot2")

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
  } else {
    update.packages(package, ask = FALSE)
  }
}

# 逐一检查并安装/更新包
sapply(required_packages, install_if_missing)

# 加载必要的包
library(foreign)
library(rms)
library(car)
library(ggplot2)


```
##0.1读数据
```{r}
# 读取数据
dataCP <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CPp.csv")
```

##0.2计算cut-off值

```{r}


# 加载所需的包
library(pROC)

# 计算ROC曲线
roc_obj <- roc(dataCP$Con.Paratracheal.LNM, dataCP$age)

# 根据最大Youden指数选择最佳cut-off值
best_cutoff <- coords(roc_obj, "best", ret="threshold")

# 打印最佳cut-off值
print(best_cutoff)
```

#1.总侧区的传统预测模型的建立
##1.1.1森林图之1:按or值排序的森林图
```{r}
library(car)
library(ggplot2)

# 读取数据
data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总编码后_插补.csv", check.names = FALSE)

# 设置列名
colnames(data1) <- c("Age", "Sex", "BMI", "Tumor border", "Aspect ratio", "Composition", "Internal echo pattern", "Internal echo homogeneous", "Calcification", "Tumor internal vascularization", "Tumor Peripheral blood flow", "Size", "Location", "Mulifocality", "Hashimoto", "Extrathyroidal extension", "Side of position", "Prelaryngeal LNM", "Pretracheal LNM", "Paratracheal LNM", "Con-Paratracheal LNM", "LNM-prRLN", "Total Central Lymph Node Metastasis", "age", "bmi", "size", "Prelaryngeal LNMR", "Prelaryngeal NLNM", "Pretracheal LNMR", "Pretracheal NLNM", "Paratracheal LNMR", "Paratracheal NLNM", "Con-Paratracheal LNMR", "Con-Paratracheal NLNM", "LNMR-prRLN", "NLNM-prRLN", "TCLNMR", "TCNLNM")

# 构建初始模型
initial_model <- glm(`Total Central Lymph Node Metastasis` ~ Age + Sex + `Tumor border` + `Aspect ratio` + `Internal echo pattern` + Calcification + `Tumor internal vascularization` + `Tumor Peripheral blood flow` + Size + Location + Mulifocality + `Extrathyroidal extension` + `Side of position`, data = data1, family = binomial())

# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Total Central Lymph Node Metastasis` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data1, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "black") +
  geom_point(aes(color = p_value < 0.05), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""), 
                x = -20, hjust = -0.1), size = 2.5) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1), size = 2.5) +
  coord_cartesian(xlim = c(-20, 20)) +
  scale_color_manual(values = c("black", "#BA3E45"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Total Central Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal()

# 保存图像函数
save_forest_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot(forest_plot, file.path(output_folder, "2.0.TCLNM-forest-1"))

print(forest_plot)

# 排序系数数据框
coef_df_sorted <- coef_df[order(coef_df$odds_ratio), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#BA3E45"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Total Central Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#BA3E45", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)

# 导出结果到CSV文件并反转顺序
write.csv(coef_df[nrow(coef_df):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.0.TCLNM-forest-2.csv", row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.0.TCLNM-forest-2.csv", row.names = FALSE)

# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.0.TCLNM-forest-2"))
print(forest_plot_sorted)


```

##1.1.2森林图之2:按输入顺序排列的森林图

```{r}
library(car)
library(ggplot2)

# 构建初始模型
initial_model <- glm(`Total Central Lymph Node Metastasis` ~Age+Sex+`Tumor border`+`Aspect ratio`+`Internal echo pattern` +Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Location+Mulifocality+`Extrathyroidal extension`+`Side of position`, data = data1, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Total Central Lymph Node Metastasis` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data1, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 计算95%置信区间
coef_df$LL <- coef_df$ci_lower
coef_df$UL <- coef_df$ci_upper

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 32.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#BA3E45"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Total Central Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df$p_value < 0.05,"#BA3E45","black")))

# 显示初始森林图
print(forest_plot)


coef_df_sorted <- coef_df[order(coef_df$variable), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#BA3E45"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Total Central Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#BA3E45", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)


# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.0.TCLNM-forest-1"))

# 保存CSV文件
write.csv(coef_df[nrow(coef_df):1, ], file = file.path(output_folder, "2.0.TCLNM-forest-1.csv"), row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = file.path(output_folder, "2.0.TCLNM-forest-1.csv"), row.names = FALSE)

print(forest_plot_sorted)


```

##1.2.1列线图以及验证曲线
```{r}
install.packages("foreign")
install.packages("rms")
library(foreign) 
library(rms)
library(car)
library(ggplot2)

# 读取数据
data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总编码后_插补.csv")

```

```{r}
data2$Age<-factor(data2$Age,levels = c(0,1),labels = c("Age≤45","Age>45"))
data2$Sex<-factor(data2$Sex,levels = c(0,1),labels = c("Female","Male"))
data2$BMI<-factor(data2$BMI,levels = c(0,1,2),labels = c("Underweight","Normal","Overweight"))

data2$Tumor.border<-factor(data2$Tumor.border,levels = c(0,1,2),labels = c("smooth or borderless","irregular shape or lsharpobed","extrandular invasion"))
data2$Aspect.ratio<-factor(data2$Aspect.ratio,levels = c(0,1),labels = c("≤1",">1"))
 data2$Composition<-factor(data2$Composition,levels = c(0,1,2),labels = c("cystic/cavernous","Mixed cystic and solid","solid"))
 data2$Internal.echo.pattern<-factor(data2$Internal.echo.pattern,levels = c(0,1,2,3),labels = c("echoless","high/isoechoic","hypoechoic","very hypoechoic"))
 data2$Internal.echo.homogeneous<-factor(data2$Internal.echo.homogeneous,levels = c(0,1),labels = c("Non-uniform","Uniform"))
 data2$Calcification<-factor(data2$Calcification,levels = c(0,1,2,3),labels = c("no or large comet tail", "coarse calcification","peripheral calcification","Microcalcification"))
data2$Tumor.internal.vascularization<-factor(data2$Tumor.internal.vascularization,levels = c(0,1),labels = c("Without","Abundant"))
data2$Tumor.Peripheral.blood.flow<-factor(data2$Tumor.Peripheral.blood.flow,levels = c(0,1),labels = c("Without","Abundant"))
data2$Size<-factor(data2$Size,levels = c(0,1),labels = c("≤5", ">5"))
data2$Location<-factor(data2$Location,levels = c(0,1),labels = c("Non-upper","Upper"))
data2$Mulifocality<-factor(data2$Mulifocality,levels = c(1,0),labels = c("Abundant", "Without"))
data2$Hashimoto<-factor(data2$Hashimoto,levels = c(1,0),labels = c("Abundant", "Without"))
data2$Extrathyroidal.extension<-factor(data2$Extrathyroidal.extension,levels = c(1,0),labels = c("Abundant", "Without"))
data2$Side.of.position<-factor(data2$Side.of.position,levels = c(0,1,2,3),labels = c("left","right","bilateral" ,"isthmus"))




data2$Prelaryngeal.LNM<-factor(data2$Prelaryngeal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data2$Pretracheal.LNM<-factor(data2$Pretracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data2$Paratracheal.LNM<-factor(data2$Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data2$Con.Paratracheal.LNM<-factor(data2$Con.Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data2$LNM.prRLN<-factor(data2$LNM.prRLN,levels = c(0,1),labels = c("No", "Yes"))
data2$Total.Central.Lymph.Node.Metastasis<-factor(data2$Total.Central.Lymph.Node.Metastasis,levels = c(0,1),labels = c("No", "Yes"))

```

```{r}
# 加载必要的包
library(rms)

# 准备数据
x <- as.data.frame(data2)
dd <- datadist(data2)
options(datadist = 'dd')

# 拟合逻辑回归模型并指定 x=TRUE 和 y=TRUE
fit1 <- lrm(Total.Central.Lymph.Node.Metastasis ~ Age + Sex + Tumor.border + Aspect.ratio + Calcification + Tumor.Peripheral.blood.flow + Size + Mulifocality + Extrathyroidal.extension, data = data2, x = TRUE, y = TRUE)

# 查看模型摘要
summary(fit1)

# 创建列线图
nom1 <- nomogram(fit1, fun = plogis, fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999), lp = FALSE, funlabel = "Total Central Lymph Node Metastasis")
plot(nom1)

# 验证曲线
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))

# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.0.TCLNM-nomogram.tiff", width = 10, height = 8, units = "in", res = 600, compression = "lzw")
plot(nom1)
dev.off()

# 保存验证曲线为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.0.TCLNM-calibration.tiff", width = 10, height = 8, units = "in", res = 600, compression = "lzw")
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
dev.off()


```

```{r}
# 改变尺寸的列线图
par(mar = c(1, 2, 2, 2))  # 调整绘图边距

# 创建 nomogram
nom2 <- nomogram(fit1, fun = plogis, fun.at = c(0.001, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.999),
                 lp = FALSE, funlabel="Total Central Lymph Node Metastasis")

# 绘制 nomogram
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)


# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.0.TCLNM-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)
dev.off()

```
##1.2.2传统预测模型的Roc曲线

```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Total.Central.Lymph.Node.Metastasis ~ Age + Sex + Tumor.border + Aspect.ratio + Calcification + Tumor.Peripheral.blood.flow + Size + Mulifocality + Extrathyroidal.extension,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Total.Central.Lymph.Node.Metastasis
test_response <- test_data$Total.Central.Lymph.Node.Metastasis
val_response1 <- val_data1$Total.Central.Lymph.Node.Metastasis
val_response2 <- val_data2$Total.Central.Lymph.Node.Metastasis
# 创建ROC对象
train_roc <- roc(train_response, train_probs)

test_roc <- roc(test_response, test_probs)
val_roc1 <- roc(val_response1, val_probs1)
val_roc2 <- roc(val_response2, val_probs2)

# 提取ROC曲线的坐标点
train_roc_data <- coords(train_roc, "all", ret = c("specificity", "sensitivity"))
test_roc_data <- coords(test_roc, "all", ret = c("specificity", "sensitivity"))
val_roc_data1 <- coords(val_roc1, "all", ret = c("specificity", "sensitivity"))
val_roc_data2 <- coords(val_roc2, "all", ret = c("specificity", "sensitivity"))

# 转换为数据框
train_roc_data <- as.data.frame(train_roc_data)
test_roc_data <- as.data.frame(test_roc_data)
val_roc_data1 <- as.data.frame(val_roc_data1)
val_roc_data2 <- as.data.frame(val_roc_data2)

# 绘制ROC曲线
roc_plot <- ggplot() +
  geom_line(data = train_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#9A4942", size = 0.6) +
  geom_line(data = test_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#BA3E45", size = 0.6) +
  geom_line(data = val_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#EABFBB", size = 0.6) +
  geom_line(data = val_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#EAB", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC for Total Central Lymph Node Metastasis Nomogram Prediction",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.4, label = paste("Train set AUC =", round(auc(train_roc), 3)), size = 4, color = "#9A4942")  +
  annotate("text", x = 0.7, y = 0.3, label = paste("Test set AUC =", round(auc(test_roc), 3)), size = 4, color = "#BA3E45")+
  annotate("text", x = 0.7, y = 0.2, label = paste("Validation set1 AUC =", round(auc(val_roc1), 3)), size = 4, color = "#EABFBB")+
  annotate("text", x = 0.7, y = 0.1, label = paste("Validation set2 AUC =", round(auc(val_roc2), 3)), size = 4, color = "#EAB")
# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.0.TCLNM-roc_curve.tiff", plot = roc_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")


```
##1.2.3传统预测模型的dca曲线
```{r}


# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Total.Central.Lymph.Node.Metastasis
test_response <- test_data$Total.Central.Lymph.Node.Metastasis
val_response1 <- val_data1$Total.Central.Lymph.Node.Metastasis
val_response2 <- val_data2$Total.Central.Lymph.Node.Metastasis


# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits <- sapply(thresholds, function(x) net_benefit(train_probs, train_response, x))
test_net_benefits <- sapply(thresholds, function(x) net_benefit(test_probs, test_response, x))
val_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs1, val_response1, x))
val_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs2, val_response2, x))


# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(val_response1)), val_response1, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 找到最大净收益点
train_max_nb <- max(train_net_benefits)
train_max_nb_threshold <- thresholds[which.max(train_net_benefits)]
test_max_nb <- max(test_net_benefits)
test_max_nb_threshold <- thresholds[which.max(test_net_benefits)]
val_max_nb1 <- max(val_net_benefits1)
val_max_nb_threshold1 <- thresholds[which.max(val_net_benefits1)]
val_max_nb2 <- max(val_net_benefits2)
val_max_nb_threshold2 <- thresholds[which.max(val_net_benefits2)]




# 绘制DCA曲线
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit = train_net_benefits,
  test_net_benefit = test_net_benefits,
  val_net_benefit1 = val_net_benefits1,
  val_net_benefit2 = val_net_benefits2,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data, aes(x = threshold)) +
  geom_line(aes(y = train_net_benefit, color = "Train set"), size = 0.6) +
  geom_line(aes(y = test_net_benefit, color = "Test set"), size = 0.6) +
  geom_line(aes(y = val_net_benefit1, color = "Validation set1"), size = 0.6) +
  geom_line(aes(y = val_net_benefit2, color = "Validation set2"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Total Central Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(values = c("Train set" = "#9A4942", "Test set" = "#BA3E45", "Validation set1" = "#EABFBB", "Validation set2" = "#EAB","All" = "grey", "None" = "black")) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotate("text", x = 0.4, y = 0.15, label = "Train set", size = 4, color = "#9A4942") +
  annotate("text", x = 0.4, y = 0.1, label = "Test set", size = 4, color = "#BA3E45") +
  annotate("text", x = 0.4, y = 0.01, label = "Validation set1", size = 4, color = "#EABFBB") +
  annotate("text", x = 0.4, y = 0.05, label = "Validation set2", size = 4, color = "#EAB") +
  annotate("text", x = train_max_nb_threshold, y = train_max_nb, label = sprintf("Max: %.3f", train_max_nb), color = "#9A4942", vjust = -1) +
  annotate("text", x = test_max_nb_threshold, y = test_max_nb, label = sprintf("Max: %.3f", test_max_nb), color = "#BA3E45", vjust = -1) +
   annotate("text", x = val_max_nb_threshold1, y = val_max_nb1, label = sprintf("Max: %.3f", val_max_nb1), color = "#EABFBB", vjust = -1) +
   annotate("text", x = val_max_nb_threshold2, y = val_max_nb2, label = sprintf("Max: %.3f", val_max_nb2), color = "#EAB", vjust = -1) +
  coord_cartesian(ylim = c(-0.05, 0.4), xlim = c(0, 0.8))


# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.0.TCLNM-dca_curve.tiff", plot = dca_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")

print(dca_plot)

```

##1.2.4 保存胜率概率
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- lrm(Total.Central.Lymph.Node.Metastasis ~ Age + Sex + Tumor.border + Aspect.ratio + Calcification + Tumor.Peripheral.blood.flow + Size + Mulifocality + Extrathyroidal.extension,
            data = tra_data,  x = TRUE, y = TRUE)

#删掉了一些
nom1 <- predict(fit1, type = "fitted")

# 导出预测结果
nomogram_predictions <- data.frame(nomogram_prediction = nom1)
write.csv(nomogram_predictions, '/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/5.对比/1.TCLNM.nomogram_predictions.csv', row.names = FALSE)

```









#2.喉前的传统预测模型的建立
##2.1.1森林图之1:按or值排序的森林图
```{r}
library(car)
library(ggplot2)


# 读取数据
data21 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补.csv", check.names = FALSE)

# 设置列名
colnames(data21) <- c("Age","Sex","BMI","Tumor border","Aspect ratio","Composition","Internal echo pattern","Internal echo homogeneous","Calcification","Tumor internal vascularization","Tumor Peripheral blood flow","Size","Location","Mulifocality","Hashimoto","Extrathyroidal extension","Side of position","Prelaryngeal LNM","Pretracheal LNM","Paratracheal LNM","Con-Paratracheal LNM","LNM-prRLN","Total Central Lymph Node Metastasis","age","bmi","size","Prelaryngeal LNMR","Prelaryngeal NLNM","Pretracheal LNMR","Pretracheal NLNM","Paratracheal LNMR","Paratracheal NLNM","Con-Paratracheal LNMR","Con-Paratracheal NLNM","LNMR-prRLN","NLNM-prRLN","TCLNMR","TCNLNM")



# 构建初始模型
initial_model <- glm(`Prelaryngeal LNM` ~Age+Sex+`Tumor border`+`Internal echo pattern` +Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Location+Mulifocality+Hashimoto+`Extrathyroidal extension`+`Pretracheal LNM`+`Paratracheal LNM`+`Con-Paratracheal LNM`+`LNM-prRLN`, data = data21, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Prelaryngeal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data21, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "black") +
  geom_point(aes(color = p_value < 0.05), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""), 
                x = -20, hjust = -0.1), size = 2.5) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1), size = 2.5) +
  coord_cartesian(xlim = c(-20, 20)) +
  scale_color_manual(values = c("black", "#D2431C"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Prelaryngeal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal()

# 保存图像函数
save_forest_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)

}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot(forest_plot, file.path(output_folder, "2.1.Prelaryngeal-forest-1"))

print(forest_plot)



coef_df_sorted <- coef_df[order(coef_df$odds_ratio), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#D2431C"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Prelaryngeal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#D2431C", "black")))


# 显示排序后的森林图
print(forest_plot_sorted)

# 导出结果到CSV文件并反转顺序
write.csv(coef_df[nrow(coef_df):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.1.Prelaryngeal-forest-2.csv", row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.1.Prelaryngeal-forest-2.csv", row.names = FALSE)

# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.1.Prelaryngeal-forest-2"))
print(forest_plot_sorted)
```

##2.1.2森林图之2:按输入顺序排列的森林图

```{r}
library(car)
library(ggplot2)

# 构建初始模型
initial_model <- glm(`Prelaryngeal LNM` ~Age+Sex+`Tumor border`+`Internal echo pattern` +Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Location+Mulifocality+Hashimoto+`Extrathyroidal extension`+`Pretracheal LNM`+`Paratracheal LNM`+`Con-Paratracheal LNM`+`LNM-prRLN`, data = data21, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 10])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Prelaryngeal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data21, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 计算95%置信区间
coef_df$LL <- coef_df$ci_lower
coef_df$UL <- coef_df$ci_upper

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 32.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#D2431C"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Prelaryngeal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df$p_value < 0.05,"#D2431C","black")))

# 显示初始森林图
print(forest_plot)


coef_df_sorted <- coef_df[order(coef_df$variable), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#D2431C"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Prelaryngeal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#D2431C", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)


# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.1.Prelaryngeal-forest-1"))

# 保存CSV文件
write.csv(coef_df[nrow(coef_df):1, ], file = file.path(output_folder, "2.1.Prelaryngeal-forest-1.csv"), row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = file.path(output_folder, "2.1.Prelaryngeal-forest-1.csv"), row.names = FALSE)


print(forest_plot_sorted)


```

##2.2.1列线图以及验证曲线
```{r}
# 读取数据
data21 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补.csv")

```

```{r}
data21$Age<-factor(data21$Age,levels = c(0,1),labels = c("Age≤45","Age>45"))
data21$Sex<-factor(data21$Sex,levels = c(0,1),labels = c("Female","Male"))
data21$BMI<-factor(data21$BMI,levels = c(0,1,2),labels = c("Underweight","Normal","Overweight"))

data21$Tumor.border<-factor(data21$Tumor.border,levels = c(0,1,2),labels = c("smooth or borderless","irregular shape or lsharpobed","extrandular invasion"))
data21$Aspect.ratio<-factor(data21$Aspect.ratio,levels = c(0,1),labels = c("≤1",">1"))
 data21$Composition<-factor(data21$Composition,levels = c(0,1,2),labels = c("cystic/cavernous","Mixed cystic and solid","solid"))
 data21$Internal.echo.pattern<-factor(data21$Internal.echo.pattern,levels = c(0,1,2,3),labels = c("echoless","high/isoechoic","hypoechoic","very hypoechoic"))
 data21$Internal.echo.homogeneous<-factor(data21$Internal.echo.homogeneous,levels = c(0,1),labels = c("Non-uniform","Uniform"))
 data21$Calcification<-factor(data21$Calcification,levels = c(0,1,2,3),labels = c("no or large comet tail", "coarse calcification","peripheral calcification","Microcalcification"))
data21$Tumor.internal.vascularization<-factor(data21$Tumor.internal.vascularization,levels = c(0,1),labels = c("Without","Abundant"))
data21$Tumor.Peripheral.blood.flow<-factor(data21$Tumor.Peripheral.blood.flow,levels = c(0,1),labels = c("Without","Abundant"))
data21$Size<-factor(data21$Size,levels = c(0,1),labels = c("≤5", ">5"))
data21$Location<-factor(data21$Location,levels = c(0,1),labels = c("Non-upper","Upper"))
data21$Mulifocality<-factor(data21$Mulifocality,levels = c(1,0),labels = c("Abundant", "Without"))
data21$Hashimoto<-factor(data21$Hashimoto,levels = c(1,0),labels = c("Abundant", "Without"))
data21$Extrathyroidal.extension<-factor(data21$Extrathyroidal.extension,levels = c(1,0),labels = c("Abundant", "Without"))
data21$Side.of.position<-factor(data21$Side.of.position,levels = c(0,1,2,3),labels = c("left","right","bilateral" ,"isthmus"))




data21$Prelaryngeal.LNM<-factor(data21$Prelaryngeal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data21$Pretracheal.LNM<-factor(data21$Pretracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data21$Paratracheal.LNM<-factor(data21$Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data21$Con.Paratracheal.LNM<-factor(data21$Con.Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data21$LNM.prRLN<-factor(data21$LNM.prRLN,levels = c(0,1),labels = c("No", "Yes"))
data21$Total.Central.Lymph.Node.Metastasis<-factor(data21$Total.Central.Lymph.Node.Metastasis,levels = c(0,1),labels = c("No", "Yes"))

```

```{r}
# 加载必要的包
library(rms)

# 准备数据
x <- as.data.frame(data21)
dd <- datadist(data21)
options(datadist = 'dd')

# 拟合逻辑回归模型并指定 x=TRUE 和 y=TRUE
fit1 <- lrm(Prelaryngeal.LNM ~Location+Hashimoto+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN, data = data21, x = TRUE, y = TRUE)

# 查看模型摘要
summary(fit1)

# 创建列线图
nom1 <- nomogram(fit1, fun = plogis, fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999), lp = FALSE, funlabel = "Prelaryngeal Lymph Node Metastasis")
plot(nom1)

# 验证曲线
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))

# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.1.Prelaryngeal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom1)
dev.off()

# 保存验证曲线为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.1.Prelaryngeal-calibration.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
dev.off()


```

```{r}
# 改变尺寸的列线图
par(mar = c(1, 2, 2, 2))  # 调整绘图边距

# 创建 nomogram
nom2 <- nomogram(fit1, fun = plogis, fun.at = c(0.001, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.999),
                 lp = FALSE, funlabel="Prelaryngeal Lymph Node Metastasis")

# 绘制 nomogram
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)


# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.1.Prelaryngeal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)
dev.off()

```
##2.2.2传统预测模型的Roc曲线

```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1H编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2H编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Prelaryngeal.LNM ~Location+Hashimoto+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Prelaryngeal.LNM
test_response <- test_data$Prelaryngeal.LNM
val_response1 <- val_data1$Prelaryngeal.LNM
val_response2 <- val_data2$Prelaryngeal.LNM
# 创建ROC对象
train_roc <- roc(train_response, train_probs)
test_roc <- roc(test_response, test_probs)
val_roc1 <- roc(val_response1, val_probs1)
val_roc2 <- roc(val_response2, val_probs2)

# 提取ROC曲线的坐标点
train_roc_data <- coords(train_roc, "all", ret = c("specificity", "sensitivity"))
test_roc_data <- coords(test_roc, "all", ret = c("specificity", "sensitivity"))
val_roc_data1 <- coords(val_roc1, "all", ret = c("specificity", "sensitivity"))
val_roc_data2 <- coords(val_roc2, "all", ret = c("specificity", "sensitivity"))

# 转换为数据框
train_roc_data <- as.data.frame(train_roc_data)
test_roc_data <- as.data.frame(test_roc_data)
val_roc_data1 <- as.data.frame(val_roc_data1)
val_roc_data2 <- as.data.frame(val_roc_data2)

# 绘制ROC曲线
roc_plot <- ggplot() +
  geom_line(data = train_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#BB431C", size = 0.6) +
  geom_line(data = test_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#D2431C", size = 0.6) +
  geom_line(data = val_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#F2AB6A", size = 0.6) +
  geom_line(data = val_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#F5D18B", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC for Prelaryngeal Lymph Node Metastasis Nomogram Prediction",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.4, label = paste("Train set AUC =", round(auc(train_roc), 3)), size = 4, color = "#BB431C")  +
  annotate("text", x = 0.7, y = 0.3, label = paste("Test set AUC =", round(auc(test_roc), 3)), size = 4, color = "#D2431C")+
  annotate("text", x = 0.7, y = 0.2, label = paste("Validation set1 AUC =", round(auc(val_roc1), 3)), size = 4, color = "#F2AB6A")+
  annotate("text", x = 0.7, y = 0.1, label = paste("Validation set2 AUC =", round(auc(val_roc2), 3)), size = 4, color = "#F5D18B")
# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.1.Prelaryngeal-roc_curve.tiff", plot = roc_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")


```
##2.2.3传统预测模型的dca曲线
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补2为了dca.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1H编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2H编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.6  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Prelaryngeal.LNM ~Size+Location+Hashimoto+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Prelaryngeal.LNM
test_response <- test_data$Prelaryngeal.LNM
val_response1 <- val_data1$Prelaryngeal.LNM
val_response2 <- val_data2$Prelaryngeal.LNM


# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits <- sapply(thresholds, function(x) net_benefit(train_probs, train_response, x))
test_net_benefits <- sapply(thresholds, function(x) net_benefit(test_probs, test_response, x))
val_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs1, val_response1, x))
val_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs2, val_response2, x))


# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(val_response1)), val_response1, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 找到最大净收益点
train_max_nb <- max(train_net_benefits)
train_max_nb_threshold <- thresholds[which.max(train_net_benefits)]
test_max_nb <- max(test_net_benefits)
test_max_nb_threshold <- thresholds[which.max(test_net_benefits)]
val_max_nb1 <- max(val_net_benefits1)
val_max_nb_threshold1 <- thresholds[which.max(val_net_benefits1)]
val_max_nb2 <- max(val_net_benefits2)
val_max_nb_threshold2 <- thresholds[which.max(val_net_benefits2)]




# 绘制DCA曲线
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit = train_net_benefits,
  test_net_benefit = test_net_benefits,
  val_net_benefit1 = val_net_benefits1,
  val_net_benefit2 = val_net_benefits2,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data, aes(x = threshold)) +
  geom_line(aes(y = train_net_benefit, color = "Train set"), size = 0.6) +
  geom_line(aes(y = test_net_benefit, color = "Test set"), size = 0.6) +
  geom_line(aes(y = val_net_benefit1, color = "Validation set1"), size = 0.6) +
  geom_line(aes(y = val_net_benefit2, color = "Validation set2"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Prelaryngeal Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(values = c("Train set" = "#BB431C", "Test set" = "#D2431C", "Validation set1" = "#F2AB6A", "Validation set2" = "#F5D18B","All" = "grey", "None" = "black")) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotate("text", x = 0.2, y = 0.02, label = "Train set", size = 4, color = "#BB431C") +
  annotate("text", x = 0.2, y = 0.01, label = "Test set", size = 4, color = "#D2431C") +
  annotate("text", x = 0.2, y = 0.06, label = "Validation set1", size = 4, color = "#F2AB6A") +
  annotate("text", x = 0.2, y = 0.04, label = "Validation set2", size = 4, color = "#F5D18B") +
  annotate("text", x = train_max_nb_threshold, y = train_max_nb, label = sprintf("Max: %.3f", train_max_nb), color = "#BB431C", vjust = -1) +
  annotate("text", x = test_max_nb_threshold, y = test_max_nb, label = sprintf("Max: %.3f", test_max_nb), color = "#D2431C", vjust = -1) +
   annotate("text", x = val_max_nb_threshold1, y = val_max_nb1, label = sprintf("Max: %.3f", val_max_nb1), color = "#F2AB6A", vjust = -1) +
   annotate("text", x = val_max_nb_threshold2, y = val_max_nb2, label = sprintf("Max: %.3f", val_max_nb2), color = "#F5D18B", vjust = -1) +
  coord_cartesian(ylim = c(-0.05, 0.15), xlim = c(0, 0.6))


# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.1.Prelaryngeal-dca_curve.tiff", plot = dca_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")

print(dca_plot)

```
##2.2.4 保存胜率概率
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1H编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2H编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")



# 构建模型
fit2 <- lrm(Prelaryngeal.LNM ~Location+Hashimoto+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data,  x = TRUE, y = TRUE)

#删掉了一些
nom2 <- predict(fit2, type = "fitted")

# 导出预测结果
nomogram_predictions2 <- data.frame(nomogram_prediction = nom2)
write.csv(nomogram_predictions2, '/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/5.对比/2.Prelaryngeal.LNM.nomogram_predictions.csv', row.names = FALSE)

```

#3.气管前的传统预测模型的建立
##3.1.1森林图之1:按or值排序的森林图
```{r}
library(car)
library(ggplot2)


# 读取数据
data22 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv", check.names = FALSE)

# 设置列名
colnames(data22) <- c("Age","Sex","BMI","Tumor border","Aspect ratio","Composition","Internal echo pattern","Internal echo homogeneous","Calcification","Tumor internal vascularization","Tumor Peripheral blood flow","Size","Location","Mulifocality","Hashimoto","Extrathyroidal extension","Side of position","Prelaryngeal LNM","Pretracheal LNM","Paratracheal LNM","Con-Paratracheal LNM","LNM-prRLN","Total Central Lymph Node Metastasis","age","bmi","size","Prelaryngeal LNMR","Prelaryngeal NLNM","Pretracheal LNMR","Pretracheal NLNM","Paratracheal LNMR","Paratracheal NLNM","Con-Paratracheal LNMR","Con-Paratracheal NLNM","LNMR-prRLN","NLNM-prRLN","TCLNMR","TCNLNM")



# 构建初始模型
initial_model <- glm(`Pretracheal LNM` ~Age+Sex+`Tumor border`+`Internal echo pattern` +`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Mulifocality+`Extrathyroidal extension`+`Side of position`+`Prelaryngeal LNM`+`Paratracheal LNM`+`Con-Paratracheal LNM`+`LNM-prRLN`, data = data22, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Pretracheal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data22, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "black") +
  geom_point(aes(color = p_value < 0.05), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""), 
                x = -20, hjust = -0.1), size = 2.5) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1), size = 2.5) +
  coord_cartesian(xlim = c(-20, 20)) +
  scale_color_manual(values = c("black", "#ECAC27"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Pretracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal()

# 保存图像函数
save_forest_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)

}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot(forest_plot, file.path(output_folder, "2.2.Pretracheal-forest-1"))

print(forest_plot)



coef_df_sorted <- coef_df[order(coef_df$odds_ratio), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#ECAC27"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Pretracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#ECAC27", "black")))


# 显示排序后的森林图
print(forest_plot_sorted)

# 导出结果到CSV文件并反转顺序
write.csv(coef_df[nrow(coef_df):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.2.Pretracheal-forest-2.csv", row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.2.Pretracheal-forest-2.csv", row.names = FALSE)

# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.2.Pretracheal-forest-2"))
print(forest_plot_sorted)
```

##3.1.2森林图之2:按输入顺序排列的森林图

```{r}
library(car)
library(ggplot2)

# 构建初始模型
initial_model <- glm(`Pretracheal LNM` ~Age+Sex+`Tumor border`+`Internal echo pattern` +`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Mulifocality+`Extrathyroidal extension`+`Side of position`+`Prelaryngeal LNM`+`Paratracheal LNM`+`Con-Paratracheal LNM`+`LNM-prRLN`, data = data22, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Pretracheal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data22, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 计算95%置信区间
coef_df$LL <- coef_df$ci_lower
coef_df$UL <- coef_df$ci_upper

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 32.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#ECAC27"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Pretracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df$p_value < 0.05,"#ECAC27","black")))

# 显示初始森林图
print(forest_plot)


coef_df_sorted <- coef_df[order(coef_df$variable), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#ECAC27"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Pretracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#ECAC27", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)


# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.2.Pretracheal-forest-1"))

# 保存CSV文件
write.csv(coef_df[nrow(coef_df):1, ], file = file.path(output_folder, "2.2.Pretracheal-forest-1.csv"), row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = file.path(output_folder, "2.2.Pretracheal-forest-1.csv"), row.names = FALSE)


print(forest_plot_sorted)


```

##3.2.1列线图以及验证曲线
```{r}
# 读取数据
data22 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv")

```

```{r}
data22$Age<-factor(data22$Age,levels = c(0,1),labels = c("Age≤45","Age>45"))
data22$Sex<-factor(data22$Sex,levels = c(0,1),labels = c("Female","Male"))
data22$BMI<-factor(data22$BMI,levels = c(0,1,2),labels = c("Underweight","Normal","Overweight"))

data22$Tumor.border<-factor(data22$Tumor.border,levels = c(0,1,2),labels = c("smooth or borderless","irregular shape or lsharpobed","extrandular invasion"))
data22$Aspect.ratio<-factor(data22$Aspect.ratio,levels = c(0,1),labels = c("≤1",">1"))
 data22$Composition<-factor(data22$Composition,levels = c(0,1,2),labels = c("cystic/cavernous","Mixed cystic and solid","solid"))
 data22$Internal.echo.pattern<-factor(data22$Internal.echo.pattern,levels = c(0,1,2,3),labels = c("echoless","high/isoechoic","hypoechoic","very hypoechoic"))
 data22$Internal.echo.homogeneous<-factor(data22$Internal.echo.homogeneous,levels = c(0,1),labels = c("Non-uniform","Uniform"))
 data22$Calcification<-factor(data22$Calcification,levels = c(0,1,2,3),labels = c("no or large comet tail", "coarse calcification","peripheral calcification","Microcalcification"))
data22$Tumor.internal.vascularization<-factor(data22$Tumor.internal.vascularization,levels = c(0,1),labels = c("Without","Abundant"))
data22$Tumor.Peripheral.blood.flow<-factor(data22$Tumor.Peripheral.blood.flow,levels = c(0,1),labels = c("Without","Abundant"))
data22$Size<-factor(data22$Size,levels = c(0,1),labels = c("≤5", ">5"))
data22$Location<-factor(data22$Location,levels = c(0,1),labels = c("Non-upper","Upper"))
data22$Mulifocality<-factor(data22$Mulifocality,levels = c(1,0),labels = c("Abundant", "Without"))
data22$Hashimoto<-factor(data22$Hashimoto,levels = c(1,0),labels = c("Abundant", "Without"))
data22$Extrathyroidal.extension<-factor(data22$Extrathyroidal.extension,levels = c(1,0),labels = c("Abundant", "Without"))
data22$Side.of.position<-factor(data22$Side.of.position,levels = c(0,1,2,3),labels = c("left","right","bilateral" ,"isthmus"))




data22$Prelaryngeal.LNM<-factor(data22$Prelaryngeal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data22$Pretracheal.LNM<-factor(data22$Pretracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data22$Paratracheal.LNM<-factor(data22$Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data22$Con.Paratracheal.LNM<-factor(data22$Con.Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data22$LNM.prRLN<-factor(data22$LNM.prRLN,levels = c(0,1),labels = c("No", "Yes"))
data22$Total.Central.Lymph.Node.Metastasis<-factor(data22$Total.Central.Lymph.Node.Metastasis,levels = c(0,1),labels = c("No", "Yes"))

```

```{r}
# 加载必要的包
library(rms)

# 准备数据
x <- as.data.frame(data22)
dd <- datadist(data22)
options(datadist = 'dd')

# 拟合逻辑回归模型并指定 x=TRUE 和 y=TRUE
fit1 <- lrm(Pretracheal.LNM ~Age+Sex+Tumor.Peripheral.blood.flow+Mulifocality+Prelaryngeal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN, data = data22, x = TRUE, y = TRUE)

# 查看模型摘要
summary(fit1)

# 创建列线图
nom1 <- nomogram(fit1, fun = plogis, fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999), lp = FALSE, funlabel = "Pretracheal Lymph Node Metastasis")
plot(nom1)

# 验证曲线
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))

# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.2.Pretracheal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom1)
dev.off()

# 保存验证曲线为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.2.Pretracheal-calibration.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
dev.off()


```

```{r}
# 改变尺寸的列线图
par(mar = c(1, 2, 2, 2))  # 调整绘图边距

# 创建 nomogram
nom2 <- nomogram(fit1, fun = plogis, fun.at = c(0.001, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.999),
                 lp = FALSE, funlabel="Pretracheal Lymph Node Metastasis")

# 绘制 nomogram
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)


# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.2.Pretracheal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)
dev.off()

```
##3.2.2传统预测模型的Roc曲线

```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1Q编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2Q编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Pretracheal.LNM ~Age+Sex+Tumor.Peripheral.blood.flow+Mulifocality+Prelaryngeal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Pretracheal.LNM
test_response <- test_data$Pretracheal.LNM
val_response1 <- val_data1$Pretracheal.LNM
val_response2 <- val_data2$Pretracheal.LNM
# 创建ROC对象
train_roc <- roc(train_response, train_probs)
test_roc <- roc(test_response, test_probs)
val_roc1 <- roc(val_response1, val_probs1)
val_roc2 <- roc(val_response2, val_probs2)

# 提取ROC曲线的坐标点
train_roc_data <- coords(train_roc, "all", ret = c("specificity", "sensitivity"))
test_roc_data <- coords(test_roc, "all", ret = c("specificity", "sensitivity"))
val_roc_data1 <- coords(val_roc1, "all", ret = c("specificity", "sensitivity"))
val_roc_data2 <- coords(val_roc2, "all", ret = c("specificity", "sensitivity"))

# 转换为数据框
train_roc_data <- as.data.frame(train_roc_data)
test_roc_data <- as.data.frame(test_roc_data)
val_roc_data1 <- as.data.frame(val_roc_data1)
val_roc_data2 <- as.data.frame(val_roc_data2)

# 绘制ROC曲线
roc_plot <- ggplot() +
  geom_line(data = train_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#C9A51A", size = 0.6) +
  geom_line(data = test_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#ECAC27", size = 0.6) +
  geom_line(data = val_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#EDDE23", size = 0.6) +
  geom_line(data = val_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#FFFF66", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC for Pretracheal Lymph Node Metastasis Nomogram Prediction",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.4, label = paste("Train set AUC =", round(auc(train_roc), 3)), size = 4, color = "#C9A51A")  +
  annotate("text", x = 0.7, y = 0.3, label = paste("Test set AUC =", round(auc(test_roc), 3)), size = 4, color = "#ECAC27")+
  annotate("text", x = 0.7, y = 0.2, label = paste("Validation set1 AUC =", round(auc(val_roc1), 3)), size = 4, color = "#EDDE23")+
  annotate("text", x = 0.7, y = 0.1, label = paste("Validation set2 AUC =", round(auc(val_roc2), 3)), size = 4, color = "#FFFF66")
# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.2.Pretracheal-roc_curve.tiff", plot = roc_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")


```
##3.2.3传统预测模型的dca曲线
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1Q编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2Q编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Pretracheal.LNM ~Age+Sex+Tumor.Peripheral.blood.flow+Size+Mulifocality+Prelaryngeal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())


# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Pretracheal.LNM
test_response <- test_data$Pretracheal.LNM
val_response1 <- val_data1$Pretracheal.LNM
val_response2 <- val_data2$Pretracheal.LNM
# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits <- sapply(thresholds, function(x) net_benefit(train_probs, train_response, x))
test_net_benefits <- sapply(thresholds, function(x) net_benefit(test_probs, test_response, x))
val_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs1, val_response1, x))
val_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs2, val_response2, x))


# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(train_response)), train_response, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 找到最大净收益点
train_max_nb <- max(train_net_benefits)
train_max_nb_threshold <- thresholds[which.max(train_net_benefits)]
test_max_nb <- max(test_net_benefits)
test_max_nb_threshold <- thresholds[which.max(test_net_benefits)]
val_max_nb1 <- max(val_net_benefits1)
val_max_nb_threshold1 <- thresholds[which.max(val_net_benefits1)]
val_max_nb2 <- max(val_net_benefits2)
val_max_nb_threshold2 <- thresholds[which.max(val_net_benefits2)]




# 绘制DCA曲线
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit = train_net_benefits,
  test_net_benefit = test_net_benefits,
  val_net_benefit1 = val_net_benefits1,
  val_net_benefit2 = val_net_benefits2,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data, aes(x = threshold)) +
  geom_line(aes(y = train_net_benefit, color = "Train set"), size = 0.6) +
  geom_line(aes(y = test_net_benefit, color = "Test set"), size = 0.6) +
  geom_line(aes(y = val_net_benefit1, color = "Validation set1"), size = 0.6) +
  geom_line(aes(y = val_net_benefit2, color = "Validation set2"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Pretracheal Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(values = c("Train set" = "#C9A51A", "Test set" = "#ECAC27", "Validation set1" = "#EDDE23", "Validation set2" = "#FFFF66","All" = "grey", "None" = "black")) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotate("text", x = 0.3, y = 0.12, label = "Train set", size = 4, color = "#C9A51A") +
  annotate("text", x = 0.3, y = 0.09, label = "Test set", size = 4, color = "#ECAC27") +
  annotate("text", x = 0.3, y = 0.2, label = "Validation set1", size = 4, color = "#EDDE23") +
  annotate("text", x = 0.3, y = 0.15, label = "Validation set2", size = 4, color = "#FFFF66") +
  annotate("text", x = train_max_nb_threshold, y = train_max_nb, label = sprintf("Max: %.3f", train_max_nb), color = "#C9A51A", vjust = -1) +
  annotate("text", x = test_max_nb_threshold, y = test_max_nb, label = sprintf("Max: %.3f", test_max_nb), color = "#ECAC27", vjust = -1) +
   annotate("text", x = val_max_nb_threshold1, y = val_max_nb1, label = sprintf("Max: %.3f", val_max_nb1), color = "#EDDE23", vjust = -1) +
   annotate("text", x = val_max_nb_threshold2, y = val_max_nb2, label = sprintf("Max: %.3f", val_max_nb2), color = "#FFFF66", vjust = -1) +
  coord_cartesian(ylim = c(-0.05, 0.4), xlim = c(0, 0.5))


# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.2.Pretracheal-dca_curve.tiff", plot = dca_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")

print(dca_plot)

```
##3.2.4 保存胜率概率
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1Q编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2Q编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")




# 构建模型
fit3 <- lrm(Pretracheal.LNM ~Age+Sex+Tumor.Peripheral.blood.flow+Mulifocality+Prelaryngeal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data,  x = TRUE, y = TRUE)

#删掉了一些
nom3 <- predict(fit3, type = "fitted")

# 导出预测结果
nomogram_predictions3 <- data.frame(nomogram_prediction = nom3)
write.csv(nomogram_predictions3, '/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/5.对比/3.Pretracheal.LNM.nomogram_predictions.csv', row.names = FALSE)

```

#4.同侧气管旁的传统预测模型的建立
##4.1.1森林图之1:按or值排序的森林图
```{r}
library(car)
library(ggplot2)


# 读取数据
data23 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv", check.names = FALSE)

# 设置列名
colnames(data23) <- c("Age","Sex","BMI","Tumor border","Aspect ratio","Composition","Internal echo pattern","Internal echo homogeneous","Calcification","Tumor internal vascularization","Tumor Peripheral blood flow","Size","Location","Mulifocality","Hashimoto","Extrathyroidal extension","Side of position","Prelaryngeal LNM","Pretracheal LNM","Paratracheal LNM","Con-Paratracheal LNM","LNM-prRLN","Total Central Lymph Node Metastasis","age","bmi","size","Prelaryngeal LNMR","Prelaryngeal NLNM","Pretracheal LNMR","Pretracheal NLNM","Paratracheal LNMR","Paratracheal NLNM","Con-Paratracheal LNMR","Con-Paratracheal NLNM","LNMR-prRLN","NLNM-prRLN","TCLNMR","TCNLNM")



# 构建初始模型
initial_model <- glm(`Paratracheal LNM` ~Age+Sex+`Tumor border`+`Aspect ratio`+`Internal echo pattern` +Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Mulifocality+`Extrathyroidal extension`+`Side of position`+`Prelaryngeal LNM`+`Pretracheal LNM`+`Con-Paratracheal LNM`+`LNM-prRLN`, data = data23, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Paratracheal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data23, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "black") +
  geom_point(aes(color = p_value < 0.05), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""), 
                x = -20, hjust = -0.1), size = 2.5) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1), size = 2.5) +
  coord_cartesian(xlim = c(-20, 20)) +
  scale_color_manual(values = c("black", "#79902D"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal()

# 保存图像函数
save_forest_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)

}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot(forest_plot, file.path(output_folder, "2.3.Paratracheal-forest-1"))

print(forest_plot)



coef_df_sorted <- coef_df[order(coef_df$odds_ratio), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#79902D"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#79902D", "black")))


# 显示排序后的森林图
print(forest_plot_sorted)

# 导出结果到CSV文件并反转顺序
write.csv(coef_df[nrow(coef_df):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.3.Paratracheal-forest-2.csv", row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.3.Paratracheal-forest-2.csv", row.names = FALSE)

# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.3.Paratracheal-forest-2"))
print(forest_plot_sorted)
```

##4.1.2森林图之2:按输入顺序排列的森林图

```{r}
library(car)
library(ggplot2)

# 构建初始模型
initial_model <- glm(`Paratracheal LNM` ~Age+Sex+`Tumor border`+`Aspect ratio`+`Internal echo pattern` +Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Mulifocality+`Extrathyroidal extension`+`Side of position`+`Prelaryngeal LNM`+`Pretracheal LNM`+`Con-Paratracheal LNM`+`LNM-prRLN`, data = data23, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Paratracheal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data23, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 计算95%置信区间
coef_df$LL <- coef_df$ci_lower
coef_df$UL <- coef_df$ci_upper

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 32.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#79902D"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df$p_value < 0.05,"#79902D","black")))

# 显示初始森林图
print(forest_plot)


coef_df_sorted <- coef_df[order(coef_df$variable), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#79902D"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#79902D", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)


# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.3.Paratracheal-forest-1"))

# 保存CSV文件
write.csv(coef_df[nrow(coef_df):1, ], file = file.path(output_folder, "2.3.Paratracheal-forest-1.csv"), row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = file.path(output_folder, "2.3.Paratracheal-forest-1.csv"), row.names = FALSE)


print(forest_plot_sorted)


```

##4.2.1列线图以及验证曲线
```{r}
# 读取数据
data23 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv")

```

```{r}
data23$Age<-factor(data23$Age,levels = c(0,1),labels = c("Age≤45","Age>45"))
data23$Sex<-factor(data23$Sex,levels = c(0,1),labels = c("Female","Male"))
data23$BMI<-factor(data23$BMI,levels = c(0,1,2),labels = c("Underweight","Normal","Overweight"))

data23$Tumor.border<-factor(data23$Tumor.border,levels = c(0,1,2),labels = c("smooth or borderless","irregular shape or lsharpobed","extrandular invasion"))
data23$Aspect.ratio<-factor(data23$Aspect.ratio,levels = c(0,1),labels = c("≤1",">1"))
 data23$Composition<-factor(data23$Composition,levels = c(0,1,2),labels = c("cystic/cavernous","Mixed cystic and solid","solid"))
 data23$Internal.echo.pattern<-factor(data23$Internal.echo.pattern,levels = c(0,1,2,3),labels = c("echoless","high/isoechoic","hypoechoic","very hypoechoic"))
 data23$Internal.echo.homogeneous<-factor(data23$Internal.echo.homogeneous,levels = c(0,1),labels = c("Non-uniform","Uniform"))
 data23$Calcification<-factor(data23$Calcification,levels = c(0,1,2,3),labels = c("no or large comet tail", "coarse calcification","peripheral calcification","Microcalcification"))
data23$Tumor.internal.vascularization<-factor(data23$Tumor.internal.vascularization,levels = c(0,1),labels = c("Without","Abundant"))
data23$Tumor.Peripheral.blood.flow<-factor(data23$Tumor.Peripheral.blood.flow,levels = c(0,1),labels = c("Without","Abundant"))
data23$Size<-factor(data23$Size,levels = c(0,1),labels = c("≤5", ">5"))
data23$Location<-factor(data23$Location,levels = c(0,1),labels = c("Non-upper","Upper"))
data23$Mulifocality<-factor(data23$Mulifocality,levels = c(1,0),labels = c("Abundant", "Without"))
data23$Hashimoto<-factor(data23$Hashimoto,levels = c(1,0),labels = c("Abundant", "Without"))
data23$Extrathyroidal.extension<-factor(data23$Extrathyroidal.extension,levels = c(1,0),labels = c("Abundant", "Without"))
data23$Side.of.position<-factor(data23$Side.of.position,levels = c(0,1,2,3),labels = c("left","right","bilateral" ,"isthmus"))




data23$Prelaryngeal.LNM<-factor(data23$Prelaryngeal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data23$Pretracheal.LNM<-factor(data23$Pretracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data23$Paratracheal.LNM<-factor(data23$Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data23$Con.Paratracheal.LNM<-factor(data23$Con.Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data23$LNM.prRLN<-factor(data23$LNM.prRLN,levels = c(0,1),labels = c("No", "Yes"))
data23$Total.Central.Lymph.Node.Metastasis<-factor(data23$Total.Central.Lymph.Node.Metastasis,levels = c(0,1),labels = c("No", "Yes"))

```

```{r}
# 加载必要的包
library(rms)

# 准备数据
x <- as.data.frame(data23)
dd <- datadist(data23)
options(datadist = 'dd')

# 拟合逻辑回归模型并指定 x=TRUE 和 y=TRUE
fit1 <- lrm(Paratracheal.LNM ~Sex+Tumor.border+Aspect.ratio+Size+Extrathyroidal.extension+Prelaryngeal.LNM+Pretracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN, data = data23, x = TRUE, y = TRUE)

# 查看模型摘要
summary(fit1)

# 创建列线图
nom1 <- nomogram(fit1, fun = plogis, fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999), lp = FALSE, funlabel = "Paratracheal Lymph Node Metastasis")
plot(nom1)

# 验证曲线
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))

# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.3.Paratracheal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom1)
dev.off()

# 保存验证曲线为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.3.Paratracheal-calibration.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
dev.off()


```

```{r}
# 改变尺寸的列线图
par(mar = c(1, 2, 2, 2))  # 调整绘图边距

# 创建 nomogram
nom2 <- nomogram(fit1, fun = plogis, fun.at = c(0.001, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.999),
                 lp = FALSE, funlabel="Paratracheal Lymph Node Metastasis")

# 绘制 nomogram
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)


# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.3.Paratracheal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)
dev.off()

```
##4.2.2传统预测模型的Roc曲线

```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1P编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2P编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Paratracheal.LNM ~Sex+Tumor.border+Aspect.ratio+Size+Extrathyroidal.extension+Prelaryngeal.LNM+Pretracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Paratracheal.LNM
test_response <- test_data$Paratracheal.LNM
val_response1 <- val_data1$Paratracheal.LNM
val_response2 <- val_data2$Paratracheal.LNM
# 创建ROC对象
train_roc <- roc(train_response, train_probs)
test_roc <- roc(test_response, test_probs)
val_roc1 <- roc(val_response1, val_probs1)
val_roc2 <- roc(val_response2, val_probs2)

# 提取ROC曲线的坐标点
train_roc_data <- coords(train_roc, "all", ret = c("specificity", "sensitivity"))
test_roc_data <- coords(test_roc, "all", ret = c("specificity", "sensitivity"))
val_roc_data1 <- coords(val_roc1, "all", ret = c("specificity", "sensitivity"))
val_roc_data2 <- coords(val_roc2, "all", ret = c("specificity", "sensitivity"))

# 转换为数据框
train_roc_data <- as.data.frame(train_roc_data)
test_roc_data <- as.data.frame(test_roc_data)
val_roc_data1 <- as.data.frame(val_roc_data1)
val_roc_data2 <- as.data.frame(val_roc_data2)

# 绘制ROC曲线
roc_plot <- ggplot() +
  geom_line(data = train_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#3D5714", size = 0.6) +
  geom_line(data = test_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#79902D", size = 0.6) +
  geom_line(data = val_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#5AB682", size = 0.6) +
  geom_line(data = val_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#C8E4D2", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC for Paratracheal Lymph Node Metastasis Nomogram Prediction",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.4, label = paste("Train set AUC =", round(auc(train_roc), 3)), size = 4, color = "#3D5714")  +
  annotate("text", x = 0.7, y = 0.3, label = paste("Test set AUC =", round(auc(test_roc), 3)), size = 4, color = "#79902D")+
  annotate("text", x = 0.7, y = 0.2, label = paste("Validation set1 AUC =", round(auc(val_roc1), 3)), size = 4, color = "#5AB682")+
  annotate("text", x = 0.7, y = 0.1, label = paste("Validation set2 AUC =", round(auc(val_roc2), 3)), size = 4, color = "#C8E4D2")
# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.3.Paratracheal-roc_curve.tiff", plot = roc_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")


```
##4.2.3传统预测模型的dca曲线
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1P编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2P编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Paratracheal.LNM ~Sex+Tumor.border+Aspect.ratio+Size+Extrathyroidal.extension+Prelaryngeal.LNM+Pretracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())


# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Paratracheal.LNM
test_response <- test_data$Paratracheal.LNM
val_response1 <- val_data1$Paratracheal.LNM
val_response2 <- val_data2$Paratracheal.LNM
# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits <- sapply(thresholds, function(x) net_benefit(train_probs, train_response, x))
test_net_benefits <- sapply(thresholds, function(x) net_benefit(test_probs, test_response, x))
val_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs1, val_response1, x))
val_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs2, val_response2, x))


# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(train_response)), train_response, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 找到最大净收益点
train_max_nb <- max(train_net_benefits)
train_max_nb_threshold <- thresholds[which.max(train_net_benefits)]
test_max_nb <- max(test_net_benefits)
test_max_nb_threshold <- thresholds[which.max(test_net_benefits)]
val_max_nb1 <- max(val_net_benefits1)
val_max_nb_threshold1 <- thresholds[which.max(val_net_benefits1)]
val_max_nb2 <- max(val_net_benefits2)
val_max_nb_threshold2 <- thresholds[which.max(val_net_benefits2)]




# 绘制DCA曲线
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit = train_net_benefits,
  test_net_benefit = test_net_benefits,
  val_net_benefit1 = val_net_benefits1,
  val_net_benefit2 = val_net_benefits2,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data, aes(x = threshold)) +
  geom_line(aes(y = train_net_benefit, color = "Train set"), size = 0.6) +
  geom_line(aes(y = test_net_benefit, color = "Test set"), size = 0.6) +
  geom_line(aes(y = val_net_benefit1, color = "Validation set1"), size = 0.6) +
  geom_line(aes(y = val_net_benefit2, color = "Validation set2"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Paratracheal Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(values = c("Train set" = "#3D5714", "Test set" = "#79902D", "Validation set1" = "#5AB682", "Validation set2" = "#C8E4D2","All" = "grey", "None" = "black")) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotate("text", x = 0.3, y = 0.09, label = "Train set", size = 4, color = "#3D5714") +
  annotate("text", x = 0.3, y = 0.12, label = "Test set", size = 4, color = "#79902D") +
  annotate("text", x = 0.3, y = 0.2, label = "Validation set1", size = 4, color = "#5AB682") +
  annotate("text", x = 0.3, y = 0.15, label = "Validation set2", size = 4, color = "#C8E4D2") +
  annotate("text", x = train_max_nb_threshold, y = train_max_nb, label = sprintf("Max: %.3f", train_max_nb), color = "#3D5714", vjust = -1) +
  annotate("text", x = test_max_nb_threshold, y = test_max_nb, label = sprintf("Max: %.3f", test_max_nb), color = "#79902D", vjust = -1) +
   annotate("text", x = val_max_nb_threshold1, y = val_max_nb1, label = sprintf("Max: %.3f", val_max_nb1), color = "#5AB682", vjust = -1) +
   annotate("text", x = val_max_nb_threshold2, y = val_max_nb2, label = sprintf("Max: %.3f", val_max_nb2), color = "#C8E4D2", vjust = -1) +
  coord_cartesian(ylim = c(-0.05, 0.4), xlim = c(0, 0.5))


# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.3.Paratracheal-dca_curve.tiff", plot = dca_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")

print(dca_plot)

```
##4.2.4 保存胜率概率
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1P编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2P编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit4 <- lrm(Paratracheal.LNM ~Sex+Tumor.border+Aspect.ratio+Size+Extrathyroidal.extension+Prelaryngeal.LNM+Pretracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data,  x = TRUE, y = TRUE)

#删掉了一些
nom4 <- predict(fit4, type = "fitted")

# 导出预测结果
nomogram_predictions4 <- data.frame(nomogram_prediction = nom4)
write.csv(nomogram_predictions4, '/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/5.对比/4.Paratracheal.LNM.nomogram_predictions.csv', row.names = FALSE)

```

#5.对侧气管旁的传统预测模型的建立
##5.1.1森林图之1:按or值排序的森林图
```{r}
library(car)
library(ggplot2)


# 读取数据
data24 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv", check.names = FALSE)

# 设置列名
colnames(data24) <- c("Age","Sex","BMI","Tumor border","Aspect ratio","Composition","Internal echo pattern","Internal echo homogeneous","Calcification","Tumor internal vascularization","Tumor Peripheral blood flow","Size","Location","Mulifocality","Hashimoto","Extrathyroidal extension","Side of position","Prelaryngeal LNM","Pretracheal LNM","Paratracheal LNM","Con-Paratracheal LNM","LNM-prRLN","Total Central Lymph Node Metastasis","age","bmi","size","Prelaryngeal LNMR","Prelaryngeal NLNM","Pretracheal LNMR","Pretracheal NLNM","Paratracheal LNMR","Paratracheal NLNM","Con-Paratracheal LNMR","Con-Paratracheal NLNM","LNMR-prRLN","NLNM-prRLN","TCLNMR","TCNLNM")



# 构建初始模型
initial_model <- glm(`Con-Paratracheal LNM` ~Age+Sex+`Tumor border`+Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+`Extrathyroidal extension`+`Side of position`+`Prelaryngeal LNM`+`Pretracheal LNM`+`Paratracheal LNM`+`LNM-prRLN`, data = data24, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Con-Paratracheal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data24, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "black") +
  geom_point(aes(color = p_value < 0.05), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""), 
                x = -20, hjust = -0.1), size = 2.5) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1), size = 2.5) +
  coord_cartesian(xlim = c(-20, 20)) +
  scale_color_manual(values = c("black", "#4E6691"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Con-Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal()

# 保存图像函数
save_forest_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)

}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot(forest_plot, file.path(output_folder, "2.4.Con-Paratracheal-forest-1"))

print(forest_plot)



coef_df_sorted <- coef_df[order(coef_df$odds_ratio), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#4E6691"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Con-Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#4E6691", "black")))


# 显示排序后的森林图
print(forest_plot_sorted)

# 导出结果到CSV文件并反转顺序
write.csv(coef_df[nrow(coef_df):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.4.Con-Paratracheal-forest-2.csv", row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.4.Con-Paratracheal-forest-2.csv", row.names = FALSE)

# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.4.Con-Paratracheal-forest-2"))
print(forest_plot_sorted)
```

##5.1.2森林图之2:按输入顺序排列的森林图

```{r}
library(car)
library(ggplot2)

# 构建初始模型
initial_model <- glm(`Con-Paratracheal LNM` ~Age+Sex+`Tumor border`+Calcification+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+`Extrathyroidal extension`+`Side of position`+`Prelaryngeal LNM`+`Pretracheal LNM`+`Paratracheal LNM`+`LNM-prRLN`, data = data24, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`Con-Paratracheal LNM` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data24, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 计算95%置信区间
coef_df$LL <- coef_df$ci_lower
coef_df$UL <- coef_df$ci_upper

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 32.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#4E6691"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Con-Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df$p_value < 0.05,"#4E6691","black")))

# 显示初始森林图
print(forest_plot)


coef_df_sorted <- coef_df[order(coef_df$variable), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#4E6691"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Con-Paratracheal Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#4E6691", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)


# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.4.Con-Paratracheal-forest-1"))

# 保存CSV文件
write.csv(coef_df[nrow(coef_df):1, ], file = file.path(output_folder, "2.4.Con-Paratracheal-forest-1.csv"), row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = file.path(output_folder, "2.4.Con-Paratracheal-forest-1.csv"), row.names = FALSE)


print(forest_plot_sorted)


```

##5.2.1列线图以及验证曲线
```{r}
# 读取数据
data24 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")

```

```{r}
data24$Age<-factor(data24$Age,levels = c(0,1),labels = c("Age≤45","Age>45"))
data24$Sex<-factor(data24$Sex,levels = c(0,1),labels = c("Female","Male"))
data24$BMI<-factor(data24$BMI,levels = c(0,1,2),labels = c("Underweight","Normal","Overweight"))

data24$Tumor.border<-factor(data24$Tumor.border,levels = c(0,1,2),labels = c("smooth or borderless","irregular shape or lsharpobed","extrandular invasion"))
data24$Aspect.ratio<-factor(data24$Aspect.ratio,levels = c(0,1),labels = c("≤1",">1"))
 data24$Composition<-factor(data24$Composition,levels = c(0,1,2),labels = c("cystic/cavernous","Mixed cystic and solid","solid"))
 data24$Internal.echo.pattern<-factor(data24$Internal.echo.pattern,levels = c(0,1,2,3),labels = c("echoless","high/isoechoic","hypoechoic","very hypoechoic"))
 data24$Internal.echo.homogeneous<-factor(data24$Internal.echo.homogeneous,levels = c(0,1),labels = c("Non-uniform","Uniform"))
 data24$Calcification<-factor(data24$Calcification,levels = c(0,1,2,3),labels = c("no or large comet tail", "coarse calcification","peripheral calcification","Microcalcification"))
data24$Tumor.internal.vascularization<-factor(data24$Tumor.internal.vascularization,levels = c(0,1),labels = c("Without","Abundant"))
data24$Tumor.Peripheral.blood.flow<-factor(data24$Tumor.Peripheral.blood.flow,levels = c(0,1),labels = c("Without","Abundant"))
data24$Size<-factor(data24$Size,levels = c(0,1),labels = c("≤5", ">5"))
data24$Location<-factor(data24$Location,levels = c(0,1),labels = c("Non-upper","Upper"))
data24$Mulifocality<-factor(data24$Mulifocality,levels = c(1,0),labels = c("Abundant", "Without"))
data24$Hashimoto<-factor(data24$Hashimoto,levels = c(1,0),labels = c("Abundant", "Without"))
data24$Extrathyroidal.extension<-factor(data24$Extrathyroidal.extension,levels = c(1,0),labels = c("Abundant", "Without"))
data24$Side.of.position<-factor(data24$Side.of.position,levels = c(0,1,2,3),labels = c("left","right","bilateral" ,"isthmus"))




data24$Prelaryngeal.LNM<-factor(data24$Prelaryngeal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data24$Pretracheal.LNM<-factor(data24$Pretracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data24$Paratracheal.LNM<-factor(data24$Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data24$Con.Paratracheal.LNM<-factor(data24$Con.Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data24$LNM.prRLN<-factor(data24$LNM.prRLN,levels = c(0,1),labels = c("No", "Yes"))
data24$Total.Central.Lymph.Node.Metastasis<-factor(data24$Total.Central.Lymph.Node.Metastasis,levels = c(0,1),labels = c("No", "Yes"))

```

```{r}
# 加载必要的包
library(rms)

# 准备数据
x <- as.data.frame(data24)
dd <- datadist(data24)
options(datadist = 'dd')

# 拟合逻辑回归模型并指定 x=TRUE 和 y=TRUE
fit1 <- lrm(Con.Paratracheal.LNM ~Side.of.position+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN, data = data24, x = TRUE, y = TRUE)

# 查看模型摘要
summary(fit1)

# 创建列线图
nom1 <- nomogram(fit1, fun = plogis, fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999), lp = FALSE, funlabel = "Con-Paratracheal Lymph Node Metastasis")
plot(nom1)

# 验证曲线
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))

# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.4.Con-Paratracheal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom1)
dev.off()

# 保存验证曲线为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.4.Con-Paratracheal-calibration.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
dev.off()


```

```{r}
# 改变尺寸的列线图
par(mar = c(1, 2, 2, 2))  # 调整绘图边距

# 创建 nomogram
nom2 <- nomogram(fit1, fun = plogis, fun.at = c(0.001, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.999),
                 lp = FALSE, funlabel="Con-Paratracheal Lymph Node Metastasis")

# 绘制 nomogram
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)


# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.4.Con-Paratracheal-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)
dev.off()

```
##5.2.2传统预测模型的Roc曲线

```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1CP编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2CP编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Con.Paratracheal.LNM ~Side.of.position+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Con.Paratracheal.LNM
test_response <- test_data$Con.Paratracheal.LNM
val_response1 <- val_data1$Con.Paratracheal.LNM
val_response2 <- val_data2$Con.Paratracheal.LNM
# 创建ROC对象
train_roc <- roc(train_response, train_probs)
test_roc <- roc(test_response, test_probs)
val_roc1 <- roc(val_response1, val_probs1)
val_roc2 <- roc(val_response2, val_probs2)

# 提取ROC曲线的坐标点
train_roc_data <- coords(train_roc, "all", ret = c("specificity", "sensitivity"))
test_roc_data <- coords(test_roc, "all", ret = c("specificity", "sensitivity"))
val_roc_data1 <- coords(val_roc1, "all", ret = c("specificity", "sensitivity"))
val_roc_data2 <- coords(val_roc2, "all", ret = c("specificity", "sensitivity"))

# 转换为数据框
train_roc_data <- as.data.frame(train_roc_data)
test_roc_data <- as.data.frame(test_roc_data)
val_roc_data1 <- as.data.frame(val_roc_data1)
val_roc_data2 <- as.data.frame(val_roc_data2)

# 绘制ROC曲线
roc_plot <- ggplot() +
  geom_line(data = train_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#82A7D1", size = 0.6) +
  geom_line(data = test_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#4E6691", size = 0.6) +
  geom_line(data = val_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#B6D7E9", size = 0.6) +
  geom_line(data = val_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#DBEAF3", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC for Con-Paratracheal Lymph Node Metastasis Nomogram Prediction",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.4, label = paste("Train set AUC =", round(auc(train_roc), 3)), size = 4, color = "#82A7D1")  +
  annotate("text", x = 0.7, y = 0.3, label = paste("Test set AUC =", round(auc(test_roc), 3)), size = 4, color = "#4E6691")+
  annotate("text", x = 0.7, y = 0.2, label = paste("Validation set1 AUC =", round(auc(val_roc1), 3)), size = 4, color = "#B6D7E9")+
  annotate("text", x = 0.7, y = 0.1, label = paste("Validation set2 AUC =", round(auc(val_roc2), 3)), size = 4, color = "#DBEAF3")
# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.4.Con-Paratracheal-roc_curve.tiff", plot = roc_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")


```
##5.2.3传统预测模型的dca曲线
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1CP编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2CP编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(Con.Paratracheal.LNM ~Side.of.position+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data, family = binomial())


# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$Con.Paratracheal.LNM
test_response <- test_data$Con.Paratracheal.LNM
val_response1 <- val_data1$Con.Paratracheal.LNM
val_response2 <- val_data2$Con.Paratracheal.LNM
# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits <- sapply(thresholds, function(x) net_benefit(train_probs, train_response, x))
test_net_benefits <- sapply(thresholds, function(x) net_benefit(test_probs, test_response, x))
val_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs1, val_response1, x))
val_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs2, val_response2, x))


# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(train_response)), train_response, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 找到最大净收益点
train_max_nb <- max(train_net_benefits)
train_max_nb_threshold <- thresholds[which.max(train_net_benefits)]
test_max_nb <- max(test_net_benefits)
test_max_nb_threshold <- thresholds[which.max(test_net_benefits)]
val_max_nb1 <- max(val_net_benefits1)
val_max_nb_threshold1 <- thresholds[which.max(val_net_benefits1)]
val_max_nb2 <- max(val_net_benefits2)
val_max_nb_threshold2 <- thresholds[which.max(val_net_benefits2)]




# 绘制DCA曲线
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit = train_net_benefits,
  test_net_benefit = test_net_benefits,
  val_net_benefit1 = val_net_benefits1,
  val_net_benefit2 = val_net_benefits2,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data, aes(x = threshold)) +
  geom_line(aes(y = train_net_benefit, color = "Train set"), size = 0.6) +
  geom_line(aes(y = test_net_benefit, color = "Test set"), size = 0.6) +
  geom_line(aes(y = val_net_benefit1, color = "Validation set1"), size = 0.6) +
  geom_line(aes(y = val_net_benefit2, color = "Validation set2"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Con-Paratracheal Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(values = c("Train set" = "#82A7D1", "Test set" = "#4E6691", "Validation set1" = "#B6D7E9", "Validation set2" = "#DBEAF3","All" = "grey", "None" = "black")) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotate("text", x = 0.2, y = 0.02, label = "Train set", size = 4, color = "#82A7D1") +
  annotate("text", x = 0.2, y = 0.05, label = "Test set", size = 4, color = "#4E6691") +
  annotate("text", x = 0.2, y = 0.08, label = "Validation set1", size = 4, color = "#B6D7E9") +
  annotate("text", x = 0.2, y = 0.13, label = "Validation set2", size = 4, color = "#DBEAF3") +
  annotate("text", x = train_max_nb_threshold, y = train_max_nb, label = sprintf("Max: %.3f", train_max_nb), color = "#82A7D1", vjust = -1) +
  annotate("text", x = test_max_nb_threshold, y = test_max_nb, label = sprintf("Max: %.3f", test_max_nb), color = "#4E6691", vjust = -1) +
   annotate("text", x = val_max_nb_threshold1, y = val_max_nb1, label = sprintf("Max: %.3f", val_max_nb1), color = "#B6D7E9", vjust = -1) +
   annotate("text", x = val_max_nb_threshold2, y = val_max_nb2, label = sprintf("Max: %.3f", val_max_nb2), color = "#DBEAF3", vjust = -1) +
  coord_cartesian(ylim = c(-0.05, 0.3), xlim = c(0, 0.5))


# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.4.Con-Paratracheal-dca_curve.tiff", plot = dca_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")

print(dca_plot)

```
##5.2.4 保存胜率概率
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1CP编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2CP编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")


# 构建模型
fit5 <- lrm(Con.Paratracheal.LNM ~Side.of.position+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data,  x = TRUE, y = TRUE)

#删掉了一些
nom5 <- predict(fit5, type = "fitted")

# 导出预测结果
nomogram_predictions5 <- data.frame(nomogram_prediction = nom5)
write.csv(nomogram_predictions5, '/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/5.对比/5.Con.Paratracheal.LNM.nomogram_predictions.csv', row.names = FALSE)

```

#6.喉返后的传统预测模型的建立
##6.1.1森林图之1:按or值排序的森林图
```{r}
library(car)
library(ggplot2)


# 读取数据
data25 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv", check.names = FALSE)

# 设置列名
colnames(data25) <- c("Age","Sex","BMI","Tumor border","Aspect ratio","Composition","Internal echo pattern","Internal echo homogeneous","Calcification","Tumor internal vascularization","Tumor Peripheral blood flow","Size","Location","Mulifocality","Hashimoto","Extrathyroidal extension","Side of position","Prelaryngeal LNM","Pretracheal LNM","Paratracheal LNM","Con-Paratracheal LNM","LNM-prRLN","Total Central Lymph Node Metastasis","age","bmi","size","Prelaryngeal LNMR","Prelaryngeal NLNM","Pretracheal LNMR","Pretracheal NLNM","Paratracheal LNMR","Paratracheal NLNM","Con-Paratracheal LNMR","Con-Paratracheal NLNM","LNMR-prRLN","NLNM-prRLN","TCLNMR","TCNLNM")



# 构建初始模型
initial_model <- glm(`LNM-prRLN` ~Age+Sex+`Tumor border`+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Mulifocality+`Extrathyroidal extension`+`Prelaryngeal LNM`+`Pretracheal LNM`+`Paratracheal LNM`+`Con-Paratracheal LNM`, data = data25, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 5])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`LNM-prRLN` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data25, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "black") +
  geom_point(aes(color = p_value < 0.05), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""), 
                x = -20, hjust = -0.1), size = 2.5) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1), size = 2.5) +
  coord_cartesian(xlim = c(-20, 20)) +
  scale_color_manual(values = c("black", "#D355FF"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal()

# 保存图像函数
save_forest_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)

}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot(forest_plot, file.path(output_folder, "2.5.LNM-prRLN-forest-1"))

print(forest_plot)



coef_df_sorted <- coef_df[order(coef_df$odds_ratio), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#D355FF"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#D355FF", "black")))


# 显示排序后的森林图
print(forest_plot_sorted)

# 导出结果到CSV文件并反转顺序
write.csv(coef_df[nrow(coef_df):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.5.LNM-prRLN-forest-2.csv", row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图/2.5.LNM-prRLN-forest-2.csv", row.names = FALSE)

# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.5.LNM-prRLN-forest-2"))
print(forest_plot_sorted)
```

##6.1.2森林图之2:按输入顺序排列的森林图

```{r}
library(car)
library(ggplot2)

# 构建初始模型
initial_model <- glm(`LNM-prRLN` ~Age+Sex+`Tumor border`+`Tumor internal vascularization`+`Tumor Peripheral blood flow`+Size+Mulifocality+`Extrathyroidal extension`+`Prelaryngeal LNM`+`Pretracheal LNM`+`Paratracheal LNM`+`Con-Paratracheal LNM`, data = data25, family = binomial())


# 计算VIF值
vif_values <- vif(initial_model)
print(vif_values)

# 移除高VIF值的变量（假设阈值为5）
selected_vars <- names(vif_values[vif_values < 10])

# 重新构建模型，消除共线性
formula <- as.formula(paste("`LNM-prRLN` ~", paste(selected_vars, collapse = " + ")))
final_model <- glm(formula, data = data25, family = binomial())

# 提取模型系数
coefficients <- coef(final_model)

# 创建系数数据框
coef_df <- data.frame(
  variable = names(coefficients),
  coefficient = coefficients,
  odds_ratio = exp(coefficients),
  p_value = summary(final_model)$coefficients[, "Pr(>|z|)"],
  ci_lower = exp(confint(final_model)[, 1]),
  ci_upper = exp(confint(final_model)[, 2])
)

# 计算95%置信区间
coef_df$LL <- coef_df$ci_lower
coef_df$UL <- coef_df$ci_upper

# 将(Intercept)标签改为Intercept
coef_df$variable[coef_df$variable == "(Intercept)"] <- "Intercept"

# 手动设置变量顺序并反转
variable_order <- c("Intercept", selected_vars)
coef_df$variable <- factor(coef_df$variable, levels = rev(variable_order))

# 创建初始森林图
forest_plot <- ggplot(coef_df, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 32.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#D355FF"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df$p_value < 0.05,"#D355FF","black")))

# 显示初始森林图
print(forest_plot)


coef_df_sorted <- coef_df[order(coef_df$variable), ]
coef_df_sorted <- rbind(coef_df_sorted[coef_df_sorted$variable != "Intercept", ], coef_df_sorted[coef_df_sorted$variable == "Intercept", ])
coef_df_sorted$variable <- factor(coef_df_sorted$variable, levels = coef_df_sorted$variable)

forest_plot_sorted <- ggplot(coef_df_sorted, aes(x = odds_ratio, y = variable)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = p_value < 0.05), height = 0.2) +
  geom_point(aes(color = p_value < 0.05), size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste(round(odds_ratio, 3), " (", round(ci_lower, 3), " - ", round(ci_upper, 3), ")", sep = ""),
                x = -5, hjust = -0.1, color = p_value < 0.05), size = 3.2) +
  geom_text(aes(label = paste("p =", round(p_value, 3)), x = 9, hjust = 1.1, color = p_value < 0.05), size = 3.2) +
  coord_cartesian(xlim = c(-5, 9)) +
  scale_color_manual(values = c("black", "#D355FF"), labels = c("p >= 0.05", "p < 0.05")) +
  labs(x = "Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Odds Ratio", y = "Variable") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0.5, color = ifelse(coef_df_sorted$p_value < 0.05, "#D355FF", "black")))

# 显示排序后的森林图
print(forest_plot_sorted)


# 保存图像函数
save_forest_plot_sorted <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/1.森林图"

# 保存
save_forest_plot_sorted(forest_plot_sorted, file.path(output_folder, "2.5.LNM-prRLN-forest-1"))

# 保存CSV文件
write.csv(coef_df[nrow(coef_df):1, ], file = file.path(output_folder, "2.5.LNM-prRLN-forest-1.csv"), row.names = FALSE)
write.csv(coef_df_sorted[nrow(coef_df_sorted):1, ], file = file.path(output_folder, "2.5.LNM-prRLN-forest-1.csv"), row.names = FALSE)


print(forest_plot_sorted)


```

##6.2.1列线图以及验证曲线
```{r}
# 读取数据
data25 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv")

```

```{r}
data25$Age<-factor(data25$Age,levels = c(0,1),labels = c("Age≤45","Age>45"))
data25$Sex<-factor(data25$Sex,levels = c(0,1),labels = c("Female","Male"))
data25$BMI<-factor(data25$BMI,levels = c(0,1,2),labels = c("Underweight","Normal","Overweight"))

data25$Tumor.border<-factor(data25$Tumor.border,levels = c(0,1,2),labels = c("smooth or borderless","irregular shape or lsharpobed","extrandular invasion"))
data25$Aspect.ratio<-factor(data25$Aspect.ratio,levels = c(0,1),labels = c("≤1",">1"))
 data25$Composition<-factor(data25$Composition,levels = c(0,1,2),labels = c("cystic/cavernous","Mixed cystic and solid","solid"))
 data25$Internal.echo.pattern<-factor(data25$Internal.echo.pattern,levels = c(0,1,2,3),labels = c("echoless","high/isoechoic","hypoechoic","very hypoechoic"))
 data25$Internal.echo.homogeneous<-factor(data25$Internal.echo.homogeneous,levels = c(0,1),labels = c("Non-uniform","Uniform"))
 data25$Calcification<-factor(data25$Calcification,levels = c(0,1,2,3),labels = c("no or large comet tail", "coarse calcification","peripheral calcification","Microcalcification"))
data25$Tumor.internal.vascularization<-factor(data25$Tumor.internal.vascularization,levels = c(0,1),labels = c("Without","Abundant"))
data25$Tumor.Peripheral.blood.flow<-factor(data25$Tumor.Peripheral.blood.flow,levels = c(0,1),labels = c("Without","Abundant"))
data25$Size<-factor(data25$Size,levels = c(0,1),labels = c("≤5", ">5"))
data25$Location<-factor(data25$Location,levels = c(0,1),labels = c("Non-upper","Upper"))
data25$Mulifocality<-factor(data25$Mulifocality,levels = c(1,0),labels = c("Abundant", "Without"))
data25$Hashimoto<-factor(data25$Hashimoto,levels = c(1,0),labels = c("Abundant", "Without"))
data25$Extrathyroidal.extension<-factor(data25$Extrathyroidal.extension,levels = c(1,0),labels = c("Abundant", "Without"))
data25$Side.of.position<-factor(data25$Side.of.position,levels = c(0,1,2,3),labels = c("left","right","bilateral" ,"isthmus"))




data25$Prelaryngeal.LNM<-factor(data25$Prelaryngeal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data25$Pretracheal.LNM<-factor(data25$Pretracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data25$Paratracheal.LNM<-factor(data25$Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data25$Con.Paratracheal.LNM<-factor(data25$Con.Paratracheal.LNM,levels = c(0,1),labels = c("No", "Yes"))
data25$LNM.prRLN<-factor(data25$LNM.prRLN,levels = c(0,1),labels = c("No", "Yes"))
data25$Total.Central.Lymph.Node.Metastasis<-factor(data25$Total.Central.Lymph.Node.Metastasis,levels = c(0,1),labels = c("No", "Yes"))

```

```{r}
# 加载必要的包
library(rms)

# 准备数据
x <- as.data.frame(data25)
dd <- datadist(data25)
options(datadist = 'dd')

# 拟合逻辑回归模型并指定 x=TRUE 和 y=TRUE
fit1 <- lrm(LNM.prRLN ~Prelaryngeal.LNM+Pretracheal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM, data = data25, x = TRUE, y = TRUE)

# 查看模型摘要
summary(fit1)

# 创建列线图
nom1 <- nomogram(fit1, fun = plogis, fun.at = c(.001, .01, .05, seq(.1, .9, by = .1), .95, .99, .999), lp = FALSE, funlabel = "Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis")
plot(nom1)

# 验证曲线
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))

# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.5.LNM-prRLN-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom1)
dev.off()

# 保存验证曲线为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.5.LNM-prRLN-calibration.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(cal1, xlim = c(0, 1.0), ylim = c(0, 1.0))
dev.off()


```

```{r}
# 改变尺寸的列线图
par(mar = c(1, 2, 2, 2))  # 调整绘图边距

# 创建 nomogram
nom2 <- nomogram(fit1, fun = plogis, fun.at = c(0.001, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.999),
                 lp = FALSE, funlabel="Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis")

# 绘制 nomogram
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)


# 保存列线图为1200 DPI的 .tiff格式
tiff("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/2.列线图和验证曲线/2.5.LNM-prRLN-nomogram.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
plot(nom2, abbreviate = FALSE, col.lines = "blue", col.points = "blue", cex.names = 0.12, cex.axis = 0.52,#这是列线图的线的字的大小
     cex.lab = 30, lwd.lines = 30, lwd.funnel = 30, cex.var = 0.6, varname.dist = 2000)
dev.off()

```
##6.2.2传统预测模型的Roc曲线

```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1F编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2F编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(LNM.prRLN ~Prelaryngeal.LNM+Pretracheal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM,
            data = tra_data, family = binomial())

# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$LNM.prRLN
test_response <- test_data$LNM.prRLN
val_response1 <- val_data1$LNM.prRLN
val_response2 <- val_data2$LNM.prRLN
# 创建ROC对象
train_roc <- roc(train_response, train_probs)
test_roc <- roc(test_response, test_probs)
val_roc1 <- roc(val_response1, val_probs1)
val_roc2 <- roc(val_response2, val_probs2)

# 提取ROC曲线的坐标点
train_roc_data <- coords(train_roc, "all", ret = c("specificity", "sensitivity"))
test_roc_data <- coords(test_roc, "all", ret = c("specificity", "sensitivity"))
val_roc_data1 <- coords(val_roc1, "all", ret = c("specificity", "sensitivity"))
val_roc_data2 <- coords(val_roc2, "all", ret = c("specificity", "sensitivity"))

# 转换为数据框
train_roc_data <- as.data.frame(train_roc_data)
test_roc_data <- as.data.frame(test_roc_data)
val_roc_data1 <- as.data.frame(val_roc_data1)
val_roc_data2 <- as.data.frame(val_roc_data2)

# 绘制ROC曲线
roc_plot <- ggplot() +
  geom_line(data = train_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#8833D5", size = 0.6) +
  geom_line(data = test_roc_data, aes(x = 1 - specificity, y = sensitivity), color = "#D355FF", size = 0.6) +
  geom_line(data = val_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#E8A4FF", size = 0.6) +
  geom_line(data = val_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#F0CCFF", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC for Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Nomogram Prediction",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.4, label = paste("Train set AUC =", round(auc(train_roc), 3)), size = 4, color = "#8833D5")  +
  annotate("text", x = 0.7, y = 0.3, label = paste("Test set AUC =", round(auc(test_roc), 3)), size = 4, color = "#D355FF")+
  annotate("text", x = 0.7, y = 0.2, label = paste("Validation set1 AUC =", round(auc(val_roc1), 3)), size = 4, color = "#E8A4FF")+
  annotate("text", x = 0.7, y = 0.1, label = paste("Validation set2 AUC =", round(auc(val_roc2), 3)), size = 4, color = "#F0CCFF")
# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.5.LNM-prRLN-roc_curve.tiff", plot = roc_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")


```
##6.2.3传统预测模型的dca曲线
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1CP编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2CP编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit1 <- glm(LNM.prRLN ~Prelaryngeal.LNM+Pretracheal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM,
            data = tra_data, family = binomial())


# 预测概率
train_probs <- predict(fit1, newdata = tra_data, type = "response")
test_probs <- predict(fit1, newdata = test_data, type = "response")
val_probs1 <- predict(fit1, newdata = val_data1, type = "response")
val_probs2 <- predict(fit1, newdata = val_data2, type = "response")


train_response <- tra_data$LNM.prRLN
test_response <- test_data$LNM.prRLN
val_response1 <- val_data1$LNM.prRLN
val_response2 <- val_data2$LNM.prRLN
# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits <- sapply(thresholds, function(x) net_benefit(train_probs, train_response, x))
test_net_benefits <- sapply(thresholds, function(x) net_benefit(test_probs, test_response, x))
val_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs1, val_response1, x))
val_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs2, val_response2, x))


# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(train_response)), train_response, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 找到最大净收益点
train_max_nb <- max(train_net_benefits)
train_max_nb_threshold <- thresholds[which.max(train_net_benefits)]
test_max_nb <- max(test_net_benefits)
test_max_nb_threshold <- thresholds[which.max(test_net_benefits)]
val_max_nb1 <- max(val_net_benefits1)
val_max_nb_threshold1 <- thresholds[which.max(val_net_benefits1)]
val_max_nb2 <- max(val_net_benefits2)
val_max_nb_threshold2 <- thresholds[which.max(val_net_benefits2)]




# 绘制DCA曲线
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit = train_net_benefits,
  test_net_benefit = test_net_benefits,
  val_net_benefit1 = val_net_benefits1,
  val_net_benefit2 = val_net_benefits2,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data, aes(x = threshold)) +
  geom_line(aes(y = train_net_benefit, color = "Train set"), size = 0.6) +
  geom_line(aes(y = test_net_benefit, color = "Test set"), size = 0.6) +
  geom_line(aes(y = val_net_benefit1, color = "Validation set1"), size = 0.6) +
  geom_line(aes(y = val_net_benefit2, color = "Validation set2"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(values = c("Train set" = "#8833D5", "Test set" = "#D355FF", "Validation set1" = "#E8A4FF", "Validation set2" = "#F0CCFF","All" = "grey", "None" = "black")) +
  theme_minimal() +
  theme(legend.position = "right") +
  annotate("text", x = 0.2, y = 0.02, label = "Train set", size = 4, color = "#8833D5") +
  annotate("text", x = 0.2, y = 0.05, label = "Test set", size = 4, color = "#D355FF") +
  annotate("text", x = 0.2, y = 0.08, label = "Validation set1", size = 4, color = "#E8A4FF") +
  annotate("text", x = 0.2, y = 0.13, label = "Validation set2", size = 4, color = "#F0CCFF") +
  annotate("text", x = train_max_nb_threshold, y = train_max_nb, label = sprintf("Max: %.3f", train_max_nb), color = "#8833D5", vjust = -1) +
  annotate("text", x = test_max_nb_threshold, y = test_max_nb, label = sprintf("Max: %.3f", test_max_nb), color = "#D355FF", vjust = -1) +
   annotate("text", x = val_max_nb_threshold1, y = val_max_nb1, label = sprintf("Max: %.3f", val_max_nb1), color = "#E8A4FF", vjust = -1) +
   annotate("text", x = val_max_nb_threshold2, y = val_max_nb2, label = sprintf("Max: %.3f", val_max_nb2), color = "#F0CCFF", vjust = -1) +
  coord_cartesian(ylim = c(-0.05, 0.2), xlim = c(0, 0.33))


# 保存ROC曲线为.tiff格式
ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/3.评价1/3.5.LNM-prRLN-dca_curve.tiff", plot = dca_plot, width = 8, height = 6, units = "in", dpi = 300, compression = "lzw")

print(dca_plot)

```
##5.2.4 保存胜率概率
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv")
val_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1F编码后_插补.csv")
val_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2F编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size <- round(tra_ratio * nrow(train_data))
tra_data <- train_data[shuffled_index[1:tra_size], ]
test_data <- train_data[shuffled_index[(tra_size + 1):nrow(train_data)], ]

cat("训练集观测数量:", nrow(tra_data), "\n")
cat("测试集观测数量:", nrow(test_data), "\n")
cat("验证集1观测数量:", nrow(val_data1), "\n")
cat("验证集2观测数量:", nrow(val_data2), "\n")

# 构建模型
fit6 <- lrm(LNM.prRLN ~Prelaryngeal.LNM+Pretracheal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM,
            data = tra_data,  x = TRUE, y = TRUE)

#删掉了一些
nom6 <- predict(fit6, type = "fitted")

# 导出预测结果
nomogram_predictions <- data.frame(nomogram_prediction = nom6)
write.csv(nomogram_predictions, '/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/5.对比/6.LNM.prRLN.nomogram_predictions.csv', row.names = FALSE)

```

#7.总的rOC曲线
##7.1加载数据和模型拟合
###7.1.1总中央区
```{r}
# 总中央区
library(pROC)
library(ggplot2)

# 读取数据
train_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总编码后_插补.csv")
val_data11 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1编码后_插补.csv")
val_data12 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data1)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size1 <- round(tra_ratio * nrow(train_data1))
tra_data1 <- train_data1[shuffled_index[1:tra_size1], ]
test_data1 <- train_data1[shuffled_index[(tra_size1 + 1):nrow(train_data1)], ]

cat("训练集观测数量:", nrow(tra_data1), "\n")
cat("测试集观测数量:", nrow(test_data1), "\n")
cat("验证集1观测数量:", nrow(val_data11), "\n")
cat("验证集2观测数量:", nrow(val_data12), "\n")

# 构建模型
fit1 <- glm(Total.Central.Lymph.Node.Metastasis ~ Age + Sex + Tumor.border + Aspect.ratio + Calcification + Tumor.Peripheral.blood.flow + Size + Mulifocality + Extrathyroidal.extension,
            data = tra_data1, family = binomial())

# 预测概率
train_probs1 <- predict(fit1, newdata = tra_data1, type = "response")
test_probs1 <- predict(fit1, newdata = test_data1, type = "response")
val_probs11 <- predict(fit1, newdata = val_data11, type = "response")
val_probs12 <- predict(fit1, newdata = val_data12, type = "response")


train_response1 <- tra_data1$Total.Central.Lymph.Node.Metastasis
test_response1 <- test_data1$Total.Central.Lymph.Node.Metastasis
val_response11 <- val_data11$Total.Central.Lymph.Node.Metastasis
val_response12 <- val_data12$Total.Central.Lymph.Node.Metastasis
# 创建ROC对象
train_roc1 <- roc(train_response1, train_probs1)
test_roc1 <- roc(test_response1, test_probs1)
val_roc11 <- roc(val_response11, val_probs11)
val_roc12 <- roc(val_response12, val_probs12)

# 提取ROC曲线的坐标点
train_roc_data1 <- coords(train_roc1, "all")
test_roc_data1 <- coords(test_roc1, "all")
val_roc_data11 <- coords(val_roc11, "all")
val_roc_data12 <- coords(val_roc12, "all")

# 转换为数据框
train_roc_data1 <- as.data.frame(train_roc_data1)
test_roc_data1 <- as.data.frame(test_roc_data1)
val_roc_data11 <- as.data.frame(val_roc_data11)
val_roc_data12 <- as.data.frame(val_roc_data12)
```


###7.1.2喉前
```{r}
# 喉前

# 读取数据
train_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补.csv")
val_data21 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1H编码后_插补.csv")
val_data22 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2H编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data2)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size2 <- round(tra_ratio * nrow(train_data2))
tra_data2 <- train_data2[shuffled_index[1:tra_size2], ]
test_data2 <- train_data2[shuffled_index[(tra_size2 + 1):nrow(train_data2)], ]

cat("训练集观测数量:", nrow(tra_data2), "\n")
cat("测试集观测数量:", nrow(test_data2), "\n")
cat("验证集1观测数量:", nrow(val_data21), "\n")
cat("验证集2观测数量:", nrow(val_data22), "\n")

# 构建模型
fit2 <- glm(Prelaryngeal.LNM ~Location+Hashimoto+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data2, family = binomial())

# 预测概率
train_probs2 <- predict(fit2, newdata = tra_data2, type = "response")
test_probs2 <- predict(fit2, newdata = test_data2, type = "response")
val_probs21 <- predict(fit2, newdata = val_data21, type = "response")
val_probs22 <- predict(fit2, newdata = val_data22, type = "response")


train_response2 <- tra_data2$Prelaryngeal.LNM
test_response2 <- test_data2$Prelaryngeal.LNM
val_response21 <- val_data21$Prelaryngeal.LNM
val_response22 <- val_data22$Prelaryngeal.LNM
# 创建ROC对象
train_roc2 <- roc(train_response2, train_probs2)
test_roc2 <- roc(test_response2, test_probs2)
val_roc21 <- roc(val_response21, val_probs21)
val_roc22 <- roc(val_response22, val_probs22)

# 提取ROC曲线的坐标点
train_roc_data2 <- coords(train_roc2, "all")
test_roc_data2 <- coords(test_roc2, "all")
val_roc_data21 <- coords(val_roc21, "all")
val_roc_data22 <- coords(val_roc22, "all")

# 转换为数据框
train_roc_data2 <- as.data.frame(train_roc_data2)
test_roc_data2 <- as.data.frame(test_roc_data2)
val_roc_data21 <- as.data.frame(val_roc_data21)
val_roc_data22 <- as.data.frame(val_roc_data22)
```
###7.1.3气管前
```{r}
# 气管前

# 读取数据
train_data3 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv")
val_data31 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1Q编码后_插补.csv")
val_data32 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2Q编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data3)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size3 <- round(tra_ratio * nrow(train_data3))
tra_data3 <- train_data3[shuffled_index[1:tra_size3], ]
test_data3 <- train_data3[shuffled_index[(tra_size3 + 1):nrow(train_data3)], ]

cat("训练集观测数量:", nrow(tra_data3), "\n")
cat("测试集观测数量:", nrow(test_data3), "\n")
cat("验证集1观测数量:", nrow(val_data31), "\n")
cat("验证集2观测数量:", nrow(val_data32), "\n")

# 构建模型
fit3 <- glm(Pretracheal.LNM ~Age+Sex+Tumor.Peripheral.blood.flow+Mulifocality+Prelaryngeal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data3, family = binomial())

# 预测概率
train_probs3 <- predict(fit3, newdata = tra_data3, type = "response")
test_probs3 <- predict(fit3, newdata = test_data3, type = "response")
val_probs31 <- predict(fit3, newdata = val_data31, type = "response")
val_probs32 <- predict(fit3, newdata = val_data32, type = "response")


train_response3 <- tra_data3$Pretracheal.LNM
test_response3 <- test_data3$Pretracheal.LNM
val_response31 <- val_data31$Pretracheal.LNM
val_response32 <- val_data32$Pretracheal.LNM
# 创建ROC对象
train_roc3 <- roc(train_response3, train_probs3)
test_roc3 <- roc(test_response3, test_probs3)
val_roc31 <- roc(val_response31, val_probs31)
val_roc32 <- roc(val_response32, val_probs32)

# 提取ROC曲线的坐标点
train_roc_data3 <- coords(train_roc3, "all")
test_roc_data3 <- coords(test_roc3, "all")
val_roc_data31 <- coords(val_roc31, "all")
val_roc_data32 <- coords(val_roc32, "all")

# 转换为数据框
train_roc_data3 <- as.data.frame(train_roc_data3)
test_roc_data3 <- as.data.frame(test_roc_data3)
val_roc_data31 <- as.data.frame(val_roc_data31)
val_roc_data32 <- as.data.frame(val_roc_data32)
```
###7.1.4同侧气管旁
```{r}
# 同侧气管旁

# 读取数据
train_data4 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv")
val_data41 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1P编码后_插补.csv")
val_data42 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2P编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data4)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size4 <- round(tra_ratio * nrow(train_data4))
tra_data4 <- train_data4[shuffled_index[1:tra_size4], ]
test_data4 <- train_data4[shuffled_index[(tra_size4 + 1):nrow(train_data4)], ]

cat("训练集观测数量:", nrow(tra_data4), "\n")
cat("测试集观测数量:", nrow(test_data4), "\n")
cat("验证集1观测数量:", nrow(val_data41), "\n")
cat("验证集2观测数量:", nrow(val_data42), "\n")

# 构建模型
fit4 <- glm(Paratracheal.LNM ~Sex+Tumor.border+Aspect.ratio+Size+Extrathyroidal.extension+Prelaryngeal.LNM+Pretracheal.LNM+Con.Paratracheal.LNM+LNM.prRLN,
            data = tra_data4, family = binomial())

# 预测概率
train_probs4 <- predict(fit4, newdata = tra_data4, type = "response")
test_probs4 <- predict(fit4, newdata = test_data4, type = "response")
val_probs41 <- predict(fit4, newdata = val_data41, type = "response")
val_probs42 <- predict(fit4, newdata = val_data42, type = "response")


train_response4 <- tra_data4$Paratracheal.LNM
test_response4 <- test_data4$Paratracheal.LNM
val_response41 <- val_data41$Paratracheal.LNM
val_response42 <- val_data42$Paratracheal.LNM
# 创建ROC对象
train_roc4 <- roc(train_response4, train_probs4)
test_roc4 <- roc(test_response4, test_probs4)
val_roc41 <- roc(val_response41, val_probs41)
val_roc42 <- roc(val_response42, val_probs42)

# 提取ROC曲线的坐标点
train_roc_data4 <- coords(train_roc4, "all")
test_roc_data4 <- coords(test_roc4, "all")
val_roc_data41 <- coords(val_roc41, "all")
val_roc_data42 <- coords(val_roc42, "all")

# 转换为数据框
train_roc_data4 <- as.data.frame(train_roc_data4)
test_roc_data4 <- as.data.frame(test_roc_data4)
val_roc_data41 <- as.data.frame(val_roc_data41)
val_roc_data42 <- as.data.frame(val_roc_data42)
```
###7.1.5对侧气管旁
```{r}
# 对侧气管旁
# 读取数据
train_data5 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")
val_data51 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1CP编码后_插补.csv")
val_data52 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2CP编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data5)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size5 <- round(tra_ratio * nrow(train_data5))
tra_data5 <- train_data5[shuffled_index[1:tra_size5], ]
test_data5 <- train_data5[shuffled_index[(tra_size5 + 1):nrow(train_data5)], ]

cat("训练集观测数量:", nrow(tra_data5), "\n")
cat("测试集观测数量:", nrow(test_data5), "\n")
cat("验证集1观测数量:", nrow(val_data51), "\n")
cat("验证集2观测数量:", nrow(val_data52), "\n")

# 构建模型
fit5 <- glm(Con.Paratracheal.LNM ~Side.of.position+Pretracheal.LNM+Paratracheal.LNM+LNM.prRLN,
            data = tra_data5, family = binomial())

# 预测概率
train_probs5 <- predict(fit5, newdata = tra_data5, type = "response")
test_probs5 <- predict(fit5, newdata = test_data5, type = "response")
val_probs51 <- predict(fit5, newdata = val_data51, type = "response")
val_probs52 <- predict(fit5, newdata = val_data52, type = "response")


train_response5 <- tra_data5$Con.Paratracheal.LNM
test_response5 <- test_data5$Con.Paratracheal.LNM
val_response51 <- val_data51$Con.Paratracheal.LNM
val_response52 <- val_data52$Con.Paratracheal.LNM
# 创建ROC对象
train_roc5 <- roc(train_response5, train_probs5)
test_roc5 <- roc(test_response5, test_probs5)
val_roc51 <- roc(val_response51, val_probs51)
val_roc52 <- roc(val_response52, val_probs52)

# 提取ROC曲线的坐标点
train_roc_data5 <- coords(train_roc5, "all")
test_roc_data5 <- coords(test_roc5, "all")
val_roc_data51 <- coords(val_roc51, "all")
val_roc_data52 <- coords(val_roc52, "all")

# 转换为数据框
train_roc_data5 <- as.data.frame(train_roc_data5)
test_roc_data5 <- as.data.frame(test_roc_data5)
val_roc_data51 <- as.data.frame(val_roc_data51)
val_roc_data52 <- as.data.frame(val_roc_data52)

```
###7.1.6后饭后
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)

# 读取数据
train_data6 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv")
val_data61 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1F编码后_插补.csv")
val_data62 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2F编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data6)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size6 <- round(tra_ratio * nrow(train_data6))
tra_data6 <- train_data6[shuffled_index[1:tra_size6], ]
test_data6 <- train_data6[shuffled_index[(tra_size6 + 1):nrow(train_data6)], ]

cat("训练集观测数量:", nrow(tra_data6), "\n")
cat("测试集观测数量:", nrow(test_data6), "\n")
cat("验证集1观测数量:", nrow(val_data61), "\n")
cat("验证集2观测数量:", nrow(val_data62), "\n")

# 构建模型
fit6 <- glm(LNM.prRLN ~Prelaryngeal.LNM+Pretracheal.LNM+Paratracheal.LNM+Con.Paratracheal.LNM,
            data = tra_data6, family = binomial())

# 预测概率
train_probs6 <- predict(fit6, newdata = tra_data6, type = "response")
test_probs6 <- predict(fit6, newdata = test_data6, type = "response")
val_probs61 <- predict(fit6, newdata = val_data61, type = "response")
val_probs62 <- predict(fit6, newdata = val_data62, type = "response")


train_response6 <- tra_data6$LNM.prRLN
test_response6 <- test_data6$LNM.prRLN
val_response61 <- val_data61$LNM.prRLN
val_response62 <- val_data62$LNM.prRLN
# 创建ROC对象
train_roc6 <- roc(train_response6, train_probs6)
test_roc6 <- roc(test_response6, test_probs6)
val_roc61 <- roc(val_response61, val_probs61)
val_roc62 <- roc(val_response62, val_probs62)

# 提取ROC曲线的坐标点
train_roc_data6 <- coords(train_roc6, "all")
test_roc_data6 <- coords(test_roc6, "all")
val_roc_data61 <- coords(val_roc61, "all")
val_roc_data62 <- coords(val_roc62, "all")

# 转换为数据框
train_roc_data6 <- as.data.frame(train_roc_data6)
test_roc_data6 <- as.data.frame(test_roc_data6)
val_roc_data61 <- as.data.frame(val_roc_data61)
val_roc_data62 <- as.data.frame(val_roc_data62)

```


##7.2评价指标
```{r}
# 安装必要的包
install.packages("pROC")
install.packages("caret")
install.packages("Metrics")

# 加载包
library(pROC)
library(caret)
library(Metrics)

```

```{r}
library(caret)
library(pROC)
library(Metrics)

# 定义函数来计算所有指标
compute_metrics <- function(response, probs, threshold = 0.5) {
  predicted <- ifelse(probs >= threshold, 1, 0)
  confusion <- confusionMatrix(factor(predicted), factor(response))
  accuracy <- confusion$overall['Accuracy']
  specificity <- confusion$byClass['Specificity']
  sensitivity <- confusion$byClass['Sensitivity']
  npv <- confusion$byClass['Neg Pred Value']
  ppv <- confusion$byClass['Pos Pred Value']
  f1 <- 2 * (ppv * sensitivity) / (ppv + sensitivity)
  auc <- roc(response, probs)$auc
  false_positive_rate <- 1 - specificity
  brier_score <- mean((probs - response)^2)
  kappa <- confusion$overall['Kappa']
  rmse_value <- rmse(response, probs)
  r2_value <- R2(response, probs)
  mae_value <- mae(response, probs)
  lift <- (sum(predicted) / length(predicted)) / (sum(response) / length(response))
  
  cat("Accuracy:", accuracy, "\n")
  cat("AUC:", auc, "\n")
  cat("Specificity:", specificity, "\n")
  cat("Sensitivity:", sensitivity, "\n")
  cat("NPV:", npv, "\n")
  cat("PPV:", ppv, "\n")
  cat("F1 Score:", f1, "\n")
  cat("False Positive Rate:", false_positive_rate, "\n")
  cat("Brier Score:", brier_score, "\n")
  cat("Kappa:", kappa, "\n")
  cat("RMSE:", rmse_value, "\n")
  cat("R2:", r2_value, "\n")
  cat("MAE:", mae_value, "\n")
  cat("Lift:", lift, "\n")
  
  return(c(accuracy, auc, specificity, sensitivity, npv, ppv, f1, false_positive_rate, lift, brier_score, kappa, rmse_value, r2_value, mae_value))
}

# 计算所有数据集上的指标
metrics_names <- c("Accuracy", "AUC", "Specificity", "Sensitivity/Recall", "Negative Predictive Value", 
                   "Positive Predictive Value/Precision", "F1 Score", "False Positive Rate", "Lift", 
                   "Brier Score", "Kappa", "RMSE", "R2", "MAE")

results <- data.frame(matrix(nrow = 24, ncol = 14))
colnames(results) <- metrics_names

# 训练集
results[1, ] <- compute_metrics(tra_data1$Total.Central.Lymph.Node.Metastasis, train_probs1)
results[2, ] <- compute_metrics(tra_data2$Prelaryngeal.LNM, train_probs2)
results[3, ] <- compute_metrics(tra_data3$Pretracheal.LNM, train_probs3)
results[4, ] <- compute_metrics(tra_data4$Paratracheal.LNM, train_probs4)
results[5, ] <- compute_metrics(tra_data5$Con.Paratracheal.LNM, train_probs5)
results[6, ] <- compute_metrics(tra_data6$LNM.prRLN, train_probs6)

# 测试集
results[7, ] <- compute_metrics(test_data1$Total.Central.Lymph.Node.Metastasis, test_probs1)
results[8, ] <- compute_metrics(test_data2$Prelaryngeal.LNM, test_probs2)
results[9, ] <- compute_metrics(test_data3$Pretracheal.LNM, test_probs3)
results[10, ] <- compute_metrics(test_data4$Paratracheal.LNM, test_probs4)
results[11, ] <- compute_metrics(test_data5$Con.Paratracheal.LNM, test_probs5)
results[12, ] <- compute_metrics(test_data6$LNM.prRLN, test_probs6)

# 验证集1
results[13, ] <- compute_metrics(val_data11$Total.Central.Lymph.Node.Metastasis, val_probs11)
results[14, ] <- compute_metrics(val_data21$Prelaryngeal.LNM, val_probs21)
results[15, ] <- compute_metrics(val_data31$Pretracheal.LNM, val_probs31)
results[16, ] <- compute_metrics(val_data41$Paratracheal.LNM, val_probs41)
results[17, ] <- compute_metrics(val_data51$Con.Paratracheal.LNM, val_probs51)
results[18, ] <- compute_metrics(val_data61$LNM.prRLN, val_probs61)

# 验证集2
results[19, ] <- compute_metrics(val_data12$Total.Central.Lymph.Node.Metastasis, val_probs12)
results[20, ] <- compute_metrics(val_data22$Prelaryngeal.LNM, val_probs22)
results[21, ] <- compute_metrics(val_data32$Pretracheal.LNM, val_probs32)
results[22, ] <- compute_metrics(val_data42$Paratracheal.LNM, val_probs42)
results[23, ] <- compute_metrics(val_data52$Con.Paratracheal.LNM, val_probs52)
results[24, ] <- compute_metrics(val_data62$LNM.prRLN, val_probs62)

# 设置行名
row.names(results) <- c(
  "Total Central Lymph Node Metastasis_train set", "Prelaryngeal LNM_train set", "Pretracheal LNM_train set", 
  "Paratracheal LNM_train set", "Con-Paratracheal LNM_train set", "LNM prRLN_train set", 
  "Total Central Lymph Node Metastasis_Test set", "Prelaryngeal LNM_Test set", "Pretracheal LNM_Test set", 
  "Paratracheal LNM_Test set", "Con-Paratracheal LNM_Test set", "LNM prRLN_Test set", 
  "Total Central Lymph Node Metastasis_validation set1", "Prelaryngeal LNM_validation set1", 
  "Pretracheal LNM_validation set1", "Paratracheal LNM_validation set1", "Con-Paratracheal LNM_validation set1", 
  "LNM prRLN_validation set1", "Total Central Lymph Node Metastasis_validation set2", "Prelaryngeal LNM_validation set2", 
  "Pretracheal LNM_validation set2", "Paratracheal LNM_validation set2", "Con-Paratracheal LNM_validation set2", 
  "LNM prRLN_validation set2"
)

# 导出到CSV文件
write.csv(results, "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/0.评价指标.csv", row.names = TRUE)

```


##7.2总ROC曲线


###7.2.1训练集
```{r}
plot2 <- ggplot() +
  geom_line(data = train_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#9A4942", size = 0.6) +
  geom_line(data = train_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#BB431C", size = 0.6) +
  geom_line(data = train_roc_data3, aes(x = 1 - specificity, y = sensitivity), color = "#C9A51A", size = 0.6) +
  geom_line(data = train_roc_data4, aes(x = 1 - specificity, y = sensitivity), color = "#3D5714", size = 0.6) +
    geom_line(data = train_roc_data5, aes(x = 1 - specificity, y = sensitivity), color = "#82A7D1", size = 0.6) +
  geom_line(data = train_roc_data6, aes(x = 1 - specificity, y = sensitivity), color = "#8833D5", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC Curves for Central Lymph Node Metastasis Nomogram Prediction (Train Set)", 
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.26, label = paste("Total Central Lymph Node Metastasis Train AUC =", round(auc(train_roc1), 3)), size = 3.5, color = "#9A4942") +
  annotate("text", x = 0.7, y = 0.21, label = paste("Prelaryngeal Lymph Node Metastasis Train AUC =", round(auc(train_roc2), 3)), size = 3.5, color = "#BB431C") +
  annotate("text", x = 0.7, y = 0.16, label = paste("Pretracheal Lymph Node Metastasis Train AUC =", round(auc(train_roc3), 3)), size = 3.5, color = "#C9A51A") +
  annotate("text", x = 0.7, y = 0.11, label = paste("Paratracheal Lymph Node Metastasis Train AUC =", round(auc(train_roc4), 3)), size = 3.5, color = "#3D5714")+
  annotate("text", x = 0.7, y = 0.06, label = paste("Con-Paratracheal Lymph Node Metastasis Train AUC =", round(auc(train_roc5), 3)), size = 3.5, color = "#82A7D1") +
  annotate("text", x = 0.6, y = 0.01, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Train AUC =", round(auc(train_roc6), 3)), size = 3.5, color = "#8833D5")



ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/2.1.overall_train-roc_plot.tiff", plot = plot2, dpi = 300, width = 10, height = 8)
print(plot2)

```
###7.2.2测试集
```{r}
plot3 <- ggplot() +
  geom_line(data =test_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#BA3E45", size = 0.6) +
  geom_line(data = test_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#D2431C", size = 0.6) +
  geom_line(data = test_roc_data3, aes(x = 1 - specificity, y = sensitivity), color = "#ECAC27", size = 0.6) +
  geom_line(data = test_roc_data4, aes(x = 1 - specificity, y = sensitivity), color = "#79902D", size = 0.6) +
    geom_line(data = test_roc_data5, aes(x = 1 - specificity, y = sensitivity), color = "#4E6691", size = 0.6) +
  geom_line(data = test_roc_data6, aes(x = 1 - specificity, y = sensitivity), color = "#D355FF", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC Curves for Central Lymph Node Metastasis Nomogram Prediction (Test Set)", 
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.26, label = paste("Total Central Lymph Node Metastasis Test AUC =", round(auc(test_roc1), 3)), size = 3.5, color = "#BA3E45") +
  annotate("text", x = 0.7, y = 0.21, label = paste("Prelaryngeal Lymph Node Metastasis Test AUC =", round(auc(test_roc2), 3)), size = 3.5, color = "#D2431C") +
  annotate("text", x = 0.7, y = 0.16, label = paste("Pretracheal Lymph Node Metastasis Test AUC =", round(auc(test_roc3), 3)), size = 3.5, color = "#ECAC27") +
  annotate("text", x = 0.7, y = 0.11, label = paste("Paratracheal Lymph Node Metastasis Test AUC =", round(auc(test_roc4), 3)), size = 3.5, color = "#79902D")+
  annotate("text", x = 0.7, y = 0.06, label = paste("Con-Paratracheal Lymph Node Metastasis Test AUC =", round(auc(test_roc5), 3)), size = 3.5, color = "#4E6691") +
  annotate("text", x = 0.6, y = 0.01, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Test AUC =", round(auc(test_roc6), 3)), size = 3.5, color = "#D355FF")



ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/2.2.overall_test-roc_plot.tiff", plot = plot3, dpi = 300, width = 10, height = 8)
print(plot3)


```
###7.2.3验证集1
```{r}
plot4 <- ggplot() +
  geom_line(data =val_roc_data11, aes(x = 1 - specificity, y = sensitivity), color = "#EABFBB", size = 0.6) +
  geom_line(data = val_roc_data21, aes(x = 1 - specificity, y = sensitivity), color = "#F2AB6A", size = 0.6) +
  geom_line(data = val_roc_data31, aes(x = 1 - specificity, y = sensitivity), color = "#EDDE23", size = 0.6) +
  geom_line(data = val_roc_data41, aes(x = 1 - specificity, y = sensitivity), color = "#5AB682", size = 0.6) +
    geom_line(data = val_roc_data51, aes(x = 1 - specificity, y = sensitivity), color = "#B6D7E9", size = 0.6) +
  geom_line(data = val_roc_data61, aes(x = 1 - specificity, y = sensitivity), color = "#E8A4FF", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC Curves for Central Lymph Node Metastasis Nomogram Prediction (Validation Set1)", 
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.26, label = paste("Total Central Lymph Node Metastasis Validatio1 AUC =", round(auc(val_roc11), 3)), size = 3.5, color = "#EABFBB") +
  annotate("text", x = 0.7, y = 0.21, label = paste("Prelaryngeal Lymph Node Metastasis Validatio1 AUC =", round(auc(val_roc21), 3)), size = 3.5, color = "#F2AB6A") +
  annotate("text", x = 0.7, y = 0.16, label = paste("Pretracheal Lymph Node Metastasis Validatio1 AUC =", round(auc(val_roc31), 3)), size = 3.5, color = "#EDDE23") +
  annotate("text", x = 0.7, y = 0.11, label = paste("Paratracheal Lymph Node Metastasis Validatio1 AUC =", round(auc(val_roc41), 3)), size = 3.5, color = "#5AB682")+
  annotate("text", x = 0.7, y = 0.06, label = paste("Con-Paratracheal Lymph Node Metastasis Validatio1 AUC =", round(auc(val_roc51), 3)), size = 3.5, color = "#B6D7E9") +
  annotate("text", x = 0.6, y = 0.01, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Validatio1 AUC =", round(auc(val_roc61), 3)), size = 3.5, color = "#E8A4FF")



ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/2.3.overall_val1-roc_plot.tiff", plot = plot4, dpi = 300, width = 10, height = 8)
print(plot4)

```

###7.2.4验证集2
```{r}
plot4 <- ggplot() +
  geom_line(data =val_roc_data12, aes(x = 1 - specificity, y = sensitivity), color = "#EAB", size = 0.6) +
  geom_line(data = val_roc_data22, aes(x = 1 - specificity, y = sensitivity), color = "#F5D18B", size = 0.6) +
  geom_line(data = val_roc_data32, aes(x = 1 - specificity, y = sensitivity), color = "#FFFF66", size = 0.6) +
  geom_line(data = val_roc_data42, aes(x = 1 - specificity, y = sensitivity), color = "#C8E4D2", size = 0.6) +
    geom_line(data = val_roc_data52, aes(x = 1 - specificity, y = sensitivity), color = "#DBEAF3", size = 0.6) +
  geom_line(data = val_roc_data62, aes(x = 1 - specificity, y = sensitivity), color = "#F0CCFF", size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC Curves for Central Lymph Node Metastasis Nomogram Prediction (Validation Set2)", 
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.26, label = paste("Total Central Lymph Node Metastasis Validatio2 AUC =", round(auc(val_roc12), 3)), size = 3.5, color = "#EAB") +
  annotate("text", x = 0.7, y = 0.21, label = paste("Prelaryngeal Lymph Node Metastasis Validatio2 AUC =", round(auc(val_roc22), 3)), size = 3.5, color = "#F2AB6A") +
  annotate("text", x = 0.7, y = 0.16, label = paste("Pretracheal Lymph Node Metastasis Validatio2 AUC =", round(auc(val_roc32), 3)), size = 3.5, color = "#FFFF66") +
  annotate("text", x = 0.7, y = 0.11, label = paste("Paratracheal Lymph Node Metastasis Validatio2 AUC =", round(auc(val_roc42), 3)), size = 3.5, color = "#C8E4D2")+
  annotate("text", x = 0.7, y = 0.06, label = paste("Con-Paratracheal Lymph Node Metastasis Validatio2 AUC =", round(auc(val_roc52), 3)), size = 3.5, color = "#DBEAF3") +
  annotate("text", x = 0.6, y = 0.01, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Validatio2 AUC =", round(auc(val_roc62), 3)), size = 3.5, color = "#F0CCFF")



ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/2.4.overall_val2-roc_plot.tiff", plot = plot4, dpi = 300, width = 10, height = 8)
print(plot4)

```


###7.2.5全24rco曲线
```{r}
library(ggplot2)


# 生成整体的图形
plot1 <- ggplot() +
  geom_line(data = train_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#9A4942", size = 0.6) +
  geom_line(data = test_roc_data1, aes(x = 1 - specificity, y = sensitivity), color = "#BA3E45", size = 0.6) +
  geom_line(data = val_roc_data11, aes(x = 1 - specificity, y = sensitivity), color = "#EABFBB", size = 0.6) +
  geom_line(data = val_roc_data12, aes(x = 1 - specificity, y = sensitivity), color = "#EAB", size = 0.6) +
  
  geom_line(data = train_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#BB431C", size = 0.6) +
  geom_line(data = test_roc_data2, aes(x = 1 - specificity, y = sensitivity), color = "#D2431C", size = 0.6) +
  geom_line(data = val_roc_data21, aes(x = 1 - specificity, y = sensitivity), color = "#F2AB6A", size = 0.6) +
  geom_line(data = val_roc_data22, aes(x = 1 - specificity, y = sensitivity), color = "#F5D18B", size = 0.6) +
  
  geom_line(data = train_roc_data3, aes(x = 1 - specificity, y = sensitivity), color = "#C9A51A", size = 0.6) +
  geom_line(data = test_roc_data3, aes(x = 1 - specificity, y = sensitivity), color = "#ECAC27", size = 0.6) +
  geom_line(data = val_roc_data31, aes(x = 1 - specificity, y = sensitivity), color = "#EDDE23", size = 0.6) +
  geom_line(data = val_roc_data32, aes(x = 1 - specificity, y = sensitivity), color = "#FFFF66", size = 0.6) +
  
  geom_line(data = train_roc_data4, aes(x = 1 - specificity, y = sensitivity), color = "#3D5714", size = 0.6) +
  geom_line(data = test_roc_data4, aes(x = 1 - specificity, y = sensitivity), color = "#79902D", size = 0.6) +
  geom_line(data = val_roc_data41, aes(x = 1 - specificity, y = sensitivity), color = "#5AB682", size = 0.6) +
  geom_line(data = val_roc_data42, aes(x = 1 - specificity, y = sensitivity), color = "#C8E4D2", size = 0.6) +
  
  geom_line(data = train_roc_data5, aes(x = 1 - specificity, y = sensitivity), color = "#82A7D1", size = 0.6) +
  geom_line(data = test_roc_data5, aes(x = 1 - specificity, y = sensitivity), color = "#4E6691", size = 0.6) +
  geom_line(data = val_roc_data51, aes(x = 1 - specificity, y = sensitivity), color = "#B6D7E9", size = 0.6) +
  geom_line(data = val_roc_data52, aes(x = 1 - specificity, y = sensitivity), color = "#DBEAF3", size = 0.6) +
  
  geom_line(data = train_roc_data6, aes(x = 1 - specificity, y = sensitivity), color = "#8833D5", size = 0.6) +
  geom_line(data = test_roc_data6, aes(x = 1 - specificity, y = sensitivity), color = "#D355FF", size = 0.6) +
  geom_line(data = val_roc_data61, aes(x = 1 - specificity, y = sensitivity), color = "#E8A4FF", size = 0.6) +
  geom_line(data = val_roc_data62, aes(x = 1 - specificity, y = sensitivity), color = "#F0CCFF", size = 0.6) + 
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC Curves for Central Lymph Node Metastasis Nomogram Prediction (All Sets)", size = 3,
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = 0.8, y = 0.72, label = paste("Total Central Lymph Node Metastasis Train AUC =", round(auc(train_roc1), 3)), size = 3, color = "#9A4942") +
  annotate("text", x = 0.8, y = 0.69, label = paste("Total Central Lymph Node Metastasis Test AUC =", round(auc(test_roc1), 3)), size = 3, color = "#BA3E45") +
  annotate("text", x = 0.8, y = 0.66, label = paste("Total Central Lymph Node Metastasis Validation1 AUC =", round(auc(val_roc11), 3)), size = 3, color = "#EABFBB") +
  annotate("text", x = 0.8, y = 0.63, label = paste("Total Central Lymph Node Metastasis Validation2 AUC =", round(auc(val_roc12), 3)), size = 3, color = "#EAB") +
  
  annotate("text", x = 0.8, y = 0.60, label = paste("Prelaryngeal Lymph Node Metastasis Train AUC =", round(auc(train_roc2), 3)), size = 3, color = "#BB431C") +
  annotate("text", x = 0.8, y = 0.57, label = paste("Prelaryngeal Lymph Node Metastasis Test AUC =", round(auc(test_roc2), 3)), size = 3, color = "#D2431C") +
  annotate("text", x = 0.8, y = 0.54, label = paste("Prelaryngeal Lymph Node Metastasis Validation1 AUC =", round(auc(val_roc21), 3)), size = 3, color = "#F2AB6A") +
  annotate("text", x = 0.8, y = 0.51, label = paste("Prelaryngeal Lymph Node Metastasis Validation2 AUC =", round(auc(val_roc22), 3)), size = 3, color = "#F5D18B") +
  
  annotate("text", x = 0.8, y = 0.48, label = paste("Pretracheal Lymph Node Metastasis Train AUC =", round(auc(train_roc3), 3)), size = 3, color = "#C9A51A") +
  annotate("text", x = 0.8, y = 0.45, label = paste("Pretracheal Lymph Node Metastasis Test AUC =", round(auc(test_roc3), 3)), size = 3, color = "#ECAC27") +
  annotate("text", x = 0.8, y = 0.42, label = paste("Pretracheal Lymph Node Metastasis Validation1 AUC =", round(auc(val_roc31), 3)), size = 3, color = "#EDDE23") + 
  annotate("text", x = 0.8, y = 0.39, label = paste("Pretracheal Lymph Node Metastasis Validation2 AUC =", round(auc(val_roc32), 3)), size = 3, color = "#FFFF66") + 
  
  annotate("text", x = 0.8, y = 0.36, label = paste("Paratracheal Lymph Node Metastasis Train AUC =", round(auc(train_roc4), 3)), size = 3, color = "#3D5714") +
  annotate("text", x = 0.8, y = 0.33, label = paste("Paratracheal Lymph Node Metastasis Test AUC =", round(auc(test_roc4), 3)), size = 3, color = "#79902D") + 
  annotate("text", x = 0.8, y = 0.3, label = paste("Paratracheal Lymph Node Metastasis Validation1 AUC =", round(auc(val_roc41), 3)), size = 3, color = "#5AB682") + 
  annotate("text", x = 0.8, y = 0.27, label = paste("Paratracheal Lymph Node Metastasis Validation2 AUC =", round(auc(val_roc42), 3)), size = 3, color = "#C8E4D2") +
  
  annotate("text", x = 0.8, y = 0.24, label = paste("Con-Paratracheal Lymph Node Metastasis Train AUC =", round(auc(train_roc5), 3)), size = 3, color = "#82A7D1") +
  annotate("text", x = 0.8, y = 0.21, label = paste("Con-Paratracheal Lymph Node Metastasis Test AUC =", round(auc(test_roc5), 3)), size = 3, color = "#4E6691") + 
  annotate("text", x = 0.8, y = 0.18, label = paste("Con-Paratracheal Lymph Node Metastasis Validation1 AUC =", round(auc(val_roc51), 3)), size = 3, color = "#B6D7E9") + 
  annotate("text", x = 0.8, y = 0.15, label = paste("Con-Paratracheal Lymph Node Metastasis Validation2 AUC =", round(auc(val_roc52), 3)), size = 3, color = "#DBEAF3") +
  
  annotate("text", x = 0.8, y = 0.11, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Train AUC =", round(auc(train_roc6), 3)), size = 3, color = "#8833D5") +
  annotate("text", x = 0.8, y = 0.07, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis Test AUC =", round(auc(test_roc6), 3)), size = 3, color = "#D355FF") + 
  annotate("text", x = 0.8, y = 0.04, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis1 Validation AUC =", round(auc(val_roc61), 3)), size = 3, color = "#E8A4FF") + 
  annotate("text", x = 0.8, y = 0.01, label = paste("Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis2 Validation AUC =", round(auc(val_roc62), 3)), size = 3, color = "#F0CCFF")

ggsave("/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/2.5.overall_roc_plot.tiff", plot = plot1, dpi = 300, width = 10, height = 8)

print(plot1)


```




##7.3总的DCA曲线
####把图标注放置在图里面适合曲线较少的
```{r}
# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits1 <- sapply(thresholds, function(x) net_benefit(train_probs1, train_response1, x))
train_net_benefits2 <- sapply(thresholds, function(x) net_benefit(train_probs2, train_response2, x))
train_net_benefits3 <- sapply(thresholds, function(x) net_benefit(train_probs3, train_response3, x))
train_net_benefits4 <- sapply(thresholds, function(x) net_benefit(train_probs4, train_response4, x))
train_net_benefits5 <- sapply(thresholds, function(x) net_benefit(train_probs3, train_response5, x))
train_net_benefits6 <- sapply(thresholds, function(x) net_benefit(train_probs4, train_response6, x))
# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(train_response6)), train_response6, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 绘制DCA曲线
dca_data_train <- data.frame(
  threshold = thresholds,
  train_net_benefit1 = train_net_benefits1,
  train_net_benefit2 = train_net_benefits2,
  train_net_benefit3 = train_net_benefits3,
  train_net_benefit4 = train_net_benefits4,
  train_net_benefit5 = train_net_benefits5,
  train_net_benefit6 = train_net_benefits6,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data_train, aes(x = threshold)) + 
  geom_line(aes(y = train_net_benefit1, color = "TCLNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit2, color = "Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit3, color = "Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit4, color = "Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit5, color = "Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit6, color = "LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Central Lymph Node Metastasis Nomogram Prediction (Train Set)",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(
    breaks = c("TCLNM", "Prelaryngeal LNM", "Pretracheal LNM", "Paratracheal LNM", "Con-Paratracheal LNM", "LNM-prRLN", "All", "None"),
    values = c("TCLNM" = "#9A4942", 
               "Prelaryngeal LNM" = "#BB431C", 
               "Pretracheal LNM" = "#C9A51A", 
               "Paratracheal LNM" = "#3D5714", 
               "Con-Paratracheal LNM" = "#82A7D1", 
               "LNM-prRLN" = "#8833D5",                                
               "All" = "grey",
               "None" = "black")
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(-0.05, 0.45), xlim = c(0, 0.8)) +
  annotate("text", x = 0.6, y = 0.4, label = "TCLNM", color = "#9A4942", size = 3.5) +
  annotate("text", x = 0.6, y = 0.35, label = "Prelaryngeal LNM", color = "#BB431C", size = 3.5) +
  annotate("text", x = 0.6, y = 0.3, label = "Pretracheal LNM", color = "#C9A51A", size = 3.5) +
  annotate("text", x = 0.6, y = 0.25, label = "Paratracheal LNM", color = "#3D5714", size = 3.5) +
  annotate("text", x = 0.6, y = 0.2, label = "Con-Paratracheal LNM", color = "#82A7D1", size = 3.5) +
  annotate("text", x = 0.6, y = 0.15, label = "LNM-prRLN", color = "#8833D5", size = 3.5) +
  annotate("text", x = 0.6, y = 0.1, label = "All", color = "grey", size = 3.5) +
  annotate("text", x = 0.6, y = 0.05, label = "None", color = "black", size = 3.5)

# 保存图像函数
save_dca_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

# 保存DCA曲线图
save_dca_plot(dca_plot, file.path(output_folder, "3.1.overall_train-dca_plot.tiff"))

print(dca_plot)
```

###7.3.1训练集
```{r}
# 定义净收益计算函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 计算不同阈值下的净收益
thresholds <- seq(0, 1, by = 0.01)
train_net_benefits1 <- sapply(thresholds, function(x) net_benefit(train_probs1, train_response1, x))
train_net_benefits2 <- sapply(thresholds, function(x) net_benefit(train_probs2, train_response2, x))
train_net_benefits3 <- sapply(thresholds, function(x) net_benefit(train_probs3, train_response3, x))
train_net_benefits4 <- sapply(thresholds, function(x) net_benefit(train_probs4, train_response4, x))
train_net_benefits5 <- sapply(thresholds, function(x) net_benefit(train_probs3, train_response5, x))
train_net_benefits6 <- sapply(thresholds, function(x) net_benefit(train_probs4, train_response6, x))
# 计算所有人都进行干预时的净收益
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(train_response6)), train_response6, x))

# 计算没有人进行干预时的净收益
none_net_benefit <- rep(0, length(thresholds))

# 绘制DCA曲线
dca_data_train <- data.frame(
  threshold = thresholds,
  train_net_benefit1 = train_net_benefits1,
  train_net_benefit2 = train_net_benefits2,
  train_net_benefit3 = train_net_benefits3,
  train_net_benefit4 = train_net_benefits4,
  train_net_benefit5 = train_net_benefits5,
  train_net_benefit6 = train_net_benefits6,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

dca_plot <- ggplot(dca_data_train, aes(x = threshold)) + 
  geom_line(aes(y = train_net_benefit1, color = "TCLNM"), size = 0.6)+
  geom_line(aes(y = train_net_benefit2, color = "Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit3, color = "Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit4, color = "Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit5, color = "Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit6, color = "LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Central Lymph Node Metastasis Nomogram Prediction (Train Set)",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(
    breaks = c("TCLNM", "Prelaryngeal LNM", "Pretracheal LNM", "Paratracheal LNM", "Con-Paratracheal LNM", "LNM-prRLN", "All", "None"),
    values = c("TCLNM" = "#9A4942", 
               "Prelaryngeal LNM" = "#BB431C", 
               "Pretracheal LNM" = "#C9A51A", 
               "Paratracheal LNM" = "#3D5714", 
               "Con-Paratracheal LNM" = "#82A7D1", 
               "LNM-prRLN" = "#8833D5",                                
               "All" = "grey",
               "None" = "black")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.05, 0.45), xlim = c(0, 0.8))

# 保存图像函数
save_dca_plot <- function(plot, filename_prefix) {
  # Save as TIFF at 300 DPI
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

# 保存DCA曲线图
save_dca_plot(dca_plot, file.path(output_folder, "3.1.overall_train-dca_plot.tiff"))

print(dca_plot)
```

###7.3.2测试集
```{r}
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

thresholds <- seq(0, 1, by = 0.01)
test_net_benefits1 <- sapply(thresholds, function(x) net_benefit(test_probs1, test_response1, x))
test_net_benefits2 <- sapply(thresholds, function(x) net_benefit(test_probs2, test_response2, x))
test_net_benefits3 <- sapply(thresholds, function(x) net_benefit(test_probs3, test_response3, x))
test_net_benefits4 <- sapply(thresholds, function(x) net_benefit(test_probs4, test_response4, x))
test_net_benefits5 <- sapply(thresholds, function(x) net_benefit(test_probs5, test_response5, x))
test_net_benefits6 <- sapply(thresholds, function(x) net_benefit(test_probs6, test_response6, x))
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(test_response6)), test_response6, x))
none_net_benefit <- rep(0, length(thresholds))

dca_data_test <- data.frame(
  threshold = thresholds,
  test_net_benefit1 = test_net_benefits1,
  test_net_benefit2 = test_net_benefits2,
  test_net_benefit3 = test_net_benefits3,
  test_net_benefit4 = test_net_benefits4,
  test_net_benefit5 = test_net_benefits5,
  test_net_benefit6 = test_net_benefits6,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

library(ggplot2)

dca_plot_test <- ggplot(dca_data_test, aes(x = threshold)) + 
  geom_line(aes(y = test_net_benefit1, color = "TCLNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit2, color = "Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit3, color = "Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit4, color = "Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit5, color = "Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit6, color = "LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Central Lymph Node Metastasis Nomogram Prediction (Test Set)",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(
    breaks = c("TCLNM", "Prelaryngeal LNM", "Pretracheal LNM", "Paratracheal LNM", "Con-Paratracheal LNM", "LNM-prRLN", "All", "None"),
    values = c("TCLNM" = "#BA3E45", 
               "Prelaryngeal LNM" = "#D2431C", 
               "Pretracheal LNM" = "#ECAC27", 
               "Paratracheal LNM" = "#79902D", 
               "Con-Paratracheal LNM" = "#4E6691", 
               "LNM-prRLN" = "#D355FF",                                
               "All" = "grey",
               "None" = "black")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.05, 0.45), xlim = c(0, 0.8))

# 保存图像
save_dca_plot <- function(plot, filename_prefix) {
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

# 保存 DCA 曲线图
save_dca_plot(dca_plot_test, file.path(output_folder, "3.2.overall_test-dca_plot"))

print(dca_plot_test)

```
###7.3.3验证集
```{r}
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

thresholds <- seq(0, 1, by = 0.01)
val1_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs11, val_response11, x))
val1_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs21, val_response21, x))
val1_net_benefits3 <- sapply(thresholds, function(x) net_benefit(val_probs31, val_response31, x))
val1_net_benefits4 <- sapply(thresholds, function(x) net_benefit(val_probs41, val_response41, x))
val1_net_benefits5 <- sapply(thresholds, function(x) net_benefit(val_probs51, val_response51, x))
val1_net_benefits6 <- sapply(thresholds, function(x) net_benefit(val_probs61, val_response61, x))
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(val_response61)), val_response61, x))
none_net_benefit <- rep(0, length(thresholds))

dca_data_val1 <- data.frame(
  threshold = thresholds,
  val1_net_benefit1 = val1_net_benefits1,
  val1_net_benefit2 = val1_net_benefits2,
  val1_net_benefit3 = val1_net_benefits3,
  val1_net_benefit4 = val1_net_benefits4,
  val1_net_benefit5 = val1_net_benefits5,
  val1_net_benefit6 = val1_net_benefits6,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

library(ggplot2)

dca_plot_val1 <- ggplot(dca_data_val1, aes(x = threshold)) + 
  geom_line(aes(y = val1_net_benefit1, color = "TCLNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit2, color = "Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit3, color = "Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit4, color = "Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit5, color = "Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit6, color = "LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Central Lymph Node Metastasis Nomogram Prediction (Validation Set 1)",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(
    breaks = c("TCLNM", "Prelaryngeal LNM", "Pretracheal LNM", "Paratracheal LNM", "Con-Paratracheal LNM", "LNM-prRLN", "All", "None"),
    values = c("TCLNM" = "#EABFBB", 
               "Prelaryngeal LNM" = "#F2AB6A", 
               "Pretracheal LNM" = "#EDDE23", 
               "Paratracheal LNM" = "#5AB682", 
               "Con-Paratracheal LNM" = "#B6D7E9", 
               "LNM-prRLN" = "#E8A4FF",                                
               "All" = "grey",
               "None" = "black")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.05, 0.45), xlim = c(0, 0.8))

# 保存图像
save_dca_plot <- function(plot, filename_prefix) {
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

# 保存 DCA 曲线图
save_dca_plot(dca_plot_val1, file.path(output_folder, "3.3.overall_val1-dca_plot"))

print(dca_plot_val1)


```

###7.3.4验证集2
```{r}
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

thresholds <- seq(0, 1, by = 0.01)
val2_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs12, val_response12, x))
val2_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs22, val_response22, x))
val2_net_benefits3 <- sapply(thresholds, function(x) net_benefit(val_probs32, val_response32, x))
val2_net_benefits4 <- sapply(thresholds, function(x) net_benefit(val_probs42, val_response42, x))
val2_net_benefits5 <- sapply(thresholds, function(x) net_benefit(val_probs52, val_response52, x))
val2_net_benefits6 <- sapply(thresholds, function(x) net_benefit(val_probs62, val_response62, x))
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(val_response62)), val_response62, x))
none_net_benefit <- rep(0, length(thresholds))

dca_data_val2 <- data.frame(
  threshold = thresholds,
  val2_net_benefit1 = val2_net_benefits1,
  val2_net_benefit2 = val2_net_benefits2,
  val2_net_benefit3 = val2_net_benefits3,
  val2_net_benefit4 = val2_net_benefits4,
  val2_net_benefit5 = val2_net_benefits5,
  val2_net_benefit6 = val2_net_benefits6,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

library(ggplot2)

dca_plot_val2 <- ggplot(dca_data_val2, aes(x = threshold)) + 
  geom_line(aes(y = val2_net_benefit1, color = "TCLNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit2, color = "Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit3, color = "Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit4, color = "Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit5, color = "Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit6, color = "LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Central Lymph Node Metastasis Nomogram Prediction (Validation Set 2)",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(
    breaks = c("TCLNM", "Prelaryngeal LNM", "Pretracheal LNM", "Paratracheal LNM", "Con-Paratracheal LNM", "LNM-prRLN", "All", "None"),
    values = c("TCLNM" = "#EAB", 
               "Prelaryngeal LNM" = "#F5D18B", 
               "Pretracheal LNM" = "#FFFF66", 
               "Paratracheal LNM" = "#C8E4D2", 
               "Con-Paratracheal LNM" = "#DBEAF3", 
               "LNM-prRLN" = "#F0CCFF",                                
               "All" = "grey",
               "None" = "black")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.05, 0.45), xlim = c(0, 0.8))

# 保存图像
save_dca_plot <- function(plot, filename_prefix) {
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

# 保存 DCA 曲线图
save_dca_plot(dca_plot_val2, file.path(output_folder, "3.4.overall_val2-dca_plot"))

print(dca_plot_val2)
```


###7.3.5全24dca曲线
```{r}
# 计算net_benefit函数
net_benefit <- function(probs, outcome, threshold) {
  tp <- sum(outcome == 1 & probs >= threshold)
  fp <- sum(outcome == 0 & probs >= threshold)
  total_population <- length(outcome)
  
  if (total_population == 0) {
    return(0)
  }
  
  net_benefit <- (tp / total_population) - ((fp / total_population) * (threshold / (1 - threshold)))
  return(net_benefit)
}

# 定义阈值范围
thresholds <- seq(0, 1, by = 0.01)

# 计算训练集的net_benefit
train_net_benefits1 <- sapply(thresholds, function(x) net_benefit(train_probs1, train_response1, x))
train_net_benefits2 <- sapply(thresholds, function(x) net_benefit(train_probs2, train_response2, x))
train_net_benefits3 <- sapply(thresholds, function(x) net_benefit(train_probs3, train_response3, x))
train_net_benefits4 <- sapply(thresholds, function(x) net_benefit(train_probs4, train_response4, x))
train_net_benefits5 <- sapply(thresholds, function(x) net_benefit(train_probs5, train_response5, x))
train_net_benefits6 <- sapply(thresholds, function(x) net_benefit(train_probs6, train_response6, x))

# 计算测试集的net_benefit
test_net_benefits1 <- sapply(thresholds, function(x) net_benefit(test_probs1, test_response1, x))
test_net_benefits2 <- sapply(thresholds, function(x) net_benefit(test_probs2, test_response2, x))
test_net_benefits3 <- sapply(thresholds, function(x) net_benefit(test_probs3, test_response3, x))
test_net_benefits4 <- sapply(thresholds, function(x) net_benefit(test_probs4, test_response4, x))
test_net_benefits5 <- sapply(thresholds, function(x) net_benefit(test_probs5, test_response5, x))
test_net_benefits6 <- sapply(thresholds, function(x) net_benefit(test_probs6, test_response6, x))

# 计算验证集1的net_benefit
val1_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs11, val_response11, x))
val1_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs21, val_response21, x))
val1_net_benefits3 <- sapply(thresholds, function(x) net_benefit(val_probs31, val_response31, x))
val1_net_benefits4 <- sapply(thresholds, function(x) net_benefit(val_probs41, val_response41, x))
val1_net_benefits5 <- sapply(thresholds, function(x) net_benefit(val_probs51, val_response51, x))
val1_net_benefits6 <- sapply(thresholds, function(x) net_benefit(val_probs61, val_response61, x))

# 计算验证集2的net_benefit
val2_net_benefits1 <- sapply(thresholds, function(x) net_benefit(val_probs12, val_response12, x))
val2_net_benefits2 <- sapply(thresholds, function(x) net_benefit(val_probs22, val_response22, x))
val2_net_benefits3 <- sapply(thresholds, function(x) net_benefit(val_probs32, val_response32, x))
val2_net_benefits4 <- sapply(thresholds, function(x) net_benefit(val_probs42, val_response42, x))
val2_net_benefits5 <- sapply(thresholds, function(x) net_benefit(val_probs52, val_response52, x))
val2_net_benefits6 <- sapply(thresholds, function(x) net_benefit(val_probs62, val_response62, x))

# 计算所有和无的net_benefit
all_net_benefit <- sapply(thresholds, function(x) net_benefit(rep(1, length(val_response62)), val_response62, x))
none_net_benefit <- rep(0, length(thresholds))

# 合并所有数据
dca_data <- data.frame(
  threshold = thresholds,
  train_net_benefit1 = train_net_benefits1,
  train_net_benefit2 = train_net_benefits2,
  train_net_benefit3 = train_net_benefits3,
  train_net_benefit4 = train_net_benefits4,
  train_net_benefit5 = train_net_benefits5,
  train_net_benefit6 = train_net_benefits6,
  test_net_benefit1 = test_net_benefits1,
  test_net_benefit2 = test_net_benefits2,
  test_net_benefit3 = test_net_benefits3,
  test_net_benefit4 = test_net_benefits4,
  test_net_benefit5 = test_net_benefits5,
  test_net_benefit6 = test_net_benefits6,
  val1_net_benefit1 = val1_net_benefits1,
  val1_net_benefit2 = val1_net_benefits2,
  val1_net_benefit3 = val1_net_benefits3,
  val1_net_benefit4 = val1_net_benefits4,
  val1_net_benefit5 = val1_net_benefits5,
  val1_net_benefit6 = val1_net_benefits6,
  val2_net_benefit1 = val2_net_benefits1,
  val2_net_benefit2 = val2_net_benefits2,
  val2_net_benefit3 = val2_net_benefits3,
  val2_net_benefit4 = val2_net_benefits4,
  val2_net_benefit5 = val2_net_benefits5,
  val2_net_benefit6 = val2_net_benefits6,
  all_net_benefit = all_net_benefit,
  none_net_benefit = none_net_benefit
)

library(ggplot2)

# 生成 DCA 曲线图
dca_plot <- ggplot(dca_data, aes(x = threshold)) + 
  geom_line(aes(y = train_net_benefit1, color = "Train TCLNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit2, color = "Train Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit3, color = "Train Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit4, color = "Train Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit5, color = "Train Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = train_net_benefit6, color = "Train LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = test_net_benefit1, color = "Test TCLNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit2, color = "Test Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit3, color = "Test Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit4, color = "Test Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit5, color = "Test Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = test_net_benefit6, color = "Test LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit1, color = "Val1 TCLNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit2, color = "Val1 Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit3, color = "Val1 Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit4, color = "Val1 Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit5, color = "Val1 Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val1_net_benefit6, color = "Val1 LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit1, color = "Val2 TCLNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit2, color = "Val2 Prelaryngeal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit3, color = "Val2 Pretracheal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit4, color = "Val2 Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit5, color = "Val2 Con-Paratracheal LNM"), size = 0.6) +
  geom_line(aes(y = val2_net_benefit6, color = "Val2 LNM-prRLN"), size = 0.6) +
  geom_line(aes(y = all_net_benefit, color = "All"), linetype = "dotted", size = 0.6) +
  geom_line(aes(y = none_net_benefit, color = "None"), linetype = "solid", size = 0.6) +
  labs(title = "DCA for Central Lymph Node Metastasis Nomogram Prediction",
       x = "Threshold Probability", y = "Net Benefit") +
  scale_color_manual(
    breaks = c("Train TCLNM", "Train Prelaryngeal LNM", "Train Pretracheal LNM", "Train Paratracheal LNM", "Train Con-Paratracheal LNM", "Train LNM-prRLN",
               "Test TCLNM", "Test Prelaryngeal LNM", "Test Pretracheal LNM", "Test Paratracheal LNM", "Test Con-Paratracheal LNM", "Test LNM-prRLN",
               "Val1 TCLNM", "Val1 Prelaryngeal LNM", "Val1 Pretracheal LNM", "Val1 Paratracheal LNM", "Val1 Con-Paratracheal LNM", "Val1 LNM-prRLN",
               "Val2 TCLNM", "Val2 Prelaryngeal LNM", "Val2 Pretracheal LNM", "Val2 Paratracheal LNM", "Val2 Con-Paratracheal LNM", "Val2 LNM-prRLN",
               "All", "None"),
    values = c(
      "Train TCLNM" = "#9A4942", 
      "Train Prelaryngeal LNM" = "#BB431C", 
      "Train Pretracheal LNM" = "#C9A51A", 
      "Train Paratracheal LNM" = "#3D5714", 
      "Train Con-Paratracheal LNM" = "#82A7D1", 
      "Train LNM-prRLN" = "#8833D5",
      "Test TCLNM" = "#BA3E45", 
      "Test Prelaryngeal LNM" = "#D2431C", 
      "Test Pretracheal LNM" = "#ECAC27", 
      "Test Paratracheal LNM" = "#79902D", 
      "Test Con-Paratracheal LNM" = "#4E6691", 
      "Test LNM-prRLN" = "#D355FF",
      "Val1 TCLNM" = "#EABFBB", 
      "Val1 Prelaryngeal LNM" = "#F2AB6A", 
      "Val1 Pretracheal LNM" = "#EDDE23", 
      "Val1 Paratracheal LNM" = "#5AB682", 
      "Val1 Con-Paratracheal LNM" = "#B6D7E9", 
      "Val1 LNM-prRLN" = "#E8A4FF",
      "Val2 TCLNM" = "#EAB", 
      "Val2 Prelaryngeal LNM" = "#F5D18B", 
      "Val2 Pretracheal LNM" = "#FFFF66", 
      "Val2 Paratracheal LNM" = "#CBE4D2", 
      "Val2 Con-Paratracheal LNM" = "#DBEAF3", 
      "Val2 LNM-prRLN" = "#F0CCFF",
      "All" = "grey",
      "None" = "black"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(-0.05, 0.45), xlim = c(0, 0.8))

# 保存图像
save_dca_plot <- function(plot, filename_prefix) {
  ggsave(paste0(filename_prefix, "_300dpi.tiff"), plot, dpi = 300, width = 10, height = 10)
}

# 文件夹路径
output_folder <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

# 保存 DCA 曲线图
save_dca_plot(dca_plot, file.path(output_folder, "3.5.overall_dca_plot"))

print(dca_plot)


```


##7.4学习曲线

###7.4.1加载及划分数据
```{r}
# 加载必要的库
library(pROC)
library(ggplot2)
library(dplyr)


# 读取数据
train_data1 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总编码后_插补.csv")
val_data11 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1编码后_插补.csv")
val_data12 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2编码后_插补.csv")
# 读取数据
train_data2 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总H编码后_插补.csv")
val_data21 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1H编码后_插补.csv")
val_data22 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2H编码后_插补.csv")
# 读取数据
train_data3 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总Q编码后_插补.csv")
val_data31 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1Q编码后_插补.csv")
val_data32 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2Q编码后_插补.csv")
# 读取数据
train_data4 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总P编码后_插补.csv")
val_data41 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1P编码后_插补.csv")
val_data42 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2P编码后_插补.csv")
# 读取数据
train_data5 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总CP编码后_插补.csv")
val_data51 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1CP编码后_插补.csv")
val_data52 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2CP编码后_插补.csv")
# 读取数据
train_data6 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv")
val_data61 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1F编码后_插补.csv")
val_data62 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2F编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data1)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size1 <- round(tra_ratio * nrow(train_data1))
tra_data1 <- train_data1[shuffled_index[1:tra_size1], ]
test_data1 <- train_data1[shuffled_index[(tra_size1 + 1):nrow(train_data1)], ]

cat("训练集观测数量:", nrow(tra_data1), "\n")
cat("测试集观测数量:", nrow(test_data1), "\n")
cat("验证集1观测数量:", nrow(val_data11), "\n")
cat("验证集2观测数量:", nrow(val_data12), "\n")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data2)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size2 <- round(tra_ratio * nrow(train_data2))
tra_data2 <- train_data2[shuffled_index[1:tra_size2], ]
test_data2 <- train_data2[shuffled_index[(tra_size2 + 1):nrow(train_data2)], ]

cat("训练集观测数量:", nrow(tra_data2), "\n")
cat("测试集观测数量:", nrow(test_data2), "\n")
cat("验证集1观测数量:", nrow(val_data21), "\n")
cat("验证集2观测数量:", nrow(val_data22), "\n")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data3)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size3 <- round(tra_ratio * nrow(train_data3))
tra_data3 <- train_data3[shuffled_index[1:tra_size3], ]
test_data3 <- train_data3[shuffled_index[(tra_size3 + 1):nrow(train_data3)], ]

cat("训练集观测数量:", nrow(tra_data3), "\n")
cat("测试集观测数量:", nrow(test_data3), "\n")
cat("验证集1观测数量:", nrow(val_data31), "\n")
cat("验证集2观测数量:", nrow(val_data32), "\n")

# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data4)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size4 <- round(tra_ratio * nrow(train_data4))
tra_data4 <- train_data4[shuffled_index[1:tra_size4], ]
test_data4 <- train_data4[shuffled_index[(tra_size4 + 1):nrow(train_data4)], ]

cat("训练集观测数量:", nrow(tra_data4), "\n")
cat("测试集观测数量:", nrow(test_data4), "\n")
cat("验证集1观测数量:", nrow(val_data41), "\n")
cat("验证集2观测数量:", nrow(val_data42), "\n")

# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data5)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size5 <- round(tra_ratio * nrow(train_data5))
tra_data5 <- train_data5[shuffled_index[1:tra_size5], ]
test_data5 <- train_data5[shuffled_index[(tra_size5 + 1):nrow(train_data5)], ]

cat("训练集观测数量:", nrow(tra_data5), "\n")
cat("测试集观测数量:", nrow(test_data5), "\n")
cat("验证集1观测数量:", nrow(val_data51), "\n")
cat("验证集2观测数量:", nrow(val_data52), "\n")
# 读取数据
train_data6 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/1.总F编码后_插补.csv")
val_data61 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/2.V1F编码后_插补.csv")
val_data62 <- read.csv("/Users/zj/Desktop/4.机器学习新/5.ptmc/0.数据/3.V2F编码后_插补.csv")
# 将训练数据按照7:3的比例随机分为训练集和验证集
# 设置随机数种子,以确保结果可复现
set.seed(123)
index <- 1:nrow(train_data6)
shuffled_index <- sample(index)
tra_ratio <- 0.7  # 训练集比例
tra_size6 <- round(tra_ratio * nrow(train_data6))
tra_data6 <- train_data6[shuffled_index[1:tra_size6], ]
test_data6 <- train_data6[shuffled_index[(tra_size6 + 1):nrow(train_data6)], ]

cat("训练集观测数量:", nrow(tra_data6), "\n")
cat("测试集观测数量:", nrow(test_data6), "\n")
cat("验证集1观测数量:", nrow(val_data61), "\n")
cat("验证集2观测数量:", nrow(val_data62), "\n")


```

###7.4.2 训练模型

```{r}
# Define the function for training and evaluation
train_and_evaluate <- function(train_data, test_data, val_data1, val_data2, formula) {
  sample_sizes <- seq(100, nrow(train_data), by = 100)
  results <- data.frame()
  
  for (size in sample_sizes) {
    subset_train_data <- train_data[1:size, ]
    
    # Train the model
    model <- glm(formula, data = subset_train_data, family = binomial())
    
    # Predict probabilities
    train_probs <- predict(model, newdata = subset_train_data, type = "response")
    test_probs <- predict(model, newdata = test_data, type = "response")
    val_probs1 <- predict(model, newdata = val_data1, type = "response")
    val_probs2 <- predict(model, newdata = val_data2, type = "response")
    
    # Calculate ROC and AUC
    train_roc <- roc(subset_train_data[[as.character(formula[[2]])]], train_probs)
    test_roc <- roc(test_data[[as.character(formula[[2]])]], test_probs)
    val_roc1 <- roc(val_data1[[as.character(formula[[2]])]], val_probs1)
    val_roc2 <- roc(val_data2[[as.character(formula[[2]])]], val_probs2)
    
    train_auc <- auc(train_roc)
    test_auc <- auc(test_roc)
    val_auc1 <- auc(val_roc1)
    val_auc2 <- auc(val_roc2)
    
    train_ci <- ci.auc(train_roc)
    test_ci <- ci.auc(test_roc)
    val_ci1 <- ci.auc(val_roc1)
    val_ci2 <- ci.auc(val_roc2)
    
    # Store results
    results <- rbind(results, data.frame(
      Sample_Size = size,
      AUC = c(train_auc, test_auc, val_auc1, val_auc2),
      CI_Lower = c(train_ci[1], test_ci[1], val_ci1[1], val_ci2[1]),
      CI_Upper = c(train_ci[3], test_ci[3], val_ci1[3], val_ci2[3]),
      Dataset = factor(c("Train", "Test", "Validation1", "Validation2"), levels = c("Train", "Test", "Validation1", "Validation2"))
    ))
  }
  
  return(results)
}

# Model formulas
formulas <- list(
  Total.Central.Lymph.Node.Metastasis ~ Age + Sex + Tumor.border + Aspect.ratio + Calcification + Tumor.Peripheral.blood.flow + Size + Mulifocality + Extrathyroidal.extension,
  Prelaryngeal.LNM ~ Location + Hashimoto + Pretracheal.LNM + Paratracheal.LNM + LNM.prRLN,
  Pretracheal.LNM ~ Age + Sex + Tumor.Peripheral.blood.flow + Mulifocality + Prelaryngeal.LNM + Paratracheal.LNM + Con.Paratracheal.LNM + LNM.prRLN,
  Paratracheal.LNM ~ Sex + Tumor.border + Aspect.ratio + Size + Extrathyroidal.extension + Prelaryngeal.LNM + Pretracheal.LNM + Con.Paratracheal.LNM + LNM.prRLN,
  Con.Paratracheal.LNM ~ Side.of.position + Pretracheal.LNM + Paratracheal.LNM + LNM.prRLN,
  LNM.prRLN ~ Prelaryngeal.LNM + Pretracheal.LNM + Paratracheal.LNM + Con.Paratracheal.LNM
)

# Calculate AUC for all models
results <- lapply(formulas, function(formula) {
  train_and_evaluate(tra_data1, test_data1, val_data11, val_data12, formula)
})

# Plotting function
plot_learning_curve <- function(data, title, save_path) {
  avg_aucs <- data %>%
    group_by(Dataset) %>%
    summarize(Avg_AUC = mean(AUC))
  
  new_labels <- setNames(paste(levels(data$Dataset), "(mean AUC =", round(avg_aucs$Avg_AUC, 3), ")"), levels(data$Dataset))
  
  p <- ggplot(data, aes(x = Sample_Size, y = AUC, color = Dataset, fill = Dataset)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2) +
    labs(title = title, x = "Training Examples", y = "AUC Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") +
    scale_color_manual(values = c("Train" = "#FF9999",  # Light red
                                  "Test" = "#99FF99",  # Light green
                                  "Validation1" = "#9999FF",  # Light blue
                                  "Validation2" = "#FFFF99"), # Light yellow
                       labels = new_labels, name = "Dataset") +
    scale_fill_manual(values = c("Train" = "#FF9999",  # Light red
                                  "Test" = "#99FF99",  # Light green
                                  "Validation1" = "#9999FF",  # Light blue
                                  "Validation2" = "#FFFF99"), # Light yellow
                      labels = new_labels, name = "Dataset")
  
  # Add AUC value labels on the plot
  p <- p + geom_text(aes(label = round(AUC, 3)), hjust = 1.5, vjust = -0.5, size = 2.5)
  
  # Save plots
  ggsave(filename = paste0(save_path, title, ".tiff"), plot = p, dpi = 300, width = 10, height = 8)
}

# Save plots
save_path <- "/Users/zj/Desktop/4.机器学习新/5.ptmc/2.结果/1.R语言的结果/4.评价2/"

plot_learning_curve(results[[1]], "Learning curve for Total Central Lymph Node Metastasis: Nomogram Prediction", save_path)
plot_learning_curve(results[[2]], "Learning curve for Prelaryngeal Lymph Node Metastasis: Nomogram Prediction", save_path)
plot_learning_curve(results[[3]], "Learning curve for Pretracheal Lymph Node Metastasis: Nomogram Prediction", save_path)
plot_learning_curve(results[[4]], "Learning curve for Paratracheal Lymph Node Metastasis: Nomogram Prediction", save_path)
plot_learning_curve(results[[5]], "Learning curve for Con.Paratracheal Lymph Node Metastasis: Nomogram Prediction", save_path)
plot_learning_curve(results[[6]], "Learning curve for Posterior Recurrent Laryngeal Nerve Lymph Node Metastasis: Nomogram Prediction", save_path)


```


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
