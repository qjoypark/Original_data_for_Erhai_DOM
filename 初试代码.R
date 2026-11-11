packages <- c(
  "tidyverse",
  "readxl",
  "lubridate",
  "janitor",
  "ggplot2",
  "patchwork",
  "lme4",
  "lmerTest",
  "emmeans",
  "broom.mixed",
  "vegan",
  "ggrepel",
  "performance",
  "DHARMa"
)

install_if_missing <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
}

invisible(lapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

file <- "E:/paper/霍达/paper/数据汇总.xlsx"

raw <- read_excel(file)

# 查看结构
str(raw)
dim(raw)
names(raw)

dat <- raw %>%filter(!is.na(Location), !is.na(Month))
dat <- dat %>%mutate(Location = factor(Location),Month = as.integer(Month))
# 检查
dim(dat)
table(dat$Month)
table(dat$Location)
# 每个站点 × 每个月应该只有1个观测
dat %>%count(Location, Month) %>%count(n)

# 2. Create sampling date / seasonal variables
dat <- dat %>%mutate(Year = if_else(Month >= 11, 2023L, 2024L),Date = as.Date(paste(Year,Month,
        "15",sep = "-")),Month_label = factor(month(Date, label = TRUE, abbr = TRUE),
      levels = month(seq(as.Date("2023-11-01"),by = "month",length.out = 12),label = TRUE,abbr = TRUE)))
dat <- dat %>%mutate(Month_label = factor(Month,levels = c(11, 12, 1:10),labels = c(
        "Nov", "Dec", "Jan", "Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug", "Sep", "Oct")))

# ============================================================
# 3. Missing data
vars <- c("WT", "Ph", "TP", "PO4", "TN", "NH4", "NO3","Chla", "DOC", "BIX", "FI", "Fn355","S275-295", "SR", "SUVA254", "M", "a(440)")
missing_summary <- dat %>%summarise(across(all_of(vars),~ sum(is.na(.)))) %>%
  pivot_longer(everything(),names_to = "Variable",values_to = "Missing") %>%
  mutate(Missing_percent = Missing / nrow(dat) * 100) %>%arrange(desc(Missing))
missing_summary
#缺失月份
missing_by_month <- dat %>%group_by(Month_label) %>%summarise(across(all_of(vars),~ sum(is.na(.))),.groups = "drop")
missing_by_month
#站点
missing_by_site <- dat %>%group_by(Location) %>%summarise(across(all_of(vars),~ sum(is.na(.))),.groups = "drop")
missing_by_site



#月均值
# 4. Monthly descriptive statistics
monthly_summary <- dat %>%group_by(Month_label) %>%summarise(
    across(all_of(vars),list(mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),median = ~ median(.x, na.rm = TRUE)),.names = "{.col}_{.fn}"),
    .groups = "drop")
monthly_summary
# 将数据保存为 csv 文件，row.names = FALSE 可避免生成无用的行号列
output_path_csv <- "E:/paper/霍达/paper/monthly_summary.csv"
write.csv(monthly_summary, output_path_csv, row.names = FALSE)

#全部变量的全年热图
core_vars <- c("WT","TN","TP","NH4","NO3","PO4","Chla","DOC")
monthly_core <- dat %>%group_by(Month_label) %>%summarise(across(all_of(core_vars),~ mean(.x, na.rm = TRUE)),
    .groups = "drop")
#转成长格式
monthly_core_long <- monthly_core %>%pivot_longer(cols = -Month_label,names_to = "Variable",values_to = "Value")
#标准化
monthly_core_long <- monthly_core_long %>%group_by(Variable) %>%mutate(Z = as.numeric(scale(Value))) %>%ungroup()
#可视化
p_core_heatmap <- ggplot(monthly_core_long,aes(x = Month_label,y = Variable,fill = Z)) +
  geom_tile(color = "white", linewidth = 0.3) +scale_fill_gradient2(midpoint = 0,name = "Standardized\nvalue") +
  labs(x = NULL,y = NULL) +theme_test(base_size = 12) +theme(
    panel.grid = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
p_core_heatmap

#每个站点的全年那你CHla
ggplot(dat,aes(x = Month_label,y = Chla,group = Location)) +geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +theme_test(base_size = 12) +labs(x = NULL,y = expression(Chl-a))

ggplot(dat,aes(x = Month_label,y = DOC,group = Location))+geom_line(alpha = 0.5) +geom_point(size = 1.5) +
  stat_summary(aes(group = 1),fun = mean,geom = "line",linewidth = 1.2) +
  stat_summary(aes(group = 1),fun = mean,geom = "point",size = 2.5) +theme_test(base_size = 12) +
  labs(x = NULL,y = "DOC")


# 5. Spatial heterogeneity
spatial_stats <- dat %>%group_by(Month_label) %>%
  summarise(across(all_of(core_vars),list(mean = ~ mean(.x, na.rm = TRUE),sd = ~ sd(.x, na.rm = TRUE),cv = ~ sd(.x, na.rm = TRUE) /
          mean(.x, na.rm = TRUE) * 100),.names = "{.col}_{.fn}"),
    .groups = "drop")
spatial_stats

ggplot(spatial_stats,aes(x = Month_label,y = DOC_cv,group = 1)) +
  geom_line(linewidth = 1) +geom_point(size = 2.5) +theme_test(base_size = 12) +
  labs(x = NULL,y = "Spatial CV of DOC (%)")
ggplot(spatial_stats,aes(x = Month_label,y = Chla_cv,group = 1)) +geom_line(linewidth = 1) +
  geom_point(size = 2.5) +theme_test(base_size = 12) +
  labs(x = NULL,y = "Spatial CV of Chl-a (%)")

#PCA
env_vars <- c("WT","TN","TP","NH4","Chla")
env_dat <- dat %>%select(Location,Month_label,all_of(env_vars)) %>%drop_na()
env_pca <- prcomp(env_dat %>%select(all_of(env_vars)),center = TRUE,scale. = TRUE)
summary(env_pca)
env_loadings <- as.data.frame(env_pca$rotation)
env_loadings
round(env_pca$rotation[, 1:2], 3)
#月尺度轨迹
env_scores <- as.data.frame(env_pca$x)
env_scores <- bind_cols(env_dat %>%select(Location, Month_label),env_scores)
ggplot(env_scores,aes(x = PC1,y = PC2,color = Month_label)) +
  geom_point(size = 3,alpha = 0.75) +theme_test(base_size = 12) +
  labs(x = paste0("PC1 (",round(summary(env_pca)$importance[2,1] * 100, 1),"%)"),
    y = paste0("PC2 (",round(summary(env_pca)$importance[2,2] * 100, 1),"%)"),color = "Month")
monthly_pca <- env_scores %>%group_by(Month_label) %>%summarise(PC1 = mean(PC1),PC2 = mean(PC2),.groups = "drop" )
monthly_pca
ggplot(monthly_pca,aes(x = PC1,y = PC2)) +
  geom_path(arrow = arrow(length = unit(0.15, "cm"),type = "closed"),linewidth = 0.8) +
  geom_point(size = 3) +geom_text_repel(aes(label = Month_label),size = 3.5) +
  theme_test(base_size = 12) +labs(x = "PC1",y = "PC2")