
library(tidyverse)
library(vegan)
library(patchwork)
library(ggrepel)
library(emmeans)
library(grid)
library(scales)


############################################################
# 1. GLOBAL SETTINGS
############################################################

theme_pub <- theme_test(
  base_size = 10.5,
  base_family = "Arial"
) +
  theme(
    axis.title = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    legend.title = element_text(
      size = 9.5
    ),
    
    legend.text = element_text(
      size = 8.5
    ),
    
    strip.text = element_text(
      size = 9.5,
      face = "bold"
    ),
    
    panel.grid = element_blank(),
    
    plot.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      5.5, 7, 5.5, 5.5
    )
  )

theme_set(theme_pub)


############################################################
# 2. MONTH ORDER
############################################################

month_levels <- c(
  "Nov",
  "Dec",
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct"
)


############################################################
# 3. MONTH COLORS
############################################################

month_cols <- c(
  "Nov" = "#4575B4",
  "Dec" = "#4575B4",
  "Jan" = "#4575B4",
  
  "Feb" = "#74ADD1",
  "Mar" = "#74ADD1",
  "Apr" = "#74ADD1",
  
  "May" = "#F6BD60",
  "Jun" = "#F6BD60",
  "Jul" = "#F6BD60",
  
  "Aug" = "#D73027",
  "Sep" = "#D73027",
  "Oct" = "#D73027"
)


############################################################
# 4. DATA PREPARATION
############################################################

dat <- dat %>%
  mutate(
    Month_label = factor(
      Month_label,
      levels = month_levels
    ),
    
    Location = factor(Location),
    
    DIN = NH4 + NO3,
    
    DIP = PO4
  )


############################################################
# 5. OUTPUT DIRECTORY
############################################################

if (!dir.exists("Figures_Final")) {
  dir.create("Figures_Final")
}


fig1_vars <- c(
  "WT",
  "TN",
  "TP",
  "NH4",
  "NO3",
  "PO4",
  "Chla",
  "DOC"
)

env_monthly <- dat %>%
  group_by(Month_label) %>%
  summarise(
    across(
      all_of(fig1_vars),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


############################################################
# 1B. Standardized monthly environmental state
############################################################

env_heat <- env_monthly %>%
  pivot_longer(
    cols = all_of(fig1_vars),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(Variable) %>%
  mutate(
    Z = as.numeric(
      scale(Value)
    )
  ) %>%
  ungroup()


############################################################
# Variable labels
############################################################

env_heat <- env_heat %>%
  mutate(
    Variable = factor(
      Variable,
      levels = rev(fig1_vars),
      labels = rev(
        c(
          "Water temperature",
          "TN",
          "TP",
          "NH4",
          "NO3",
          "PO4",
          "Chl-a",
          "DOC"
        )
      )
    )
  )


############################################################
# 1C. Heatmap
############################################################

p1a <- ggplot(
  env_heat,
  aes(
    x = Month_label,
    y = Variable,
    fill = Z
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +
  
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "Monthly\nz-score"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 8.5
    ),
    
    axis.text.y = element_text(
      size = 8.5
    ),
    
    legend.position = "right"
  )


############################################################
# 1D. Key variables – actual values
############################################################

key_monthly <- env_monthly %>%
  select(
    Month_label,
    WT,
    Chla,
    DOC
  ) %>%
  pivot_longer(
    cols = c(
      WT,
      Chla,
      DOC
    ),
    names_to = "Variable",
    values_to = "Value"
  )


############################################################
# Separate standardized trend plots
############################################################

p1b <- ggplot(
  key_monthly,
  aes(
    x = Month_label,
    y = Value,
    group = Variable
  )
) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.1
  ) +
  
  facet_wrap(
    ~ Variable,
    scales = "free_y",
    ncol = 1
  ) +
  
  labs(
    x = NULL,
    y = NULL
  )


############################################################
# 1E. Figure 1
############################################################

fig1 <- p1a / p1b +
  
  plot_annotation(
    tag_levels = "A"
  )


############################################################
# Save
############################################################

ggsave(
  "Figures_Final/Fig1_Seasonal_environmental_state.pdf",
  fig1,
  width = 8.8,
  height = 8.5
)

ggsave(
  "Figures_Final/Fig1_Seasonal_environmental_state.tiff",
  fig1,
  width = 8.8,
  height = 8.5,
  dpi = 600,
  compression = "lzw"
)

fig1


############################################################
############################################################
# FIGURE 2
# ENVIRONMENTAL PCA
############################################################
############################################################

############################################################
# 2A. Environmental dataset
############################################################

env_vars <- c(
  "WT",
  "TN",
  "TP",
  "NH4",
  "NO3",
  "PO4",
  "Chla"
)

env_complete <- dat %>%
  select(
    Month_label,
    Location,
    all_of(env_vars)
  ) %>%
  drop_na()


############################################################
# PCA
############################################################

env_pca <- prcomp(
  env_complete %>%
    select(all_of(env_vars)),
  center = TRUE,
  scale. = TRUE
)


############################################################
# Scores
############################################################

env_scores <- as.data.frame(
  env_pca$x
) %>%
  mutate(
    Month_label = env_complete$Month_label,
    Location = env_complete$Location
  )


############################################################
# Loadings
############################################################

env_loadings <- as.data.frame(
  env_pca$rotation
) %>%
  rownames_to_column(
    "Variable"
  )


############################################################
# Explained variance
############################################################

env_var <- summary(
  env_pca
)$importance[2, ]


############################################################
# 2B. Sample PCA
############################################################

p2a <- ggplot(
  env_scores,
  aes(
    x = PC1,
    y = PC2,
    colour = Month_label
  )
) +
  
  geom_point(
    size = 1.9,
    alpha = 0.60
  ) +
  
  stat_ellipse(
    aes(
      group = Month_label
    ),
    type = "t",
    linewidth = 0.35,
    linetype = 2,
    alpha = 0.35
  ) +
  
  scale_colour_manual(
    values = month_cols
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(env_var[1] * 100, 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(env_var[2] * 100, 1),
      "%)"
    ),
    
    colour = "Month"
  )


############################################################
# 2C. Loading plot
############################################################

p2b <- ggplot(
  env_loadings,
  aes(
    x = PC1,
    y = PC2
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey60"
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    colour = "grey60"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2
    ),
    arrow = arrow(
      length = unit(
        0.15,
        "cm"
      )
    ),
    linewidth = 0.5
  ) +
  
  geom_text_repel(
    aes(
      label = Variable
    ),
    size = 3.1,
    max.overlaps = Inf
  ) +
  
  coord_equal() +
  
  labs(
    x = paste0(
      "PC1 (",
      round(env_var[1] * 100, 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(env_var[2] * 100, 1),
      "%)"
    )
  )


############################################################
# 2D. Monthly centroids
############################################################

env_centroid <- env_scores %>%
  group_by(Month_label) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  ) %>%
  mutate(
    Month_order = match(
      Month_label,
      month_levels
    )
  ) %>%
  arrange(Month_order)


############################################################
# Monthly trajectory
############################################################

p2c <- ggplot(
  env_centroid,
  aes(
    x = PC1,
    y = PC2
  )
) +
  
  geom_path(
    linewidth = 0.65,
    arrow = arrow(
      length = unit(
        0.14,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_point(
    aes(
      fill = Month_label
    ),
    shape = 21,
    size = 3.2,
    colour = "black",
    stroke = 0.35
  ) +
  
  geom_text_repel(
    aes(
      label = Month_label
    ),
    size = 3.0,
    box.padding = 0.30
  ) +
  
  scale_fill_manual(
    values = month_cols
  ) +
  
  labs(
    x = paste0(
      "Environmental PC1 (",
      round(env_var[1] * 100, 1),
      "%)"
    ),
    
    y = paste0(
      "Environmental PC2 (",
      round(env_var[2] * 100, 1),
      "%)"
    ),
    
    fill = "Month"
  )


############################################################
# 2E. Figure 2
############################################################

fig2 <- (
  
  p2a | p2b
  
) /
  
  p2c +
  
  plot_annotation(
    tag_levels = "A"
  )


############################################################
# Save
############################################################

ggsave(
  "Figures_Final/Fig2_Environmental_PCA.pdf",
  fig2,
  width = 9.2,
  height = 8.8
)

ggsave(
  "Figures_Final/Fig2_Environmental_PCA.tiff",
  fig2,
  width = 9.2,
  height = 8.8,
  dpi = 600,
  compression = "lzw"
)

fig2

############################################################
# FIGURE 3
# Seasonal dynamics of Chl-a and its associations
# with nutrient availability
############################################################

library(tidyverse)
library(broom)
library(patchwork)
library(scales)

############################################################
# 3.0 Data preparation
############################################################

chla_dat <- dat %>%
  mutate(
    DIN = NH4 + NO3,
    DIP = PO4
  ) %>%
  select(
    Chla,
    Month_label,
    WT,
    DIN,
    DIP,
    Location
  ) %>%
  drop_na() %>%
  mutate(
    Month_label = factor(
      Month_label,
      levels = c(
        "Nov", "Dec", "Jan", "Feb", "Mar", "Apr",
        "May", "Jun", "Jul", "Aug", "Sep", "Oct"
      )
    ),
    Location = factor(Location)
  )

############################################################
# Final Chl-a model
############################################################

m_chla <- lm(
  log1p(Chla) ~ Month_label + WT + DIN + DIP + Location,
  data = chla_dat
)

summary(m_chla)

############################################################
# 3A. Observed monthly Chl-a
############################################################

chla_month <- chla_dat %>%
  group_by(Month_label) %>%
  summarise(
    mean_Chla = mean(Chla, na.rm = TRUE),
    sd_Chla = sd(Chla, na.rm = TRUE),
    n = n(),
    se_Chla = sd_Chla / sqrt(n),
    .groups = "drop"
  )

p3A <- ggplot(
  chla_dat,
  aes(x = Month_label, y = Chla)
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = "grey92",
    colour = "grey35",
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.4,
    alpha = 0.55
  ) +
  geom_point(
    data = chla_month,
    aes(y = mean_Chla),
    size = 2.5
  ) +
  labs(
    x = NULL,
    y = expression("Chl-a ("*mu*"g L"^-1*")")
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title.y = element_text(
      margin = margin(r = 8)
    )
  )

############################################################
# Function for model-adjusted partial effects
############################################################

partial_effect <- function(
    model,
    data,
    focal,
    n = 100
) {
  
  # coefficient
  beta <- coef(model)[focal]
  
  # variance of focal coefficient
  beta_se <- sqrt(
    vcov(model)[focal, focal]
  )
  
  # range of focal variable
  x_seq <- seq(
    min(data[[focal]], na.rm = TRUE),
    max(data[[focal]], na.rm = TRUE),
    length.out = n
  )
  
  # centered around mean
  x_mean <- mean(
    data[[focal]],
    na.rm = TRUE
  )
  
  # fitted contribution on model scale
  fit <- beta * (x_seq - x_mean)
  
  # 95% CI
  crit <- qt(
    0.975,
    df = df.residual(model)
  )
  
  lower <- fit - crit * beta_se * abs(x_seq - x_mean)
  upper <- fit + crit * beta_se * abs(x_seq - x_mean)
  
  tibble(
    x = x_seq,
    fit = fit,
    lower = lower,
    upper = upper,
    variable = focal
  )
}

############################################################
# 3B–3D. Partial effects
############################################################

effect_DIN <- partial_effect(
  model = m_chla,
  data = chla_dat,
  focal = "DIN"
)

effect_DIP <- partial_effect(
  model = m_chla,
  data = chla_dat,
  focal = "DIP"
)

effect_WT <- partial_effect(
  model = m_chla,
  data = chla_dat,
  focal = "WT"
)
############################################################
# Model coefficients for annotation
############################################################

coef_tab <- broom::tidy(m_chla)

beta_DIN <- coef_tab %>%
  filter(term == "DIN") %>%
  pull(estimate)

p_DIN <- coef_tab %>%
  filter(term == "DIN") %>%
  pull(p.value)

beta_DIP <- coef_tab %>%
  filter(term == "DIP") %>%
  pull(estimate)

p_DIP <- coef_tab %>%
  filter(term == "DIP") %>%
  pull(p.value)

beta_WT <- coef_tab %>%
  filter(term == "WT") %>%
  pull(estimate)

p_WT <- coef_tab %>%
  filter(term == "WT") %>%
  pull(p.value)

############################################################
# 3B. Partial effect of DIN
############################################################

p3B <- ggplot(
  effect_DIN,
  aes(x = x, y = fit)
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey55"
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "β = %.3f\np = %.3f",
      beta_DIN,
      p_DIN
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 3.2
  ) +
  labs(
    x = expression("DIN ("*mu*"mol L"^-1*")"),
    y = expression("Partial effect on log1p(Chl-a)")
  ) +
  theme_classic(base_size = 10.5)
############################################################
# 3C. Partial effect of DIP
############################################################

p3C <- ggplot(
  effect_DIP,
  aes(x = x, y = fit)
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey55"
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "β = %.3f\np = %.3f",
      beta_DIP,
      p_DIP
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 3.2
  ) +
  labs(
    x = expression("DIP ("*mu*"mol L"^-1*")"),
    y = expression("Partial effect on log1p(Chl-a)")
  ) +
  theme_classic(base_size = 10.5)
############################################################
# 3D. Partial effect of WT
############################################################

p3D <- ggplot(
  effect_WT,
  aes(x = x, y = fit)
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey55"
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "β = %.3f\np = %.3f",
      beta_WT,
      p_WT
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 3.2
  ) +
  labs(
    x = "Water temperature (°C)",
    y = expression("Partial effect on log1p(Chl-a)")
  ) +
  theme_classic(base_size = 10.5)
############################################################
# Combine Fig. 3
############################################################

p3 <- p3A /
  (p3B | p3C | p3D) +
  plot_layout(
    heights = c(1.15, 1)
  ) &
  theme(
    plot.margin = margin(5, 8, 5, 8)
  )

p3


############################################################
############################################################
# FIGURE 4
# LATE-SEASON DOC REGIME
############################################################
############################################################

############################################################
# 4A. Monthly DOC
############################################################

doc_month <- dat %>%
  group_by(Month_label) %>%
  summarise(
    mean_DOC = mean(DOC, na.rm = TRUE),
    se_DOC = sd(DOC, na.rm = TRUE) /
      sqrt(sum(!is.na(DOC))),
    .groups = "drop"
  )


############################################################
# 4B. Late-season definition
############################################################

doc_dat <- dat %>%
  filter(
    complete.cases(
      DOC,
      Chla,
      WT,
      TN,
      TP,
      Month_label,
      Location
    )
  ) %>%
  mutate(
    LateSeason = if_else(
      Month_label %in%
        c("Aug", "Sep", "Oct"),
      "Aug–Oct",
      "Other months"
    ),
    
    LateSeason = factor(
      LateSeason,
      levels = c(
        "Other months",
        "Aug–Oct"
      )
    )
  )


############################################################
# 4C. DOC model
############################################################

m_doc <- lm(
  log1p(DOC) ~
    Month_label +
    log1p(Chla) +
    WT +
    TN +
    TP +
    Location,
  data = doc_dat
)

m_doc_late <- lm(
  log1p(DOC) ~
    LateSeason +
    Location,
  data = doc_dat
)


############################################################
# 4D. Monthly DOC trend
############################################################

p4a <- ggplot(
  doc_month,
  aes(
    x = Month_label,
    y = mean_DOC,
    group = 1
  )
) +
  
  geom_ribbon(
    aes(
      ymin = mean_DOC - se_DOC,
      ymax = mean_DOC + se_DOC
    ),
    alpha = 0.15
  ) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  annotate(
    "rect",
    xmin = 9.5,
    xmax = 12.5,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.06
  ) +
  
  annotate(
    "text",
    x = 10.5,
    y = max(doc_month$mean_DOC) * 0.95,
    label = "Late-season\nhigh-DOC state",
    size = 3.2,
    hjust = 0.5
  ) +
  
  labs(
    x = NULL,
    y = "DOC (mg L⁻¹)"
  )


############################################################
# 4E. Late-season effect
############################################################

late_summary <- doc_dat %>%
  group_by(LateSeason) %>%
  summarise(
    mean = mean(log1p(DOC)),
    se = sd(log1p(DOC)) / sqrt(n()),
    .groups = "drop"
  )


p4b <- ggplot(
  late_summary,
  aes(
    x = LateSeason,
    y = mean
  )
) +
  
  geom_col(
    width = 0.58,
    fill = "grey75",
    colour = "black"
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.12,
    linewidth = 0.5
  ) +
  
  labs(
    x = NULL,
    y = expression(
      Mean~log[ e ](1 + DOC)
    )
  )


############################################################
# 4F. DOC vs Chl-a
############################################################

p4c <- ggplot(
  doc_dat,
  aes(
    x = log1p(Chla),
    y = log1p(DOC)
  )
) +
  
  geom_point(
    aes(
      fill = LateSeason
    ),
    shape = 21,
    size = 2.4,
    alpha = 0.70,
    colour = "black",
    stroke = 0.25
  ) +
  
  scale_fill_manual(
    values = c(
      "Other months" = "grey75",
      "Aug–Oct" = "#D73027"
    )
  ) +
  
  labs(
    x = expression(
      log[ e ](1 + Chl-a)
    ),
    y = expression(
      log[ e ](1 + DOC)
    ),
    fill = NULL
  )+theme(legend.position = "none")


############################################################
# 4G. DOC vs WT
############################################################

p4d <- ggplot(
  doc_dat,
  aes(
    x = WT,
    y = log1p(DOC)
  )
) +
  
  geom_point(
    aes(
      fill = LateSeason
    ),
    shape = 21,
    size = 2.4,
    alpha = 0.70,
    colour = "black",
    stroke = 0.25
  ) +
  
  scale_fill_manual(
    values = c(
      "Other months" = "grey75",
      "Aug–Oct" = "#D73027"
    )
  ) +
  
  labs(
    x = "Water temperature (°C)",
    y = expression(
      log[ e ](1 + DOC)
    ),
    fill = NULL
  )


############################################################
# 4H. Figure 4
############################################################

fig4 <- (
  
  p4a | p4b
  
) /
  
  (
    
    p4c | p4d
    
  ) +
  
  plot_annotation(
    tag_levels = "A"
  )


############################################################
# Save
############################################################

ggsave(
  "Figures_Final/Fig4_Late_season_DOC.pdf",
  fig4,
  width = 9.2,
  height = 8.0
)

ggsave(
  "Figures_Final/Fig4_Late_season_DOC.tiff",
  fig4,
  width = 9.2,
  height = 8.0,
  dpi = 600,
  compression = "lzw"
)

fig4


############################################################
############################################################
# FIGURE 5
# DOM COMPOSITIONAL REORGANIZATION
############################################################
############################################################

############################################################
# 5A. DOM variables
############################################################

dom_vars <- c(
  "DOC",
  "BIX",
  "FI",
  "Fn355",
  "S275-295",
  "SR",
  "SUVA254",
  "M",
  "a.440."
)


############################################################
# 5B. DOM PCA
############################################################

############################################################
# 5B. DOM PCA
############################################################

dom_vars <- c(
  "DOC",
  "BIX",
  "FI",
  "Fn355",
  "S275-295",
  "SR",
  "SUVA254",
  "M",
  "a(440)"
)

dom_complete <- dat %>%
  select(
    Month_label,
    Location,
    all_of(dom_vars)
  ) %>%
  drop_na()

# DOM PCA
dom_pca <- prcomp(
  dom_complete %>%
    select(all_of(dom_vars)),
  center = TRUE,
  scale. = TRUE
)

# PCA scores
dom_scores <- dom_complete %>%
  mutate(
    PC1 = dom_pca$x[, 1],
    PC2 = dom_pca$x[, 2],
    PC3 = dom_pca$x[, 3]
  )


#########################################################
# Loadings
############################################################

dom_loadings <- as.data.frame(
  dom_pca$rotation
) %>%
  rownames_to_column(
    "Variable"
  )


############################################################
# Explained variance
############################################################

dom_var <- summary(
  dom_pca
)$importance[2, ]


############################################################
# 5C. DOM sample space
############################################################

p5a <- ggplot(
  dom_scores,
  aes(
    x = PC1,
    y = PC2,
    colour = Month_label
  )
) +
  
  geom_point(
    size = 1.9,
    alpha = 0.60
  ) +
  
  scale_colour_manual(
    values = month_cols
  ) +
  
  labs(
    x = paste0(
      "DOM PC1 (",
      round(dom_var[1] * 100, 1),
      "%)"
    ),
    
    y = paste0(
      "DOM PC2 (",
      round(dom_var[2] * 100, 1),
      "%)"
    ),
    
    colour = "Month"
  )


############################################################
# 5D. DOM loadings
############################################################

p5b <- ggplot(
  dom_loadings,
  aes(
    x = PC1,
    y = PC2
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey60"
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    colour = "grey60"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2
    ),
    arrow = arrow(
      length = unit(
        0.15,
        "cm"
      )
    ),
    linewidth = 0.5
  ) +
  
  geom_text_repel(
    aes(
      label = Variable
    ),
    size = 3.0,
    max.overlaps = Inf
  ) +
  
  coord_equal() +
  
  labs(
    x = paste0(
      "PC1 (",
      round(dom_var[1] * 100, 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(dom_var[2] * 100, 1),
      "%)"
    )
  )


############################################################
# 5E. Monthly DOM centroids
############################################################

dom_centroid <- dom_scores %>%
  group_by(Month_label) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    PC3 = mean(PC3),
    .groups = "drop"
  ) %>%
  mutate(
    Month_order = match(
      Month_label,
      month_levels
    )
  ) %>%
  arrange(Month_order)


############################################################
# 5F. Monthly DOM trajectory
############################################################

p5c <- ggplot(
  dom_centroid,
  aes(
    x = PC1,
    y = PC2
  )
) +
  
  geom_path(
    linewidth = 0.65,
    arrow = arrow(
      length = unit(
        0.14,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_point(
    aes(
      fill = Month_label
    ),
    shape = 21,
    size = 3.2,
    colour = "black",
    stroke = 0.35
  ) +
  
  geom_text_repel(
    aes(
      label = Month_label
    ),
    size = 3.0,
    box.padding = 0.30
  ) +
  
  scale_fill_manual(
    values = month_cols
  ) +
  
  labs(
    x = paste0(
      "DOM PC1 (",
      round(dom_var[1] * 100, 1),
      "%)"
    ),
    
    y = paste0(
      "DOM PC2 (",
      round(dom_var[2] * 100, 1),
      "%)"
    ),
    
    fill = "Month"
  )


############################################################
# 5G. Selected optical indices
############################################################

dom_optical <- dat %>%
  group_by(Month_label) %>%
  summarise(
    BIX = mean(BIX, na.rm = TRUE),
    FI = mean(FI, na.rm = TRUE),
    SUVA254 = mean(SUVA254, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      BIX,
      FI,
      SUVA254
    ),
    names_to = "Indicator",
    values_to = "Value"
  )


p5d <- ggplot(
  dom_optical,
  aes(
    x = Month_label,
    y = Value,
    group = Indicator
  )
) +
  
  geom_line(
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 1.8
  ) +
  
  facet_wrap(
    ~ Indicator,
    scales = "free_y",
    ncol = 1
  ) +
  
  labs(
    x = NULL,
    y = NULL
  )


############################################################
# 5H. Figure 5
############################################################

fig5 <- (
  
  p5a | p5b
  
) /
  
  (
    
    p5c | p5d
    
  ) +
  
  plot_annotation(
    tag_levels = "A"
  )


############################################################
# Save
############################################################

ggsave(
  "Figures_Final/Fig5_DOM_composition.pdf",
  fig5,
  width = 9.2,
  height = 9.0
)

ggsave(
  "Figures_Final/Fig5_DOM_composition.tiff",
  fig5,
  width = 9.2,
  height = 9.0,
  dpi = 600,
  compression = "lzw"
)

fig5


############################################################
############################################################
# FIGURE 6
# ENVIRONMENTAL–DOM COUPLING
############################################################
############################################################

############################################################
# 6A. Environmental monthly centroids
############################################################

env_centroid6 <- env_scores %>%
  group_by(Month_label) %>%
  summarise(
    Env_PC1 = mean(PC1),
    Env_PC2 = mean(PC2),
    .groups = "drop"
  )


############################################################
# 6B. DOM monthly centroids
############################################################

dom_centroid6 <- dom_scores %>%
  group_by(Month_label) %>%
  summarise(
    DOM_PC1 = mean(PC1),
    DOM_PC2 = mean(PC2),
    DOM_PC3 = mean(PC3),
    .groups = "drop"
  )


############################################################
# 6C. Combined monthly state
############################################################

monthly_state <- env_centroid6 %>%
  left_join(
    dom_centroid6,
    by = "Month_label"
  ) %>%
  mutate(
    Month_order = match(
      Month_label,
      month_levels
    )
  ) %>%
  arrange(Month_order)


############################################################
# 6D. Transition distances
############################################################

transition_dat <- monthly_state %>%
  mutate(
    
    Env_transition = sqrt(
      (Env_PC1 - lag(Env_PC1))^2 +
        (Env_PC2 - lag(Env_PC2))^2
    ),
    
    DOM_transition = sqrt(
      (DOM_PC1 - lag(DOM_PC1))^2 +
        (DOM_PC2 - lag(DOM_PC2))^2 +
        (DOM_PC3 - lag(DOM_PC3))^2
    ),
    
    Transition = paste(
      lag(as.character(Month_label)),
      as.character(Month_label),
      sep = "–"
    )
    
  ) %>%
  filter(
    !is.na(Env_transition),
    !is.na(DOM_transition)
  )


############################################################
# 6E. Mantel test
############################################################

env_dist_month <- dist(
  monthly_state %>%
    select(
      Env_PC1,
      Env_PC2
    )
)


dom_dist_month <- dist(
  monthly_state %>%
    select(
      DOM_PC1,
      DOM_PC2,
      DOM_PC3
    )
)


mantel_result <- mantel(
  env_dist_month,
  dom_dist_month,
  method = "spearman",
  permutations = 999
)


############################################################
# 6F. Consecutive-transition Spearman test
############################################################

transition_cor <- cor.test(
  transition_dat$Env_transition,
  transition_dat$DOM_transition,
  method = "spearman"
)


############################################################
# 6G. State-space plot
############################################################

p6a <- ggplot(
  monthly_state,
  aes(
    x = Env_PC1,
    y = DOM_PC1
  )
) +
  
  geom_path(
    linewidth = 0.65,
    arrow = arrow(
      length = unit(
        0.14,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_point(
    aes(
      fill = Month_label
    ),
    shape = 21,
    size = 3.3,
    colour = "black",
    stroke = 0.35
  ) +
  
  geom_text_repel(
    aes(
      label = Month_label
    ),
    size = 3.0,
    box.padding = 0.3
  ) +
  
  scale_fill_manual(
    values = month_cols
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Mantel r = ",
      round(
        mantel_result$statistic,
        3
      ),
      "\np = ",
      format.pval(
        mantel_result$signif,
        digits = 2,
        eps = 0.001
      )
    ),
    hjust = 1.05,
    vjust = 1.25,
    size = 3.2
  ) +
  
  labs(
    x = "Environmental PC1",
    y = "DOM PC1",
    fill = "Month"
  )




############################################################
# 6H. Transition distance
############################################################

transition_long <- transition_dat %>%
  select(
    Transition,
    Env_transition,
    DOM_transition
  ) %>%
  pivot_longer(
    cols = c(
      Env_transition,
      DOM_transition
    ),
    names_to = "Type",
    values_to = "Distance"
  ) %>%
  mutate(
    Type = dplyr::recode(
      Type,
      "Env_transition" = "Environmental",
      "DOM_transition" = "DOM"
    )
  )


p6b <- ggplot(
  transition_long,
  aes(
    x = Transition,
    y = Distance,
    group = Type,
    linetype = Type
  )
) +
  
  geom_line(
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 1.9
  ) +
  
  labs(
    x = "Consecutive monthly transition",
    y = "Transition distance",
    linetype = NULL
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.5
    )
  )


############################################################
# 6I. Transition distance relationship
############################################################

p6c <- ggplot(
  transition_dat,
  aes(
    x = Env_transition,
    y = DOM_transition
  )
) +
  
  geom_point(
    size = 2.7
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.65
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Spearman \u03C1 = ",
      round(
        transition_cor$estimate,
        3
      ),
      "\np = ",
      format.pval(
        transition_cor$p.value,
        digits = 2,
        eps = 0.001
      )
    ),
    hjust = 1.05,
    vjust = 1.25,
    size = 3.2
  ) +
  
  labs(
    x = "Environmental transition distance",
    y = "DOM transition distance"
  )


############################################################
# 6J. Figure 6
############################################################

fig6 <- (
  
  p6a | p6b
  
) /
  
  p6c +
  
  plot_annotation(
    tag_levels = "A"
  )


############################################################
# Save
############################################################

ggsave(
  "Figures_Final/Fig6_Environmental_DOM_coupling.pdf",
  fig6,
  width = 9.2,
  height = 8.0
)

ggsave(
  "Figures_Final/Fig6_Environmental_DOM_coupling.tiff",
  fig6,
  width = 9.2,
  height = 8.0,
  dpi = 600,
  compression = "lzw"
)

fig6


############################################################
############################################################
# 7. SAVE ANALYSIS OBJECTS
############################################################
############################################################

save(
  env_pca,
  env_scores,
  env_loadings,
  env_centroid,
  dom_pca,
  dom_scores,
  dom_loadings,
  dom_centroid,
  monthly_state,
  transition_dat,
  mantel_result,
  transition_cor,
  m_chla,
  m_doc,
  m_doc_late,
  file = "Figures_Final/Figure_analysis_objects.RData"
)


############################################################
# 8. EXPORT STATISTICS FOR MANUSCRIPT
############################################################

cat(
  "\n============================================\n",
  "ENVIRONMENTAL PCA\n",
  "============================================\n"
)

print(
  summary(env_pca)
)


cat(
  "\n============================================\n",
  "DOM PCA\n",
  "============================================\n"
)

print(
  summary(dom_pca)
)


cat(
  "\n============================================\n",
  "MANTEL TEST\n",
  "============================================\n"
)

print(
  mantel_result
)


cat(
  "\n============================================\n",
  "TRANSITION CORRELATION\n",
  "============================================\n"
)

print(
  transition_cor
)


cat(
  "\n============================================\n",
  "CHL-a MODEL\n",
  "============================================\n"
)

print(
  summary(m_chla)
)


cat(
  "\n============================================\n",
  "DOC MODEL\n",
  "============================================\n"
)

print(
  summary(m_doc)
)


cat(
  "\n============================================\n",
  "LATE-SEASON DOC MODEL\n",
  "============================================\n"
)

print(
  summary(m_doc_late)
)


cat(
  "\n\nAll publication figures saved in:\n",
  "Figures_Final/\n"
)