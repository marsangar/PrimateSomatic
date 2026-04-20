# ============================================
# Title: Primate burdens script
# Author: Martín Santamarina García
# Date: 10-02-2026
# Description: This script applies assess burdens
# ============================================

### Load libraries
library(stringr)
library(dplyr)
library(gridExtra)
library("googledrive")
library("googlesheets4")

### Define palette
tissue_cols <- c(
  "Bladder" = "#F8766D",
  "Liver" = "#A3A500",
  "Oesophagus" = "#00BF7D",
  "Prostate" = "#00B0F6",
  "Skin" = "#E76BF3"
)

### Load Metadata
dataMETA<-read_sheet("https://docs.google.com/spreadsheets/d/1k3srj6Kx14Np2vBxxbhEJZQ7s6TdEdo-82_gwwgW0Jk/edit?gid=0#gid=0")
dataMETA$age<-as.numeric(dataMETA$age)
str(dataMETA)

samples<-dataMETA$sampleID
samples_long<-dataMETA$sampleID_long


### Define burden paths
burden_paths<-paste0("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/burdens/", samples_long, ".mut_burden.tsv")

### Read burden files
list_burdens<-lapply(burden_paths, function(x) read.delim(x))
names(list_burdens)<-samples

str(list_burdens)

dataMETA$muts<-sapply(list_burdens, function(x) x["observed","muts"])
dataMETA$burden<-sapply(list_burdens, function(x) x["observed","burden"])

barplot(dataMETA$burdenburdens, main="Mutation burdens", las=2)

### Exploratory
cor.test(dataMETA$muts, dataMETA$burden)
plot(dataMETA$muts, dataMETA$burden)

# Step 1: prepare summary object
summary <- dataMETA %>%
  group_by(sampleID,tissue, age) 

# Step 2: fit linear model
fit <- lm(burden ~ muts, data = summary)

# Step 3: get model stats9
summary_fit <- summary(fit)
intercept <- coef(fit)[1]
slope <- coef(fit)[2]
ci <- confint(fit)  # confidence intervals
cor_val <- cor(summary$muts, summary$burden, method = "pearson")
p_val <- summary_fit$coefficients[2,4]


# 4. Scatter plot, color by tissue, shape by individual
p<-ggplot(summary, aes(x = muts, y = burden, color = tissue, label=sampleID)) +
  geom_point(size = 5) +
  geom_text(size = 3, vjust = -1, col="black") +
  geom_smooth(method = "lm", se = TRUE, inherit.aes = FALSE,
              aes(x = muts, y = burden), color = "black", linetype = "dashed") +
  #coord_cartesian(ylim = c(0, 2200)) +   # Zoom in to 0–0.1 range
  #coord_cartesian(ylim = c(0, 1)) +   # Zoom in to 0–0.1 range
  scale_y_continuous(breaks = seq(0, max(summary$burden, na.rm=TRUE), by = 0.000001))  +
  scale_x_continuous(breaks = seq(0, max(summary$muts, na.rm=TRUE), by = 100))  +
  scale_color_manual(values = tissue_cols) +
  labs(
    title = "Correlation between mutation burden & age",
    subtitle = paste0("Pearson's correlation=", signif(cor_val,2), " // p-val=", signif(p_val,2)),
    x = "Muts",
    y = "Burden",
    color = "Tissue",
    shape = "Sample"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

grid.arrange(p)


