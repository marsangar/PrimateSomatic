library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Read data + reshape
df <- read.table("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/counts_dupcaller.tsv", header = TRUE, sep = "\t")

df_long <- df %>%
  pivot_longer(cols = c(ax10, ax50),
               names_to = "method",
               values_to = "count")

# 2. Grouped barplot
ggplot(df_long, aes(x = sample, y = count, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Dupcaller comparison: ax10 vs ax50",
    x = "Sample",
    y = "PASS SNV count"
  ) +
  scale_fill_manual(values = c("ax10" = "#2C7BB6", "ax50" = "#D7191C"))




# 1. Read data + reshape
df <- read.table("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/counts_TDC_TNS.tsv", header = TRUE, sep = "\t")

df_long <- df %>%
  pivot_longer(cols = c(TDC, TNS),
               names_to = "method",
               values_to = "count")

# 2. Grouped barplot
ggplot(df_long, aes(x = sample, y = count, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Targeted Dupcaller vs Targeted Nanoseq",
    x = "Sample",
    y = "PASS SNV count"
  ) +
  scale_fill_manual(values = c("TDC" = "#2C7BB6", "TNS" = "#D7191C"))
