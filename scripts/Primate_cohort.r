# ============================================
# PrimateSomatic — Cohort summary plots
# --------------------------------------------
# Author: Martín Santamarina García
# Contact: ms3242@cam.ac.uk
# Date: 29/09/2026

### Script co-generated with Claude by Anthropic

#### This script is intended to draw summary plots for a primate cohort, based on the information collected as part of the cross-species-master-database.

# ==============================================================================
# Rhesus macaque cohort summary figure (tailored to the master sample sheet)
#  A  Samples per tissue, coloured by furthest processing stage reached
#  B  Sex balance per tissue
#  C  Tissue x age-group coverage (sample counts)
#  D  Age distribution of donor animals, by sex
#  E  Overall workflow progress (% of samples completing each step)
# ==============================================================================

# install.packages(c("googlesheets4", "tidyverse", "janitor", "patchwork",
#                    "scales", "glue", "ragg"))
library(googlesheets4)
library(tidyverse)
library(janitor)
library(patchwork)
library(scales)
library(glue)

# ---- 1. CONFIG ---------------------------------------------------------------
sheet_url  <- "https://docs.google.com/spreadsheets/d/1iwvh9JNUOJ4p2KDC5XDlOtklXw8zcjtLVi__Jp3gVYM/edit?gid=276995000#gid=276995000"
sheet_tab  <- 1
taxon      <- "9544"                 # Rhesus macaque; keeps only these rows
out_prefix <- "macaque_cohort_summary"

# Age bins in years (adjust freely)
age_breaks <- c(0, 1, 3, 5, 10, 15, Inf)
age_labels <- c("<1", "1-3", "3-5", "5-10", "10-15", "15+")

# ---- 2. LOAD -----------------------------------------------------------------
gs4_auth()                           # gs4_deauth() if the sheet is link-shared/public

# Read everything as text: mixed TRUE/FALSE/date/blank columns otherwise become
# list-columns or get silently coerced.
raw <- read_sheet(sheet_url, sheet = sheet_tab, col_types = "c") |>
  clean_names()                      # "Paxgene Fixation?" -> paxgene_fixation, etc.

# TRUE for TRUE/"yes"/any date or value; FALSE for blank/FALSE/NO/0
is_done <- function(x) {
  x <- str_to_upper(str_squish(as.character(x)))
  !is.na(x) & !x %in% c("", "FALSE", "F", "NO", "N", "0", "NA")
}

# ---- 3. CLEAN ----------------------------------------------------------------
df <- raw |>
  filter(taxon_id == taxon, !is.na(tissue)) |>
  transmute(
    sample = supplier_id,
    # PD_ID = animal ID + tissue letter (MQX0008l, MQX0008b) -> strip last letter
    animal = coalesce(str_remove(pd_id, "[a-zA-Z]$"),
                      str_remove(supplier_id, "_.*$")),
    tissue = str_to_title(str_squish(tissue)),
    age    = coalesce(as.numeric(age_at_collection_days) / 365.25,
                      as.numeric(age_at_collection_months) / 12,
                      as.numeric(age_at_collection_years)),
    sex    = factor(str_to_title(str_squish(sex)), levels = c("Female", "Male")),
    fixed     = is_done(paxgene_fixation),
    histology = is_done(histology_request_complete),
    lcm       = is_done(lcm_performed),
    sorted    = is_done(sorted),
    extracted = is_done(submitted_for_extraction),
    pooled    = is_done(pooled_for_library_preparation),
    sequenced = is_done(submitted_for_sequencing)
  ) |>
  mutate(
    age_group = cut(age, breaks = age_breaks, labels = age_labels,
                    right = FALSE, include.lowest = TRUE),
    # Furthest stage reached (most advanced wins)
    stage = case_when(
      sequenced ~ "Sequencing submitted",
      pooled    ~ "Library pooled",
      extracted ~ "DNA extraction",
      histology ~ "Histology complete",
      fixed     ~ "PAXgene fixed",
      TRUE      ~ "Received"
    ),
    stage = factor(stage, levels = c("Received", "PAXgene fixed",
                                     "Histology complete", "DNA extraction",
                                     "Library pooled", "Sequencing submitted"))
  ) |>
  distinct(sample, .keep_all = TRUE)

message("Samples: ", nrow(df), " | animals: ", n_distinct(df$animal),
        " | tissues: ", n_distinct(df$tissue),
        " | missing age: ", sum(is.na(df$age)),
        " | missing sex: ", sum(is.na(df$sex)))

tissue_order <- df |> count(tissue, sort = TRUE) |> pull(tissue)
df <- df |> mutate(tissue = factor(tissue, levels = rev(tissue_order)))  # biggest on top
animals <- df |> distinct(animal, .keep_all = TRUE)

# ---- 4. THEME & PALETTES -----------------------------------------------------
sex_cols <- c(Female = "#CC79A7", Male = "#0072B2")   # Okabe-Ito, colourblind-safe

theme_set(
  theme_classic(base_size = 11) +
    theme(
      axis.text        = element_text(colour = "black"),
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.4),
      legend.key.size  = unit(0.4, "cm"),
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 9),
      plot.title       = element_text(face = "bold", size = 11),
      plot.tag         = element_text(face = "bold", size = 16)
    )
)

# ---- 5. PANELS ---------------------------------------------------------------
# A: tissue x furthest stage
pA <- df |>
  count(tissue, stage) |>
  ggplot(aes(n, tissue, fill = fct_rev(stage))) +
  geom_col(width = 0.75) +
  geom_text(data = count(df, tissue), aes(n, tissue, label = n),
            inherit.aes = FALSE, hjust = -0.3, size = 3.2) +
  scale_fill_viridis_d(option = "mako", begin = 0.15, end = 0.9, direction = -1,
                       breaks = levels(df$stage), labels = levels(df$stage),
                       name = "Furthest stage", guide = guide_legend(reverse = FALSE)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Tissues & processing stage", x = "Samples (n)", y = NULL)

# B: sex balance per tissue
pB <- df |>
  filter(!is.na(sex)) |>
  count(tissue, sex, .drop = FALSE) |>
  group_by(tissue) |> mutate(prop = n / sum(n)) |> ungroup() |>
  ggplot(aes(prop, tissue, fill = sex)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = sex_cols, name = "Sex", drop = FALSE) +
  scale_x_continuous(labels = percent, breaks = c(0, .5, 1), expand = c(0, 0)) +
  labs(title = "Sex balance", x = "Samples (%)", y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# C: tissue x age group heatmap
pC <- df |>
  filter(!is.na(age_group)) |>
  count(tissue, age_group, .drop = FALSE) |>
  ggplot(aes(age_group, tissue, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(n > 0, n, "")), size = 3,
            colour = "white", fontface = "bold") +
  scale_fill_gradient(low = "#D6E4F0", high = "#08306B", name = "Samples",
                      breaks = pretty_breaks(3), limits = c(0, NA)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Age coverage", x = "Age at collection (years)", y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank())

# D: age distribution of donors by sex
pD <- animals |>
  filter(!is.na(age), !is.na(sex)) |>
  ggplot(aes(age, fill = sex)) +
  geom_histogram(binwidth = 1, boundary = 0, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = sex_cols, name = "Sex", drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(breaks = pretty_breaks(), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Donor animals: age & sex", x = "Age at collection (years)",
       y = "Animals (n)")

# E: overall workflow progress
steps <- c(fixed = "PAXgene fixed", histology = "Histology complete",
           lcm = "LCM performed", sorted = "Cells sorted",
           extracted = "DNA extraction", pooled = "Library pooled",
           sequenced = "Sequencing submitted")

pE <- df |>
  summarise(across(all_of(names(steps)), sum)) |>
  pivot_longer(everything(), names_to = "step", values_to = "done") |>
  mutate(step  = factor(steps[step], levels = rev(steps)),
         total = nrow(df), frac = done / total) |>
  ggplot(aes(y = step)) +
  geom_col(aes(x = 1), fill = "grey90", width = 0.7) +
  geom_col(aes(x = frac), fill = "#009E73", width = 0.7) +
  geom_text(aes(x = 1.02, label = glue("{done}/{total}")), hjust = 0, size = 3.2) +
  scale_x_continuous(labels = percent, breaks = c(0, .5, 1),
                     expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Workflow progress", x = "Samples completed (%)", y = NULL)

# ---- 6. ASSEMBLE -------------------------------------------------------------
n_samples <- nrow(df)
n_animals <- n_distinct(df$animal)
n_tissues <- n_distinct(df$tissue)
age_rng   <- range(animals$age, na.rm = TRUE)
n_seq     <- sum(df$sequenced)

#fig <- pA + pB + pC + pD + pE +
fig <- pA + pB + pC +  pE +
  plot_layout(design = "ABC\nDDE", widths = c(1.5, 0.8, 1.3), heights = c(1.5, 1),
              guides = "collect") +
  plot_annotation(
    title    = expression(italic("Macaca mulatta") ~ "cohort overview"),
    subtitle = glue("{n_samples} samples  |  {n_animals} animals  |  ",
                    "{n_tissues} tissues  |  ages {round(age_rng[1], 1)}-",
                    "{round(age_rng[2], 1)} years  |  {n_seq} submitted for sequencing"),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11, colour = "grey30"))
  ) &
  theme(legend.position = "right")

fig

# ---- 7. EXPORT ---------------------------------------------------------------
w <- 12; h <- 8
ggsave(paste0(out_prefix, ".pdf"),  fig, width = w, height = h, device = cairo_pdf)
ggsave(paste0(out_prefix, ".png"),  fig, width = w, height = h, dpi = 300,
       device = ragg::agg_png)
ggsave(paste0(out_prefix, ".tiff"), fig, width = w, height = h, dpi = 300,
       device = ragg::agg_tiff, compression = "lzw")





