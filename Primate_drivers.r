### Load configs
source("/Users/ms84/Research/postdoc/projects/PrimateSomatic/Primate_master.r")

### Load Mutation Data
dataSOM<-read.delim(paste0(RESDIR_VARIANTS,"/somatic_variants.txt"))

### Load Genes in the targeted panel 
dataPANEL<-read.delim(PANEL)
  
### Evaluate impact of filtering on dnds
dndsout = dndscv(dataSOM, refdb=REFCDS, cv=NULL)
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)

list_dndsout <- lapply(
  split(dataSOM, dataSOM$sample),
  function(df_sample) {
    dndscv(df_sample, refdb = REFCDS, cv = NULL)
  }
)

###
str(list_dndsout)
unique(list_dndsout$MQD0001d$sel_cv$gene_name)


## Select dnds_cv 
df<-dndsout$sel_cv[,c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind")]
dim(df)
### restrict to genes in the panel
df<-df[df$gene_name %in% dataPANEL$gene,]
dim(df)


# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "mutation_type",
    values_to = "count"
  )

# Optional: order genes by total mutation burden
gene_order <- df %>%
  mutate(total = n_syn + n_mis + n_non + n_spl + n_ind) %>%
  arrange(desc(total)) %>%
  pull(gene_name)

df_long$gene_name <- factor(df_long$gene_name, levels = gene_order)

# Rename mutation classes for prettier legend
df_long$mutation_type <- recode(
  df_long$mutation_type,
  n_syn = "Synonymous",
  n_mis = "Missense",
  n_non = "Nonsense",
  n_spl = "Splice",
  n_ind = "Indel"
)

# Stacked barplot
p1<-ggplot(df_long, aes(x = gene_name, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", width = 0.8) +
  #coord_flip() +
  labs(
    title = "Mutation landscape on the Targeted Panel",
    x = "Gene",
    y = "Mutation Count",
    fill = "Mutation Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90, size=5, face = "bold", hjust = 1),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2")

pdf(paste0(RESDIR_DRIVERS,"/drivers.pdf"),  width = 18, height = 6)
p1
dev.off()
