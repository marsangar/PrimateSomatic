# Genomic distribution of a gene panel on vertical ideograms (hg38)
#
# Install once:
#   install.packages(c("BiocManager", "ggplot2", "ggrepel"))
#   BiocManager::install(c("karyoploteR", "org.Hs.eg.db",
#                          "TxDb.Hsapiens.UCSC.hg38.knownGene"))

suppressPackageStartupMessages({
  library(karyoploteR)   # only for getCytobands()
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
})

#### Currently this code provides display for hg38

# ---- Settings ----------------------------------------------------------------
gene_file <- "/Users/ms84/Downloads/Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv"
out_file  <- "/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/genes/genes_karyotype_vertical"

highlight_genes <- c("NOTCH1", "FAT1", "TP53", "NOTCH2", "ATM", "MGA", "ASXL1", "CHEK2", "TET2", "PPM1D", "BCORL1", "ARID1A", "KDM5C", "RAC1", "EPHA2", "ZFP36L1", "DNMT3A", "AJUBA")   # <- edit, or NULL to use a 'highlight' column

col_all  <- "gray70"
col_high <- "#C0392B"
chroms   <- paste0("chr", c(1:22, "X", "Y"))

# ---- Gene list -> hg38 coordinates -------------------------------------------
genes_df <- read.delim(gene_file, stringsAsFactors = FALSE)
symbols  <- unique(trimws(genes_df$gene))
if (is.null(highlight_genes) && "highlight" %in% names(genes_df))
  highlight_genes <- genes_df$gene[as.logical(genes_df$highlight)]

entrez <- mapIds(org.Hs.eg.db, symbols, "ENTREZID", "SYMBOL", multiVals = "first")
miss <- names(entrez)[is.na(entrez)]
if (length(miss)) {
  rescued <- tryCatch(mapIds(org.Hs.eg.db, miss, "ENTREZID", "ALIAS", multiVals = "first"),
                      error = function(e) setNames(rep(NA_character_, length(miss)), miss))
  entrez[miss] <- rescued[miss]
}
entrez <- entrez[!is.na(entrez)]

all_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = FALSE)
gr <- unlist(range(all_genes[intersect(entrez, names(all_genes))]))
gr <- gr[as.character(seqnames(gr)) %in% chroms]
gr <- gr[!duplicated(names(gr))]

genes_pos <- data.frame(
  chr       = as.character(seqnames(gr)),
  mid       = (start(gr) + end(gr)) / 2 / 1e6,          # Mb
  symbol    = names(entrez)[match(names(gr), entrez)]
)
genes_pos$x         <- match(genes_pos$chr, chroms)
genes_pos$highlight <- genes_pos$symbol %in% highlight_genes
high <- genes_pos[genes_pos$highlight, ]

message(nrow(genes_pos), " of ", length(symbols), " genes placed.")
np <- setdiff(symbols, genes_pos$symbol)
if (length(np)) message("Not placed: ", paste(np, collapse = ", "))
bh <- setdiff(highlight_genes, genes_pos$symbol)
if (length(bh)) warning("Highlight genes not found: ", paste(bh, collapse = ", "))

# ---- Cytobands ----------------------------------------------------------------
cyto <- as.data.frame(getCytobands("hg38"))
cyto <- cyto[as.character(cyto$seqnames) %in% chroms, ]
cyto$x     <- match(as.character(cyto$seqnames), chroms)
cyto$start <- cyto$start / 1e6
cyto$end   <- cyto$end / 1e6

chr_len <- aggregate(end ~ x, data = cyto, FUN = max)

stain_cols <- c(gneg = "white", gpos25 = "grey85", gpos33 = "grey80",
                gpos50 = "grey65", gpos66 = "grey55", gpos75 = "grey45",
                gpos100 = "grey30", acen = "#E6B0AA", gvar = "grey90",
                stalk = "grey90")

w <- 0.09   # narrower chromosomes = more room for labels

genes_pos$face <- ifelse(genes_pos$highlight, "bold.italic", "italic")

p <- ggplot() +
  # bands
  geom_rect(data = cyto,
            aes(xmin = x - w, xmax = x + w, ymin = start, ymax = end, fill = gieStain),
            colour = NA) +
  # chromosome outline
  geom_rect(data = chr_len,
            aes(xmin = x - w, xmax = x + w, ymin = 0, ymax = end),
            fill = NA, colour = "grey20", linewidth = 0.3) +
  # ticks for all genes (red if highlighted)
  geom_segment(data = genes_pos,
               aes(x = x - w - 0.08, xend = x + w + 0.02, y = mid, yend = mid,
                   colour = highlight),
               linewidth = 0.3) +
  # labels for all genes
  geom_text_repel(data = genes_pos,
                  aes(x = x + w + 0.03, y = mid, label = symbol,
                      colour = highlight, fontface = face),
                  size = 2, hjust = 0, direction = "y", nudge_x = 0.08,
                  #size = 1.8, hjust = 0, direction = "y", nudge_x = 0.08,
                  segment.size = 0.15, segment.alpha = 0.6,
                  min.segment.length = 0, box.padding = 0.02, point.padding = 0,
                  force = 0.5, max.iter = 20000, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = c(`FALSE` = col_all, `TRUE` = col_high), guide = "none") +
  scale_fill_manual(values = stain_cols, guide = "none") +
  scale_x_continuous(breaks = seq_along(chroms), labels = sub("chr", "", chroms),
                     position = "top", expand = expansion(add = c(0.3, 0.7))) +
  scale_y_reverse(name = "Mb", breaks = seq(0, 250, 50),
                  expand = expansion(add = c(2, 2))) +
  labs(caption = sprintf("All genes (n = %d)   |   location displayed based on hg38 coordinates",
                         nrow(genes_pos), nrow(high))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid   = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face = "bold", size = 10),
        axis.text.y  = element_text(size = 8, colour = "grey40"),
        axis.title.y = element_text(size = 8, colour = "grey40"),
        plot.caption = element_text(size = 8, colour = "grey30", hjust = 0.5),
        plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

print(p)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(out_file, ".pdf"), p, width = 15, height = 7.5, bg = "transparent")
ggsave(paste0(out_file, ".png"), p, width = 15, height = 7.5, dpi = 300, bg = "transparent")

