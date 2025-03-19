library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)

# Load your phyloseq object
# Assuming your object is named ps
ps <- load("senior_thesis/plant-growth-metrics/ITS.phyloseq.RData")

# 1. Alpha Diversity
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
meta <- sample_data(ps)
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- merge(alpha_div, meta, by = "SampleID")

# Visualize alpha diversity
ggplot(alpha_div, aes(x = SoilType, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  labs(title = "Shannon Diversity by Soil Type")

# 2. Beta Diversity
# Ordination with PCoA and Bray-Curtis
ord <- ordinate(ps, method = "PCoA", distance = "bray")
plot_ordination(ps, ord, color = "SoilType") +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA of Microbial Communities")

# PERMANOVA test
bray_dist <- distance(ps, method = "bray")
adonis2(bray_dist ~ SoilType, data = meta)

# 3. Taxonomic Composition
# Aggregate by Phylum for visualization
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
plot_bar(ps_rel, fill = "Phylum") +
  theme_minimal() +
  labs(title = "Relative Abundance by Phylum")

# 4. Differential Abundance Analysis
dds <- phyloseq_to_deseq2(ps, ~ SoilType)
dds <- DESeq(dds)
res <- results(dds, contrast = c("SoilType", "Rincon_silty_clay_loam", "Delhi_loamy_sand"))
sig_res <- res[which(res$padj < 0.05), ]
head(sig_res)

# Visualize significant taxa
sig_taxa <- rownames(sig_res)
ps_sig <- prune_taxa(sig_taxa, ps)
plot_heatmap(ps_sig, sample.label = "SoilType", taxa.label = "Genus")
