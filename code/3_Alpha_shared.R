### Alpha diversity ###

library(vegan)
library(phyloseq)
library(ggplot2)

#Remove sample 50_acru_211 due to low number of reads after QC and taxonomy filtering
ps.amf.F <- readRDS("ps_amf_F.rds")
ps.amf.F.prune <- prune_samples(sample_names(ps.amf.F) != "50_acru_211", ps.amf.F)
saveRDS(ps.amf.F.prune, "ps_amf_F_prune.rds")

#Create lists of taxa names from tree types for shared phyloseq object
asv.acru <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps.amf.F.prune, sample_data(ps.amf.F.prune)$tree == "acru"))>0, ps.amf.F.prune))
asv.acsa <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps.amf.F.prune, sample_data(ps.amf.F.prune)$tree == "acsa"))>0, ps.amf.F.prune))
shared.asv <- subset(asv.acru, asv.acru %in% asv.acsa)
ps.amf.F.prune.shared <- subset_taxa(ps.amf.F.prune, taxa_names(ps.amf.F.prune) %in% shared.asv)
saveRDS(ps.amf.F.prune.shared, "ps_amf_F_prune_shared.rds")

##Rarefaction curves
ps.amf.F.prune.shared.rare <- rarefy_even_depth(ps.amf.F.prune.shared, rngseed = 49657120, sample.size=min(sample_sums(ps.amf.F.prune.shared)), replace = F)
saveRDS(ps.amf.F.prune.shared.rare, "ps_amf_F_prune_shared_rare.rds")
rarefaction.amf.F.shared <- rarecurve(as(otu_table(ps.amf.F.prune.shared.rare), "matrix"), step = 100)

#Site curves
ps.amf.F.shared.rare.site.merged <- merge_samples(ps.amf.F.prune.shared.rare, "site")
saveRDS(ps.amf.F.rare.site.merged, "ps_amf_F_shared_rare_site_merged.rds")
stand.counts <- c(5,5,10,5,5,5,8,5,10,3,4,5)
rarefaction.amf.F.shared.site <- rarecurve(as(round(otu_table(ps.amf.F.shared.rare.site.merged)/stand.counts), "matrix"), step = 100)

##Richness scatterplot regression (Chao1) by Nmin rate

data.chao.ps.amf.F.prune.shared.rare <- plot_richness(ps.amf.F.prune.shared.rare, x="site", measures = "Chao1")$data
chao.nmin.scatter.data.shared <- data.chao.ps.amf.F.prune.shared.rare[,c(1,11,22,23)]
scatter.nmin.chao.shared <- ggplot(data=chao.nmin.scatter.data.shared, aes(x=Nmin_rate, y=value)) + geom_point(color = "#454545") + geom_smooth(method="lm", color = "#458B00") + xlab(paste0("Net N mineralization rate (µg N g-1 d-1)")) + ylab(paste0("Chao1 Richness ")) + ggtitle("Chao1 species richness across the N mineralization gradient") + theme_bw()

#Stat summary
chao.nmin.shared.summary <- summary(lm(value~Nmin_rate, data=chao.nmin.scatter.data.shared))
Adj.r.square <- round(chao.nmin.shared.summary[9]$adj.r.squared, 2)
lm.pvalue <- format(round(chao.nmin.shared.summary[[4]][8],6), scientific=FALSE)
r2.pvalue <- data.frame(label = c(paste("italic(R^2) == ", Adj.r.square), ifelse(lm.pvalue>0.001, paste("italic(P) == ", lm.pvalue), paste("italic(P) < ", 0.001))))

scatter.nmin.chao.shared + annotate(geom="text", x=0.085,y=c(115,110),size=4, label=r2.pvalue$label, parse=TRUE)

ggsave("Richness_regression_shared.pdf", width = 6, height = 5)
