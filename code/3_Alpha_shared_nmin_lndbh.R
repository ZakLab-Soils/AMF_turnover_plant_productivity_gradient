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
chao.scatter.data.shared <- data.chao.ps.amf.F.prune.shared.rare[,c(1,5,11,14,15)]

#Stats for regressions for selected variables
variables <- c("Nmin_rate", "lndbh")
formulas <- list()
for(i in seq_along(variables)){
  formulas[[i]] <- lm(formula(paste("value", "~", paste(variables[i]))), data=chao.scatter.data.shared)
}
glance_results <- data.frame(bind_rows(lapply(formulas, glance)))
glance_results$variable <- variables

r2.pvalue.label <- glance_results %>% data.frame() %>% select(adj.r.squared, p.value, variable) %>% mutate(labelr2  = paste("italic(R^2) ==", sprintf("%.2f", adj.r.squared))) %>% mutate(labelpv = ifelse(p.value >0.001, paste("italic(P) == ", sprintf("%.3f", p.value)), paste("italic(P) < ", 0.001)))

#Richness scatterplot regression (Chao1) by N mineralization rate
scatter.nmin.chao.shared <- ggplot(data=chao.scatter.data.shared, aes(x=Nmin_rate, y=value)) + geom_point(color = "#454545") + geom_smooth(method="lm", color = "#458B00") + xlab(paste0("Net N mineralization rate (µg N g-1 d-1)")) + ylab(paste0("Chao1 Richness ")) + ggtitle("Chao1 species richness across the N mineralization gradient") + theme_bw() + annotate(geom="text", x=0.085,y=c(115,110),size=4, label=c(r2.pvalue.label[1,4], r2.pvalue.label[1,5]), parse=TRUE) 

scatter.nmin.chao.shared

ggsave("Richness_Nmin_regression_shared.pdf", width = 6, height = 5)

##Richness scatterplot regression (Chao1) by tree growth/size (lnDBH)

scatter.chao.lndbh <- ggplot(data=chao.lndbh, aes(x=lndbh, y=value)) + geom_point(color = "#454545") + geom_smooth(method="lm", color = "#458B00") + xlab(paste0("lnDBH")) + ylab(paste0("Chao1 Richness ")) + ggtitle("Chao1 species richness along tree size (lnDBH)") + theme_bw() + annotate(geom="text", x=2.7,y=c(90,85),size=4, label=c(r2.pvalue.label[2,4], r2.pvalue.label[2,5]), parse=TRUE)

scatter.chao.lndbh

ggsave("Richness_lnDBH_regression_shared.pdf", width = 6, height = 5)
