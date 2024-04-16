###Multivariate GLMs to examine relationship between ASV and VT groups to environmental variables###

library(mvabund)
library(ggplot)
library(phyloseq)
library(tidyverse)
library(broom)

ps.amf.F.prune.shared <- readRDS("ps_amf_F_prune_shared.rds")

#Create a taxglom version of the dataset (collapse to VT)
ps.taxglom.vt.shared <- tax_glom(ps.amf.F.prune.shared, taxrank="VTX")
taxa_names(ps.taxglom.vt.shared) <- tax_table(ps.taxglom.vt.shared)[,7] #Rename taxa to VT names

#Filtering out low abundance ASVs - only for ASVs but will not use this object to taxglom
ps.tsc <- transform_sample_counts(ps.amf.F.prune.shared, function(x)x/sum(x))
ps.tsc.ft <- filter_taxa(ps.tsc, function(x) mean(x) > .01, TRUE) #>1% of the data
top1.asv <- taxa_names(ps.tsc.ft)
ps.amf.F.prune.shared.top1 <- subset_taxa(ps.amf.F.prune.shared, taxa_names(ps.amf.F.prune.shared) %in% top1.asv)

#Grab the count data (otu table)
shared.otu.top1 <- as.data.frame(otu_table(ps.amf.F.prune.shared.top1))
shared.otu.vt <- as.data.frame(otu_table(ps.taxglom.vt.shared))

#Create object that is recognized as multivariate
otu.top1.mvabund <- mvabund(shared.otu.top1)
otu.vt.mvabund <- mvabund(shared.otu.vt)

#Sample data which has the environmental variables
samp.data <- data.frame(sample_data(ps.amf.F.prune.shared)[,c(1,2,11,14)])

#Model the multivariate otu table with glm by our environmental variable (N mineralization). 
otu.top1.manyglm <- manyglm(otu.top1.mvabund~samp.data$Nmin_rate, family="negative.binomial")
otu.vt.manyglm <- manyglm(otu.vt.mvabund~samp.data$Nmin_rate, family="negative.binomial")

plot(otu.top1.manyglm) #no pattern
plot(otu.vt.manyglm) #no pattern

#Test model with ANOVA and adjust pvalues for multiple comparisons
set.seed(125689)
glm.mod.nmin.anova.top1.puni <- anova.manyglm(otu.top1.manyglm, p.uni="adjusted", nBoot=999)

set.seed(125689)
glm.mod.nmin.anova.vt.puni <- anova.manyglm(otu.vt.manyglm, p.uni="adjusted", nBoot=999)

#Need to pull out the adjusted pvalues
top1.glm.unip <- setNames(data.frame(glm.mod.nmin.anova.top1.puni$uni.p[2,]), "p_val")
top1.glm.unip$OTU <- rownames(top1.glm.unip)
vt.unip <- setNames(data.frame(glm.mod.nmin.anova.vt.puni$uni.p[2,]), "p_val")
vt.unip$VTX <- rownames(vt.unip)

##Graphing ASVs and VTs with significant responses to Nmin rate##

#Hellinger transformed counts for figures 
otu.table.shared.top1.relab <- decostand(otu_table(ps.amf.F.prune.shared.top1), method = "hellinger")
otu_table(ps.amf.F.prune.shared.top1) <- otu_table(otu.table.shared.top1.relab, taxa_are_rows = FALSE)

ps.top1.psmelt <- psmelt(ps.amf.F.prune.shared.top1)

#Nest data and map models
psmelt.simple <- ps.top1.psmelt %>% 
  as_tibble() %>% 
  select(Sample = Sample, nmin = Nmin_rate, abund = Abundance, genus = Genera, VTX = VTX, OTU = OTU)

psmelt.simple_asv <- psmelt.simple %>% 
  group_by(OTU) %>% 
  nest() %>% 
  mutate(
    mod = map(data, ~ glm(abund ~ nmin, data = .))) %>%
  inner_join(top1.glm.unip, by="OTU")

psmelt.simple_asv %>% 
  select(-mod) %>% 
  unnest(data) %>% 
  filter(p_val <= 0.055) %>%
  ggplot(aes(x = nmin, y = abund)) +
  geom_point(color = "#454545") +
  geom_smooth(method = "glm", color = "#458B00") +
  xlab(paste0("Net N mineralization (µg N g-1 d-1)")) + ylab(paste0("Relative Abundance (%)")) + 
  ggtitle("Relative abundance of ASV's across the N mineralization gradient") +
  theme_bw() +
  facet_wrap(~ OTU, scales = "free_y")

ggsave("Relative_Abundance_Nmin_Shared_Top1_ASV.pdf", height = 11, width=8.5)

#Now examine the VT from the taxglom data
otu.table.vtx.relab.shared <- decostand(otu_table(ps.taxglom.vt.shared), method = "hellinger")
ps.taxglom.vt.shared2<- ps.taxglom.vt.shared
otu_table(ps.taxglom.vt.shared2) <- otu_table(otu.table.vtx.relab.shared, taxa_are_rows = FALSE)

ps.taxglom.vt.shared.psmelt <- psmelt(ps.taxglom.vt.shared2)

#Nest data and map models
ps.taxglom.vt.shared.psmelt.simple <- ps.taxglom.vt.shared.psmelt %>% 
  as_tibble() %>% 
  select(Sample = Sample, nmin = Nmin_rate, abund = Abundance, genus = Genera, VTX = VTX, OTU = OTU)

ps.taxglom.vt.shared.psmelt.simple_vt <- ps.taxglom.vt.shared.psmelt.simple %>% 
  group_by(VTX) %>% 
  nest() %>% 
  mutate(
    mod = map(data, ~ glm(abund ~ nmin, data = .))) %>%
  inner_join(vt.unip, by="VTX")

ps.taxglom.vt.shared.psmelt.simple_vt %>% 
  select(-mod) %>% 
  unnest(data) %>% 
  filter(p_val <= 0.05) %>%
  ggplot(aes(x = nmin, y = abund)) +
  geom_point(color = "#454545") +
  geom_smooth(method = "glm", color = "#458B00") +
  xlab(paste0("Net N mineralization (µg N g-1 d-1)")) + ylab(paste0("Relative Abundance (%)")) + 
  ggtitle("Relative abundance of VT's across the N mineralization gradient") +
  theme_bw() +
  facet_wrap(~ VTX, scales = "free_y")

ggsave("Relative_Abundance_Nmin_Shared_VT.pdf", height = 11, width=8.5)
