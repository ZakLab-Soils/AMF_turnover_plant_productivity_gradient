###Indicator taxa examination to characterize relationship between ASV and VT groups to environmental variables###

library(phyloseq)
library(tidyverse)
library(TITAN2)
library(ggplot2)
library(broom)
library(vegan)
library(viridis)
library(ggpubr)
library(tidytext)

ps.amf.F.prune.shared <- readRDS("ps_amf_F_prune_shared.rds")

samp.data <- sample_data(ps.amf.F.prune.shared)

#Create a taxglom version of the dataset (collapse into VT) and rename
ps.taxglom.vt.shared <- tax_glom(ps.amf.F.prune.shared, taxrank="VTX") #24 VT

vtnames <- tax_table(ps.taxglom.vt.shared)[,7] %>%
  gsub("s__VTX", "VT", .)
taxa_names(ps.taxglom.vt.shared) <- vtnames

#Hellinger transform and remove ASVs with less than 3 occurrences

asv.pa.occ3.names <- ps.amf.F.prune.shared %>% otu_table() %>% decostand(., method="pa") %>% data.frame %>% colSums() %>% subset(.,.>=3) %>% names()

asv.hel.occ3 <- ps.amf.F.prune.shared %>% otu_table() %>% decostand(., method="hellinger") %>% data.frame %>% .[,colnames(.) %in% asv.pa.occ3.names]

#VT do not have have less than 3 occurrences
vt.hel <- data.frame(decostand(otu_table(ps.taxglom.vt.shared), method="hellinger"))

##Run TITAN for N mineralization rate and mean BAI

 titan_results_Nmin <- titan(
   samp.data$Nmin_rate,
   asv.hel.occ3,
   minSplt = 5,
   numPerm = 1000,
   boot = TRUE,
   nBoot = 1000,
   imax = FALSE,
   ivTot = FALSE,
   pur.cut = 0.95,
   rel.cut = 0.95,
   ncpus = 4,
   memory = FALSE
 )

#TITAN cannot use missing values: remove samples that have missing BAI data and ASV < 3 occurrences 

bai.data <- samp.data[complete.cases(samp.data$baimean),]

bai.asv.names <-ps.amf.F.prune.shared  %>% prune_samples(sample_names(.) %in% bai.data$sampleID, .) %>% otu_table() %>% decostand(., method="pa") %>% data.frame %>% colSums() %>% subset(.,.>=3) %>% names()
 
asv.hel.bai.occ3 <- asv.hel.occ3 %>% .[,colnames(.) %in% bai.asv.names] %>% subset(., rownames(.) %in% bai.data$sampleID)

 titan_results_BAI <- titan(
   bai.data$baimean,
   asv.hel.bai.occ3,
   minSplt = 5,
   numPerm = 1000,
   boot = TRUE,
   nBoot = 1000,
   imax = FALSE,
   ivTot = FALSE,
   pur.cut = 0.95,
   rel.cut = 0.95,
   ncpus = 4,
   memory = FALSE
 )
 

#Run TITAN on virtual taxa (VT)
titan_results_Nmin_vt <- titan(
  samp.data$Nmin_rate,
  vt.hel,
  minSplt = 5,
  numPerm = 1000,
  boot = TRUE,
  nBoot = 1000,
  imax = FALSE,
  ivTot = FALSE,
  pur.cut = 0.95,
  rel.cut = 0.95,
  ncpus = 4,
  memory = FALSE
)

saveRDS(titan_results_Nmin_vt, "titan_results_Nmin_vt.rds")

#Remove and samples that do not have BAI data
vt.hel.bai <- vt.hel %>% subset(., rownames(.) %in% bai.data$sampleID)

titan_results_BAI_vt <- titan(
  bai.data$baimean,
  vt.hel.bai,
  minSplt = 5,
  numPerm = 1000,
  boot = TRUE,
  nBoot = 1000,
  imax = FALSE,
  ivTot = FALSE,
  pur.cut = 0.95,
  rel.cut = 0.95,
  ncpus = 4,
  memory = FALSE
)

saveRDS(titan_results_BAI_vt, "titan_results_BAI_vt.rds")

  
# Create tibble data from TITAN results
titan_Nmin <- titan_results_Nmin$sppmax %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") 
titan_BAI <- titan_results_BAI$sppmax %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") 
titan_Nmin_vt <- titan_results_Nmin_vt$sppmax %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") 
titan_BAI_vt <- titan_results_BAI_vt$sppmax %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") 

#Summary for figure creation
titan_summary <- bind_rows(lst(BAI=titan_BAI, Nmin = titan_Nmin, BAIvt=titan_BAI_vt, Nminvt = titan_Nmin_vt), .id="type") %>%
  select(type, id, z.median, filter) %>%

# Add response group
filter(filter > 0) %>%
  mutate(
    response_group = NA,
    response_group = ifelse(filter == 1, "Decreasing", response_group),
    response_group = ifelse(filter == 2, "Increasing", response_group),
    response_group = factor(response_group, levels = c("Increasing", "Decreasing"))
  ) %>%
  
  # Arrange taxa
  mutate(z.median = ifelse(filter == 1, -1 * z.median, z.median)) %>%
  group_by(type) %>%
  arrange(z.median, .by_group=TRUE) %>%
  mutate(id=reorder_within(id, z.median, type)) 

titan_summary$label <- factor(titan_summary$type, labels = c("Z-score for Mean BAI" ~cm^2~year^-1, "Z-score for Mean BAI" ~cm^2~year^-1, "Z-score for Net N mineralization rate (µg N"~g^-1~d^-1~")","Z-score for Net N mineralization rate (µg N"~g^-1~d^-1~")"
))

ASV_TITAN_FIG <- ggplot() +
  # Line
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  
  # Bars
  geom_col(
    data = titan_summary[titan_summary$type %in% c("BAI", "Nmin"),],
    aes(x = id, y = z.median, fill = response_group),
    colour = "black",
    width = 0.8
  ) +
  
  #Multiple panels
  facet_wrap(~label, scales="free", dir="v", strip.position = "bottom",
             labeller = label_parsed
  ) +
  
  ylab(NULL) +
  xlab(NULL) +
  
  scale_fill_grey(
    name = "Taxa Response",
    start = 0.8,
    end = 0
  ) +
  
  scale_x_reordered() +
  
  #Flip axes
  coord_flip() +
  
  # Format
  theme(
    
    strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=12),
    
    # Panels
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    # Titles
    plot.title = element_text(colour = "black", size = 16, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 14),
    
    # Axes
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(colour = "black", size = 9),
    
    # Legend
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 10),
    legend.position = "bottom"
    
  )

annotate_figure(ASV_TITAN_FIG, left = "AMF ASV")
ggsave("ASV_Titan.pdf", height = 11, width=8.5, dpi=600)




VT_TITAN_FIG <- ggplot() +
  # Line
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  
  # Bars
  geom_col(
    data = titan_summary[titan_summary$type %in% c("BAIvt", "Nminvt"),],
    aes(x = id, y = z.median, fill = response_group),
    colour = "black",
    width = 0.8
  ) +
  
  #Multiple panels
  facet_wrap(~label, scales="free", dir="v", strip.position = "bottom",
             labeller = label_parsed
  ) +
  
  ylab(NULL) +
  xlab(NULL) +
  
  scale_fill_grey(
    name = "Taxa Response",
    start = 0.8,
    end = 0
  ) +
  
  scale_x_reordered() +
  
  #Flip axes
  coord_flip() +
  
  # Format
  theme(
    
    strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=12),
    
    # Panels
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    # Titles
    plot.title = element_text(colour = "black", size = 16, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 14),
    
    # Axes
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(colour = "black", size = 9),
    
    # Legend
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 10),
    legend.position = "bottom"
    
  )

annotate_figure(VT_TITAN_FIG, left = "AMF VT")
ggsave("VT_Titan.pdf", height = 11, width=8.5, dpi=600)
