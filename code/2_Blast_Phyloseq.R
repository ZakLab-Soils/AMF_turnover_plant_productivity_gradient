### NCBI Blast Local Taxonomy Assignment and Phyloseq object creation###

library(stringr)
library(phyloseq)
library(Biostrings)

##Processing data for taxonomy and ASV by sequence dataframe
seqtab.nochim.F <- readRDS("seqtab_nochim_F.rds")
uniquesToFasta(getUniques(seqtab.nochim.F), "seqtab_uniquesF.fasta")

fastaToDF <- function(fastaFile){ 
  dnaSeq <- readBStringSet(fastaFile)
  fasta_df <- data.frame(header=names(dnaSeq), sequence=paste(dnaSeq))
}

seqtab.uniques.ftdf.F <- fastaToDF("seqtab_uniquesF.fasta")

asvF <- data.frame(t(seqtab.nochim.F))
colnames(asvF) <- sample.names

asvF$sequence <- rownames(asvF)
rownames(asvF) <- NULL
nrow(asvF)

#Merge the two together
asvF2 <- merge(seqtab.uniques.ftdf.F, asvF)
View(asvF2)

##Use local blast to filter suspected non-AMF and assign taxonomy in command line outside of R

#Maarjam database (non-type) edited to remove sequences with several N's and those that were too short to give meaningful overlap and used to create local database (in data folder on Github)

#>makeblastdb -in maarjam35_alignment_edited_fasta.txt -dbtype nucl

#>blastn -query seqtab_uniquesF.fasta -db maarjam35_alignment_edited_fasta -out seqtab_uniquesF_blast.tab -max_target_seqs 5 -outfmt 6

blast.output <- read.csv("~/NCBI/blast-2.15.0+/bin/seqtab_uniquesF_blast.tab", sep = "\t", header = FALSE)

colnames(blast.output) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
length(unique(blast.output$qseqid)) #number of ASVs 2361

#For our data, we determined that less than 300 bitscore may be non-AMF reads, so we filtered those out
amf.reads.bitscore <- subset(blast.output, blast.output$bitscore > 299)
length(unique(amf.reads.bitscore$qseqid)) #number of ASVs 499

amf.reads.bitscore <- amf.reads.bitscore[order(amf.reads.bitscore$bitscore, decreasing = TRUE),]

#Pull out first hit and then only grab accession number
amf.top.acc <- amf.reads.bitscore[match(unique(amf.reads.bitscore$qseqid), amf.reads.bitscore$qseqid),]
amf.top.acc <- amf.top.acc[1:2]
amf.top.acc$sseqid <- str_split_fixed(amf.top.acc$sseqid, "\\|", 5)[,4]
colnames(amf.top.acc)[2] <- "ACC"

##Clean up the output to have the Maarjam taxonomy and seqID

#Maarjam qiime file edited to have accession numbers with taxonomy
maarjam.tax <- read.csv("maarjam_database_SSU_edited.txt", header = FALSE, sep=";")
colnames(maarjam.tax) <- c("ACC", "Kingdom", "Phylum", "Class", "Order", "Family", "Genera", "VTX")
amf.tax <- merge(amf.top.acc, maarjam.tax, by="ACC")

colnames(amf.tax)[2] <- "header"
amf.tax$ACC <- NULL

amf.tax.merge <- merge(amf.tax, seqtab.uniques.ftdf.F, by="header")
amf.tax.merge$header <- NULL

# merge ASVs with Taxonomy to create a tax table
asv_tax <- merge(asvF2, amf.tax.merge, by="sequence")
View(asv_tax)

asv_tax$ASV <- asv_tax$header
asv_tax$header <- NULL
asv_tax$ASV <- gsub("sq", "ASV", str_split_fixed(asv_tax$ASV, ";", 2)[,1])

saveRDS(asv_tax, "ASV_Combine_Tax.rds")

rownames(asv_tax) <- asv_tax$ASV
asv_tax$ASV <- NULL

##Phyloseq creation
asv.refreads <- DNAStringSet(asv_tax[,1])
names(asv.refreads) <- rownames(asv_tax)
asv.tbl <- asv_tax[,2:72] #Column 1 is reads, and ASV counts 2-72
asv.tbl <- t(asv.tbl)
tax.tbl <- asv_tax[,73:length(asv_tax)] #Taxonomy in columns 73 to end
soil.metadata <- read.csv("Manistee_maples_metadata.csv")
soil.metadata <- soil.metadata[match(rownames(asv.tbl), soil.metadata$sampleID),]
rownames(soil.metadata) <- soil.metadata$sampleID

ps.amf.F <- phyloseq(otu_table(asv.tbl, taxa_are_rows = FALSE), tax_table(as.matrix(tax.tbl)), sample_data(soil.metadata), refseq(asv.refreads))

#Correcting the tax_table for ASV285 and ASV207 to correct the Genus from Glomus to Diversispora
tax_table(ps.amf.F)[102,6] <- "g__Diversispora"
tax_table(ps.amf.F)[104,6] <- "g__Diversispora"

saveRDS(ps.amf.F, "ps_amf_F.rds")