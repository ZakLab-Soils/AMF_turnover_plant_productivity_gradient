###Analyzing Manistee maple 18S AMF reads###

##Intial Quality filtering and primer removal
#We are only using forward reads as our reads do not overlap and the quality was better for our forward reads

library(dada2)
library(ShortRead)
library(Biostrings)

#Set working directory to where the reads are located
setwd("~/Manistee_AMF_reads")

fnFs <- sort(list.files(".", pattern = "_R1_001.fastq.gz", full.names = TRUE))

#Quality plots of the data for F (can be done in smaller batches for ease of reading)

plotQualityProfile(fnFs)
#plotQualityProfile(fnFs[17:20])

fnFs.filtN <- file.path(".", "filtN", basename(fnFs))
noN.ft <- filterAndTrim(fnFs, fnFs.filtN, maxN=0, multithread = TRUE)

#Remove primer contaminantion with cutadapt - must be installed locally  
cutadapt <-"~/cutadapt"
path.cut <- file.path(".", "cutadapt_processed")
if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))

FWD <- "TTGGAGGGCAAGTCTGGTGCC" #NS31
REV <- "GAACCCAAACACTTTGGTTTCC" #AML2
REV.RC <- dada2:::rc(REV)
FWD.RC <- dada2:::rc(FWD)

for(i in seq_along(fnFs)){system2(cutadapt, args = c("-a", FWD, "-a", REV.RC, "-n", 2, "-o", fnFs.cut[i], fnFs.filtN[i]))}

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#Check to see how much primer removed by sample 
#(so [[1]] needs to change number for each individual library)

rbind(FWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnFs.filtN[[1]]))

rbind(FWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnFs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))

##Filter and trim for length, quality and size; Calculate error rate.
get.sample.name <- function(fname) strsplit(basename(fname), "_R1_001.fastq.gz")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))

#Given these parameters, we have only 2 potentially problematic samples: 24_acru_290 & 24_acsa_248

filter.summary <- filterAndTrim(cutFs, filtFs, maxN = 0, minLen = 200, truncLen = 240, maxEE = 1.75, truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

head(filter.summary)

#Check to see if sequences are improved
plotQualityProfile(filtFs)

Perc_reads_retained <- filter.summary[,2]/filter.summary[,1]
saveRDS(Perc_reads_retained, "Perc_reads_retained.rds")

#Calculate error rate
errF <- learnErrors(filtFs, multithread = TRUE)

saveRDS(errF, file = "errF_AMF.rds")

errF.plot <- plotErrors(errF, nominalQ = TRUE)
errF.plot

##ASV creation and chimera removal
derepFs <- derepFastq(filtFs, verbose = TRUE)

names(derepFs) <- sample.names
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

seqtab.F <- makeSequenceTable(dadaFs)
saveRDS(seqtab.F, file = "seqtab.F.rds")

dim(seqtab.F) #gives number of samples and number of ASVs - 71 samples & 2545 ASVs

#Chimera removal (denovo)
seqtab.nochim.F <- removeBimeraDenovo(seqtab.F, method = "consensus", multithread = TRUE, verbose = TRUE) #Identified 125 bimeras out of 2545 input sequences.
saveRDS(seqtab.nochim.F, "seqtab_nochim_F.rds")

