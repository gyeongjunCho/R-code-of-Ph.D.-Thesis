require(dada2)
require(doParallel)
registerDoParallel(cores=64) # CPU threads to multithreads
library(ShortRead)
library(Biostrings)
require(DECIPHER)
require(ape)
require(vegan)
require(picante)
require(progress)
require(reshape2)
require(ggplot2)
require(tidyverse)
require(ggh4x)
require(ggalt)
require(ggnetwork)
require(RColorBrewer)
library(yyplot)
require(ggvenn)
require(ggpubr)
require(iNEXT)
require(ggrepel)
#===============
# DADA2 analysis
#===============

# Input raw files directory
path_c <- "cycle_raw"
path_a <- "amino_raw"
path_f <- "field_raw"

# Forward reads list
fnFs_c <- sort(list.files(path_c, pattern="_1.fastq.gz", full.names = TRUE))
fnFs_a <- sort(list.files(path_a, pattern="_1.fastq.gz", full.names = TRUE))
fnFs_f <- sort(list.files(path_f, pattern="_1.fastq.gz", full.names = TRUE))

# Reverse reads list
fnRs_c <- sort(list.files(path_c, pattern="_2.fastq.gz", full.names = TRUE))
fnRs_a <- sort(list.files(path_a, pattern="_2.fastq.gz", full.names = TRUE))
fnRs_f <- sort(list.files(path_f, pattern="_2.fastq.gz", full.names = TRUE))

raw_foward_Q <- plotQualityProfile(c(fnFs_c,fnFs_a,fnFs_f))
raw_foward_Q

raw_reverse_Q <- plotQualityProfile(c(fnRs_c,fnRs_a,fnRs_f))
raw_reverse_Q

# primer info
P515F <- "GTGYCAGCMGCCGCGGTAA"  
P805R <-  "GACTACHVGGGTATCTAATCC"
P806R <- "GGACTACNVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

P515F.orients <- allOrients(P515F)
P805R.orients <- allOrients(P805R)
P806R.orients <- allOrients(P806R)

dir.create(paste0(path_c, "/primer_trim"))
dir.create(paste0(path_f, "/primer_trim"))
dir.create(paste0(path_a, "/primer_trim"))

#Remove Primers
cutadapt <- "/opt/miniconda3/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut_c <- file.path(path_c, "primer_trim")
if(!dir.exists(path.cut_c)) dir.create(path.cut_c)
fnFs.cut_c <- file.path(path.cut_c, basename(fnFs_c))
fnRs.cut_c <- file.path(path.cut_c, basename(fnRs_c))

path.cut_a <- file.path(path_a, "primer_trim")
if(!dir.exists(path.cut_a)) dir.create(path.cut_a)
fnFs.cut_a <- file.path(path.cut_a, basename(fnFs_a))
fnRs.cut_a <- file.path(path.cut_a, basename(fnRs_a))

path.cut_f <- file.path(path_f, "primer_trim")
if(!dir.exists(path.cut_f)) dir.create(path.cut_f)
fnFs.cut_f <- file.path(path.cut_f, basename(fnFs_f))
fnRs.cut_f <- file.path(path.cut_f, basename(fnRs_f))


P515F.RC <- dada2:::rc(P515F)
P805R.RC <- dada2:::rc(P805R)
P806R.RC <- dada2:::rc(P806R)

# Run Cutadapt
R1.flags_c <- paste("-g", P515F, "-a", P806R.RC) 
R2.flags_c <- paste("-G", P806R, "-A", P515F.RC) 

for(i in seq_along(fnFs_c)) {
  system2(cutadapt, args = c(R1.flags_c, R2.flags_c, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut_c[i], "-p", fnRs.cut_c[i], # output files
                             fnFs_c[i], fnRs_c[i])) # input files
}

R1.flags_a <- paste("-g", P515F, "-a", P805R.RC) 
R2.flags_a <- paste("-G", P805R, "-A", P515F.RC) 

for(i in seq_along(fnFs_a)) {
  system2(cutadapt, args = c(R1.flags_a, R2.flags_a, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut_a[i], "-p", fnRs.cut_a[i], # output files
                             fnFs_a[i], fnRs_a[i])) # input files
}

R1.flags_f <- paste("-g", P515F, "-a", P805R.RC) 
R2.flags_f <- paste("-G", P805R, "-A", P515F.RC) 

for(i in seq_along(fnFs_f)) {
  system2(cutadapt, args = c(R1.flags_f, R2.flags_f, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut_f[i], "-p", fnRs.cut_f[i], # output files
                             fnFs_f[i], fnRs_f[i])) # input files
}



path_ct <- "cycle_raw/primer_trim"
path_at <- "amino_raw/primer_trim"
path_ft <- "field_raw/primer_trim"

# Forward reads list
fnFs <- sort(list.files(c(path_ct, path_at, path_ft), pattern="_1.fastq.gz", full.names = TRUE))

# Reverse reads list
fnRs <- sort(list.files(c(path_ct, path_at, path_ft), pattern="_2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Extract raw files list
filt_path <- "filtered_fastq_files"

filtFs <- file.path(filt_path , paste0(sample.names, "_F_filt.fastq.gz")) # set filtered forward path info
filtRs <- file.path(filt_path , paste0(sample.names, "_R_filt.fastq.gz")) # set filtered reverse path info

names(filtFs) <- sample.names
names(filtRs) <- sample.names

filterAndTrim(fwd = fnFs, # forward reads
              filt = filtFs, # filtered forward reads save info
              rev = fnRs, # reverse reads 
              filt.rev = filtRs, # filtered reverse reads save info
              truncLen=c(250,190), # trim to qualify score 30 or higher
              maxN = 0, maxEE = c(1,1), rm.phix = TRUE, # Other statistical quality check criteria
              n=1e8, # sampling reads number
              compress = TRUE, # compress
              multithread = FALSE) # On Windows set multithread=FALSE; An error has been reported when the operation has an unacceptable amount of RAM. In this case, "FALSE" is recommended.

# Visualization filtered forward quality
filt_foward_Q <- plotQualityProfile(filtFs)
filt_foward_Q

# Visualization filtered reverse quality
filt_reverse_Q <- plotQualityProfile(filtRs)
filt_reverse_Q

gc()

errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+11)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e+11)
errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errR, nominalQ=TRUE)
errF_plot
errR_plot

# Sample Inference using Divisive Amplicon Denoising Algorithm (DADA)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 200:300]
nrow(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
seqtab.nochim


#===================
# Clustering and Alignment
#===================
# Assign taxonomy with IDTAXA
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet
load("IDTAXA SILVA DB/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
rm(trainingSet) # remove SILVA DB on ram for ram management

sample_ord <- sample.names[c(11:15, 21:25, 16:20, 6:10, 1:5, 26:30, 35:70, 31:34 ,71:85)]
sample_ord
# No-bacterial sq removing
taxid_sqtable <- cbind(data.frame(taxid), as.data.frame(t(seqtab.nochim)[,sample_ord]))
taxid_sqtable <- subset(taxid_sqtable, domain == "Bacteria")
taxid_sqtable <- subset(taxid_sqtable, order != "Chloroplast")
taxid_sqtable <- subset(taxid_sqtable, family != "Mitochondria")
taxid_sqtable

# sq ID generation
sqID <- paste0("Sq_", c(1:nrow(taxid_sqtable)))

# OTU fasta files
sq_fasta <- foreach(i= 1:nrow(taxid_sqtable), .combine="rbind")%dopar%{
  match <- rbind(sqID[i], rownames(taxid_sqtable)[i])
  print(match)
} 

sq_fasta <- as.vector(sq_fasta)
sq_fasta <- as.data.frame(gsub("^Sq_", ">Sq_",sq_fasta))
colnames(sq_fasta) <- NA
write.table(sq_fasta, file = "metagenome_analysis_results/DADA2_results/Sq_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

# Load fast file
sq_fasta_path <- "metagenome_analysis_results/DADA2_results/Sq_sequence.fasta"
sq_seqs <- readDNAStringSet(sq_fasta_path)
sq_seqs <- OrientNucleotides(sq_seqs)

# Alignment
sq_aligned_seqs <- AlignSeqs(sq_seqs)
BrowseSeqs(sq_aligned_seqs, highlight = 0)
system("mkdir 'metagenome_analysis_results/metagenomic sequence pylogenetic tree'")
writeXStringSet(sq_aligned_seqs, file="metagenome_analysis_results/metagenomic sequence pylogenetic tree/sq_aligned_sequence.fasta")

aligned_sq <- read.csv("metagenome_analysis_results/metagenomic sequence pylogenetic tree/sq_aligned_sequence.fasta")

paste(aligned_sq[,1])

aligned_sq2 <- c()
for(i in 1:nrow(aligned_sq)){
  aligned_sq2 <- paste0(aligned_sq2, aligned_sq[i,1])
}

aligned_sq2 <- gsub("[0-9]", "", aligned_sq2)
aligned_sq2 <- strsplit(aligned_sq2, split=">Sq_")
aligned_sq2 <- unlist(aligned_sq2)

for(i in 1:length(aligned_sq2)){
aligned_sq2[i] <- substr(aligned_sq2[i], 1, 319)
}
nchar(aligned_sq2)
OTU_sq <- gsub("-" ,"" ,aligned_sq2)
taxid_sqtable$OTU <- OTU_sq

for(i in 1:length(sample_ord))

OTUtable <- taxid_sqtable[,c("OTU", sample_ord)]
OTUtable <- melt(OTUtable)

OTUtable <- OTUtable%>%
  group_by(OTU, variable)%>%
  summarise("read"=sum(value))
OTUtable <- dcast(OTUtable, OTU ~ variable)
rownames(OTUtable) <- OTUtable$OTU
OTUtable <- OTUtable[,which(colnames(OTUtable)!="OTU")]
OTUtable <- OTUtable[order(rowSums(OTUtable), decreasing = T), ]


# Assign taxonomy with IDTAXA
dna <- DNAStringSet(rownames(OTUtable)) # Create a DNAStringSet
load("IDTAXA SILVA DB/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
rm(trainingSet) # remove SILVA DB on ram for ram management

# OTUs
taxid_OTUtable <- cbind(data.frame(taxid), as.data.frame(OTUtable[,sample_ord]))
taxid_OTUtable <- subset(taxid_OTUtable, domain == "Bacteria")
taxid_OTUtable <- subset(taxid_OTUtable, order != "Chloroplast")
taxid_OTUtable <- subset(taxid_OTUtable, family != "Mitochondria")

taxid_OTUtable

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

inputF <- raw_foward_Q$layers[[6]]$data$rc
names(inputF) <- raw_foward_Q$layers[[6]]$data$label

inputR <- raw_reverse_Q$layers[[6]]$data$rc
names(inputR) <- raw_reverse_Q$layers[[6]]$data$label

filteredF <- filt_foward_Q$layers[[6]]$data$rc
names(filteredF) <- filt_foward_Q$layers[[6]]$data$label

filteredR <- filt_reverse_Q$layers[[6]]$data$rc
names(filteredR) <- filt_reverse_Q$layers[[6]]$data$label

inoutput <- data.frame(inputF,inputR,filteredF,filteredR)

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
track <- track[sample_ord,]
tax_Bacteria <- colSums(taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
track <- cbind(inputF, inputR, filteredF, filteredR, track, tax_Bacteria)
rownames(track) <- sample_ord

system("mkdir metagenome_analysis_results")
system("mkdir metagenome_analysis_results/DADA2_results")
write.csv(track, "metagenome_analysis_results/DADA2_results/Reads_track.csv")

#===========================
# DADA2 results output files
#===========================


# OTUs ID generation
OTUsID <- foreach(i = 1:nrow(taxid_OTUtable), .combine="rbind")%dopar%{
  OTUsN <- paste("OTU", i, sep="_")
  print(OTUsN)
}
OTUsID <- as.vector(OTUsID)  
OTUsID 

# OTU fasta files
fasta <- foreach(i= 1:nrow(taxid_OTUtable), .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], rownames(taxid_OTUtable)[i])
  print(match)
} 

fasta <- as.vector(fasta)
fasta <- as.data.frame(gsub("^OTU_", ">OTU_",fasta))
colnames(fasta) <- NA
write.table(fasta, file = "metagenome_analysis_results/DADA2_results/OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

# OTU counts table
OTU_counts_table <- cbind(OTUsID, taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
dim(OTU_counts_table)
row.names(OTU_counts_table) <- NULL
colnames(OTU_counts_table)[1] <- "#OTU ID"
write.table(OTU_counts_table, file = "metagenome_analysis_results/DADA2_results/OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# OTU tax table
OTU_tax_table <- cbind(OTUsID, taxid_OTUtable[,1:7])
row.names(OTU_tax_table) <- NULL
colnames(OTU_tax_table)[1] <- "#OTU ID"
write.table(OTU_tax_table, file = "metagenome_analysis_results/DADA2_results/OTU_tax_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# sequence tax counts table
OTU_seq_tax_counts_table <- cbind(as.vector(row.names(taxid_OTUtable)), OTUsID, taxid_OTUtable)
row.names(OTU_seq_tax_counts_table) <- NULL
colnames(OTU_seq_tax_counts_table)[1] <- "sequence" 
write.table(OTU_seq_tax_counts_table, file = "metagenome_analysis_results/DADA2_results/OTU_seq_tax_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ++++++++++++++++++++++++++++++++++
# Microbiome transform data anlaysis
# ++++++++++++++++++++++++++++++++++
require(phyloseq)
# OTU
OTU <- otu_table(
  t(
    read.csv("metagenome_analysis_results/DADA2_results/OTU_counts_table.tsv", sep = "\t",row.names = 1)
  ),
  taxa_are_rows = FALSE
)
rownames(OTU) <- sample_ord

# Taxa
TAX <- tax_table(
  as.matrix(
    read.csv("metagenome_analysis_results/DADA2_results/OTU_tax_table.tsv", sep = "\t",row.names = 1)
  )
)

# Sample inforamtion
SAM <- sample_data(
  read.csv("phenotype_recode/sample_information.tsv", sep="\t", row.names = 1)
)
SAM

# To make phyloseq obj
# reads base obj
phy_obj_reads <- phyloseq(OTU,TAX,SAM)

# to modified relative OTU phyloseq obj
phy_obj_relative <- phy_obj_reads
phy_obj_relative@otu_table <- phy_obj_reads@otu_table/rowSums(phy_obj_reads@otu_table)
write.csv(t(phy_obj_relative@otu_table[,paste0("OTU_",1:ncol(phy_obj_relative@otu_table))]), "metagenome_analysis_results/DADA2_results/relative_abundance_by_OTU.csv")


rarefaction_curve_data <- ggiNEXT( 
  iNEXT(as.data.frame(
    t(OTU)
  )
  )
)$data


sam_info <- data.frame(SAM)
sam_info$site <- rownames(sam_info)
rarefaction_curve_data <- merge(rarefaction_curve_data, sam_info, by="site")

rarefaction_curve_data$Cycling <- factor(rarefaction_curve_data$Cycling, levels=unique(rarefaction_curve_data$Cycling)[c(1:5,7:15,6,16)])
rarefaction_curve_data$Cycling
rarefaction_curve_max_data <- aggregate(x ~ site + lty + Cycling + cultivation + Site_ID, subset(rarefaction_curve_data, lty=="interpolated"), max)
rarefaction_curve_max_data$y <- aggregate(y ~ site + lty + Cycling + cultivation + Site_ID, subset(rarefaction_curve_data, lty=="interpolated"), max)$y
rarefaction_curve_max_data$Cycling <- factor(rarefaction_curve_max_data$Cycling, levels=unique(rarefaction_curve_data$Cycling)[c(1:5,7:15,6,16)])

rarefaction_partA <-
ggplot(data = subset(rarefaction_curve_data, cultivation =="field"), mapping = aes(x=x , y=y, color=site, linetype=lty, label=site))+
  geom_line()+
  geom_point(data = subset(rarefaction_curve_max_data, cultivation =="field"))+
  theme_light()+
  scale_x_continuous(limits=c(0,110000))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  xlab("")+
  ylab("")+
  facet_wrap(Cycling~., ncol=5)+
  labs(linetype="Method")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.text = element_text(color="black"))

rarefaction_partB <-
  ggplot(data = subset(rarefaction_curve_data,cultivation =="pot" & Cycling != "10% bulk soil from 10th"), mapping = aes(x=x , y=y, color=site, linetype=lty, label=site))+
  geom_line()+
  geom_point(data = subset(rarefaction_curve_max_data,cultivation =="pot" & Cycling != "10% bulk soil from 10th"))+
  theme_light()+
  scale_x_continuous(limits=c(0,110000))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  xlab("")+
  ylab("Inferenced OTU number")+
  ggh4x::facet_nested_wrap(Cycling~., nrow =2 )+
  labs(linetype="Method")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.text = element_text(color="black"))

rarefaction_partC <-
  ggplot(data = subset(rarefaction_curve_data, Cycling == "10% bulk soil from 10th"), mapping = aes(x=x , y=y, color=site, linetype=lty, label=site))+
  geom_line()+
  geom_point(data = subset(rarefaction_curve_max_data, Cycling == "10% bulk soil from 10th"))+
  theme_light()+
  scale_x_continuous(limits=c(0,110000))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  xlab("Reads")+
  ylab("  ")+
  ggh4x::facet_nested_wrap(Cycling+Site_ID~., nrow =1 )+
  labs(linetype="Method")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.text = element_text(color="black"))

ggarrange(rarefaction_partA, rarefaction_partB, rarefaction_partC, ncol=1, labels = c("A","B","C"), heights = c(1,1.5,1.1))
#---------
# Coverage 
#----------
require(entropart)
OTU_long_table <- data.frame(t(OTU))

# Chao's coverage
chao_coverage <-  unlist(
  foreach(i = 1:dim(OTU_long_table)[2])%dopar%{
    Coverage(OTU_long_table[,i], Estimator = "Chao")
  }
)
names(chao_coverage) <- colnames(OTU_long_table)

# Good's coverage
good_coverage <- unlist(
  foreach(i = 1:dim(OTU_long_table)[2])%dopar%{
    1 - (sum(OTU_long_table[,i] == 1) / sum(OTU_long_table[,i]))
  }
)

# Coverage summary
coverage <- data.frame(good_coverage, chao_coverage)*100
colnames(coverage) <- c("Good's coverage (%)", "Chao's coverage (%)")
write.csv(coverage, "metagenome_analysis_results/coverage.csv")

#------------------
# Alphyadiversity
#----------------
# Calcuation
alpha_diversity <- plot_richness(phy_obj_reads)$data
alpha_diversity_index  <- subset(alpha_diversity, variable=="Observed"|variable=="Shannon"|variable=="Simpson")

shapiro.test(subset(alpha_diversity_index, variable == "Observed"&design =="field")$value)

bartlett.test(value ~ design, subset(alpha_diversity_index, variable == "Observed"))
bartlett.test(value ~ design, subset(alpha_diversity_index, variable == "Shannon"))
bartlett.test(value ~ design, subset(alpha_diversity_index, variable == "Simpson"))


alpha_observed_aov <-aov(value ~ design, subset(alpha_diversity_index, variable == "Observed"))
alpha_observed_aov_summary <- summary(alpha_observed_aov)[1]

unlist(alpha_observed_aov_summary)

kruskal.test(value ~ design, subset(alpha_diversity_index, variable == "Shannon"))

kruskal.test(value ~ design, subset(alpha_diversity_index, variable == "Simpson"))

alpha_statistics <- 
data.frame(
variable = c("Observed", "Shannon", "Simpson"),
label1 = c("ANOVA", "Kruskal-Wallis test", "Kruskal-Wallis test"),
pval = c("italic(P)==5.070%*%10^-41", "italic(P)==1.112%*%10^-14", "italic(P)==7.113%*%10^-12"),
x = c(0.5,0.5,0.5),
y = c(1300, 6.5, 1.1)
)

alpha_observed_posthoc<- TukeyHSD(alpha_observed_aov)
alpha_observed_posthoc$design

require(nparcomp)
nparcomp(value ~ design, subset(alpha_diversity_index, variable == "Shannon"))
nparcomp(value ~ design, subset(alpha_diversity_index, variable == "Simpson"))

alpha_posthoc <- data.frame(
  variable = c(rep("Observed", 3), rep("Shannon", 3), rep("Simpson", 3)),
  x1 = c("pot+cycling", "pot+nitrogen", "pot+nitrogen",
         "pot+cycling", "pot+nitrogen", "pot+nitrogen",
         "pot+cycling", "pot+nitrogen", "pot+nitrogen"
         ),
  
  x2 = c("field", "field", "pot+cycling",
         "field", "field", "pot+cycling",
         "field", "field", "pot+cycling"
         ),
  y=c(3:1*60+1200,
      3:1*0.3+6,
      3:1*0.05+1)
)

alpha_posthoc$sig <- c(rep("***",6), "**", "*", "***")


# visualization
alpha_diversity_p <- ggplot(data = alpha_diversity_index , mapping = aes(x = design, y=value))+
  stat_summary(geom = "bar", fun.data = "mean_se")+
  stat_summary(geom = "errorbar", width=0.1, fun.data = "mean_se")+
  geom_text(data=alpha_statistics, aes(x=x, y=y*1.2, label=label1), hjust=0, vjust=2)+
  geom_text(data=alpha_statistics, aes(x=x, y=y*1.185, label=pval), hjust=0, vjust=2, parse = T)+
  geom_segment(data=alpha_posthoc, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=alpha_posthoc, aes(x=x1, xend=x1, y=y, yend=y*0.98))+
  geom_segment(data=alpha_posthoc, aes(x=x2, xend=x2, y=y, yend=y*0.98))+
  geom_text(data=alpha_posthoc, aes(x=c(1.5,2,2.5, 1.5,2,2.5,1.5,2,2.5), y=y, label=sig))+
  geom_jitter()+
  facet_wrap(. ~ variable, scales = "free_y", nrow=1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  ylab("Value")+
  xlab("")+
  coord_cartesian(xlim = c(0.4,3.6), expand = c(0,0))+
  scale_x_discrete(labels=c("Field\n(in this study)","Serial cultivation model\nin pot", "Nitrogen source treatment\nafter the serial model in pot"))+
  theme(legend.position = "top",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size=13, angle=45, hjust=1),
        strip.text = element_text(size=13))

alpha_diversity_p

#=================================================
# Relative abundance bar graph (beta diversity)
#=================================================
# sorting relative abundance data  
relative_abundance <- plot_bar(phy_obj_relative)$data
relative_abundance <- relative_abundance[as.character(1:dim(relative_abundance)[1]),]

# output relative abandance data and NA replacement
write.table(relative_abundance,  file = "metagenome_analysis_results/relative_abundance.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "Not applicable")
relative_abundance <- read.csv("metagenome_analysis_results/relative_abundance.tsv", sep="\t")

require(tidyverse)
# Sums data by taxa levels
phylum_relative_abundance <- aggregate(Abundance ~ phylum + Sample + Site_ID + cultivation + Cycling + soil.management + years + design, relative_abundance, sum)
class_relative_abundance <- aggregate(Abundance ~ class + Sample + Site_ID + cultivation + Cycling + soil.management + years + design, relative_abundance, sum)
order_relative_abundance <- aggregate(Abundance ~ order + Sample + Site_ID + cultivation + Cycling + soil.management + years + design, relative_abundance, sum)
family_relative_abundance <- aggregate(Abundance ~ family + Sample + Site_ID + cultivation + Cycling + soil.management + years + design, relative_abundance, sum)
genus_relative_abundance <- aggregate(Abundance ~ genus + Sample + Site_ID + cultivation + Cycling + soil.management + years + design, relative_abundance, sum)

# top 10 sorting
# phylum_level
sum_phylum <- aggregate(Abundance ~ phylum, subset(relative_abundance, cultivation=="field"), sum)
sum_phylum <- sum_phylum %>% 
  arrange(desc(Abundance))
top_phylum <- as.vector(sum_phylum[1:10,1])

# class_level
sum_class <- aggregate(Abundance ~ class, subset(relative_abundance, cultivation=="field"), sum)
sum_class <- sum_class %>% 
  arrange(desc(Abundance))
top_class <- as.vector(sum_class[1:10,1])

# order_level
sum_order <- aggregate(Abundance ~ order, subset(relative_abundance, cultivation=="field"), sum)
sum_order <- sum_order %>% 
  arrange(desc(Abundance))
top_order <- as.vector(sum_order[1:10,1])


# family_level
sum_family <- aggregate(Abundance ~ family, subset(relative_abundance, cultivation=="field"), sum)
sum_family <- sum_family %>% 
  arrange(desc(Abundance))
top_family <- as.vector(sum_family[1:10,1])


# genus_level
sum_genus <- aggregate(Abundance ~ genus, subset(relative_abundance, cultivation=="field"), sum)
sum_genus <- sum_genus %>% 
  arrange(desc(Abundance))
top_genus <- as.vector(sum_genus[1:10,1])

require(ggpubr)
# taxa top bar graph

top_bar <- function(df, fill_label="Top taxa", top_toxa_names){
  df$top_taxa <- df[,1]
  df$top_taxa <- factor(df$top_taxa, levels = c(top_toxa_names,"others"))
  df$top_taxa[is.na(df$top_taxa)] <- "others"
  df$Sample <- factor(df$Sample, levels = sample_ord) 
  df$Sample_group <- df$Cycling
  df$Sample_group <- factor(df$Sample_group, levels = unique(df$Sample_group)[c(1:5, 7:15,6, 16)])
  agg_df <- aggregate(Abundance ~ Sample + Sample_group + Site_ID + cultivation + Cycling + soil.management + years + design + top_taxa, df, sum)
  agg_df$floor <- ifelse(agg_df$Sample_group %in% c("2 years field", "3 years field", "4 years field", "5 years field", "6 years field", "1st", "2nd", "3rd"), "1", "2")
  agg_df$floor[which(agg_df$design == "pot+nitrogen")] <- "3" 
  agg_df$Site_ID <- gsub("NH4Cl", "NH₄Cl", agg_df$Site_ID)
  agg_df$Site_ID <- factor(agg_df$Site_ID, levels = c("control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
  
  p1 <- ggplot(data = subset(agg_df, floor=="1"), mapping = aes(x = Sample, y=Abundance, fill=top_taxa))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c(brewer.pal(n=10, name="Paired"),"#999999"))+
    theme_light()+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill="grey90", color = "white", size=1),
          strip.text = element_text(color="black"))+
    facet_grid(. ~ Sample_group, space = "free", scales = "free")+
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), expand = c(0,0))+ # fill option
    ylab("Relative abundance")+
    xlab("")+
    labs(fill=fill_label)
  p2 <- ggplot(data = subset(agg_df, floor=="2"), mapping = aes(x = Sample, y=Abundance, fill=top_taxa))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c(brewer.pal(n=10, name="Paired"),"#999999"))+
    theme_light()+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill="grey90", color = "white", size=1),
          strip.text = element_text(color="black"))+
    facet_grid(. ~ Sample_group, space = "free", scales = "free")+
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), expand = c(0,0))+ # fill option
    ylab("Relative abundance")+
    xlab("")+
    labs(fill=fill_label)
  p3 <- ggplot(data = subset(agg_df, floor=="3"), mapping = aes(x = Sample, y=Abundance, fill=top_taxa))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c(brewer.pal(n=10, name="Paired"),"#999999"))+
    theme_light()+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill="grey90", color = "white", size=1),
          strip.text = element_text(color="black"))+
    facet_nested(.~Sample_group+Site_ID, space = "free", scales = "free")+
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), expand = c(0,0))+ # fill option
    ylab("Relative abundance")+
    xlab("")+
    labs(fill=fill_label)
  p1
  p2
  p3
  print(ggarrange(p1, p2, p3, nrow=3, legend = "right", common.legend=T, heights = c(1,1,1.1)))
}

top_10_phylum_bar <- top_bar(phylum_relative_abundance, fill_label = "Field top 10 phylum", top_toxa_names = top_phylum )
top_10_phylum_bar 

top_10_class_bar <- top_bar(class_relative_abundance,  fill_label = "Field top 10 class", top_toxa_names = top_class )
top_10_class_bar

top_10_order_bar <- top_bar(order_relative_abundance, fill_label = "Field top 10 order", top_toxa_names = top_order )
top_10_order_bar

top_10_family_bar <- top_bar(family_relative_abundance, fill_label = "Field top 10 family", top_toxa_names = top_family )
top_10_family_bar

top_10_genus_bar <- top_bar(genus_relative_abundance, fill_label = "Field top 10 genus", top_toxa_names = top_genus )
top_10_genus_bar

nitrogene_treatment_OTU <- 
rowSums(
OTU_counts_table[,sample_ord[1:30]]
)

nitrogene_treatment_OTU <- nitrogene_treatment_OTU[nitrogene_treatment_OTU >0]


cycling_treatment_OTU <- 
  rowSums(
    OTU_counts_table[,sample_ord[31:70]]
  )

cycling_treatment_OTU <- cycling_treatment_OTU[cycling_treatment_OTU >0]


field_OTU <- 
  rowSums(
    OTU_counts_table[,sample_ord[71:85]]
  )

field_OTU <- field_OTU[field_OTU > 0]

OTU_venn_list<- list("Serial cultivation" = names(cycling_treatment_OTU),
     "Nitrogen treatment" = names(nitrogene_treatment_OTU),
     "Field"=names(field_OTU))

Venn_OTU_p <- ggvenn(OTU_venn_list, show_percentage = FALSE)

NCF_OTU <- names(nitrogene_treatment_OTU[names(nitrogene_treatment_OTU)%in%names(cycling_treatment_OTU)&names(nitrogene_treatment_OTU)%in%names(field_OTU)])

length(NCF_OTU)

relative_abundance_of_common_OTUs <- data.frame(colSums(t(OTU_rel)[NCF_OTU]))
colnames(relative_abundance_of_common_OTUs) <- "Relative abundance"
relative_abundance_of_common_OTUs$sample <- rownames(relative_abundance_of_common_OTUs)
relative_abundance_of_common_OTUs$groups <- c(rep("Nitrogen treatment", 30), rep("Serial cultivation",40), rep("Field",15))

relative_abundance_of_common_OTUs$groups <- factor(relative_abundance_of_common_OTUs$groups, levels= unique(relative_abundance_of_common_OTUs$groups)[3:1])


relative_abundance_of_common_OTU_p<- 
  ggplot(relative_abundance_of_common_OTUs, aes(x=groups, y=`Relative abundance`*100)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,100))+
  theme(panel.background = element_blank(),
        axis.line = element_line(size=1, color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Relative abundance (%)\n of The 593 common OTUs")+
  xlab("")

family_relative_abundance2 <- dcast(family_relative_abundance[,c("family", "Sample", "Abundance")], family ~ Sample)
rownames(family_relative_abundance2) <- family_relative_abundance2$family
family_relative_abundance2 <- family_relative_abundance2[, -1]


nitrogene_treatment_family <- 
  rowSums(
    family_relative_abundance2[,sample_ord[1:30]]
  )

nitrogene_treatment_family <- nitrogene_treatment_family[nitrogene_treatment_family >0]


cycling_treatment_family <- 
  rowSums(
    family_relative_abundance2 [,sample_ord[31:70]]
  )

cycling_treatment_family <- cycling_treatment_family[cycling_treatment_family >0]


field_family <- 
  rowSums(
    family_relative_abundance2[,sample_ord[71:85]]
  )

field_family <- field_family[field_family > 0]

family_venn_list <- list("Serial cultivation" = names(cycling_treatment_family),
                         "Nitrogen treatment" = names(nitrogene_treatment_family),
                         "Field"=names(field_family))

Venn_family_p <- ggvenn(family_venn_list, show_percentage = FALSE)

NCF_family <- names(nitrogene_treatment_family[names(nitrogene_treatment_family)%in%names(cycling_treatment_family)&names(nitrogene_treatment_family)%in%names(field_family)])

relative_abundance_of_common_family <- data.frame(colSums(family_relative_abundance2[NCF_family,]))
colnames(relative_abundance_of_common_family) <- "Relative abundance"
relative_abundance_of_common_family$sample <- rownames(relative_abundance_of_common_family)
relative_abundance_of_common_family$groups <- c(rep("Nitrogen treatment", 30), rep("Serial cultivation",40), rep("Field",15))

relative_abundance_of_common_family$groups <- factor(relative_abundance_of_common_family$groups, levels= unique(relative_abundance_of_common_family$groups)[3:1])

relative_abundance_of_common_family_p <- 
  ggplot(relative_abundance_of_common_family, aes(x=groups, y=`Relative abundance`*100)) + 
  geom_boxplot() + 
  coord_cartesian(ylim=c(0,100))+
  theme(panel.background = element_blank(),
        axis.line = element_line(size=1, color="black"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Relative abundance (%)\n of The 146 common families")+
  xlab("")


ggarrange(Venn_OTU_p, relative_abundance_of_common_OTU_p, Venn_family_p, relative_abundance_of_common_family_p)

#========
# PCoA
#===========
# Bray-Curtis
OTU_rel <- phy_obj_relative@otu_table

bray <- as.matrix(vegdist(OTU_rel, method="bray"))

pcoa_bray <-cmdscale(bray, eig=TRUE)
colnames(pcoa_bray$points) <- c("PCoA1", "PCoA2")
pcoa_bray$points <- cbind(pcoa_bray$points, SAM)

pcoa_bray$points$group <- c(pcoa_bray$points$Site_ID[1:30], pcoa_bray$points$Cycling[31:70], pcoa_bray$points$cultivation[71:85])
pcoa_bray$points$group <- gsub("NH4Cl", "NH₄Cl", pcoa_bray$points$group)
pcoa_bray$points$group <- gsub("control", "Control", pcoa_bray$points$group)
pcoa_bray$points$group <- gsub("field", "Field", pcoa_bray$points$group)
pcoa_bray$points$group <- factor(pcoa_bray$points$group, levels=c("Serial cultivation in pot", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th", " ", "Nitrogen treatment in pot", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val","  ", "In this study","Field"))

group_lv <- levels(pcoa_bray$points$group)
group_color <- hcl(seq(15,325,length.out=17), c=100, l=65, alpha=1)

PCoA_bray_p <- 
  ggplot(pcoa_bray$point, aes(PCoA1, PCoA2, color=group, fill=group, group=group))+
  geom_encircle(aes(fill=group, group=group), expand = 0, s_shape=1, alpha =0.1)+
  geom_point(size=3, aes(shape=group))+
  labs(x=paste0("PCoA 1 (",
                round(pcoa_bray$eig[1]/sum(pcoa_bray$eig), 3)*100,
                "%)"),
       y=paste0("PCoA 2 (",
                round(pcoa_bray$eig[2]/sum(pcoa_bray$eig), 3)*100,
                "%)"),
       color="",
       fill="",
       shape="",
       alpha=""
  )+
  scale_fill_manual(values=c("white", group_color[1:10],
                             "white","white", group_color[11:16],
                             "white","white", group_color[17]),
                    drop=FALSE)+
  scale_color_manual(values=c("white", group_color[1:10],
                              "white","white", group_color[11:16],
                              "white","white", group_color[17]),
                     drop=FALSE)+
  scale_shape_manual(values=c(NA, 1:10, NA, NA, 11:16, NA, NA, 17), drop=FALSE)+
  scale_alpha_manual(values=c(0, rep(1, 10),
                              0, 0, rep(1, 6),
                              0, 0, 1),
                     drop=FALSE)+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  guides(fill=guide_legend(ncol=1), color=guide_legend(ncol=1), shape=guide_legend(ncol=1))+
  coord_fixed()+
  theme_light()
  
PCoA_bray_p

#PICRUSt2 PW
PW <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/pathways_out/path_abun_unstrat.tsv.gz"), sep="\t", row.names = 1)

PW_bray <- as.matrix(vegdist(t(PW), method="bray"))
PW_pcoa_bray <-cmdscale(PW_bray, eig=TRUE)
colnames(PW_pcoa_bray$points) <- c("PCoA1", "PCoA2")
PW_pcoa_bray$point <- data.frame(PW_pcoa_bray$point)

PW_pcoa_bray$points <- cbind(PW_pcoa_bray$points, SAM)
PW_pcoa_bray$points$group <- c(PW_pcoa_bray$points$Site_ID[1:30], PW_pcoa_bray$points$Cycling[31:70], PW_pcoa_bray$points$cultivation[71:85])
PW_pcoa_bray$points$group <- gsub("NH4Cl", "NH₄Cl", PW_pcoa_bray$points$group)
PW_pcoa_bray$points$group <- gsub("control", "Control", PW_pcoa_bray$points$group)
PW_pcoa_bray$points$group <- gsub("field", "Field", PW_pcoa_bray$points$group)
PW_pcoa_bray$points$group <- factor(PW_pcoa_bray$points$group, levels=c("Serial cultivation in pot", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th", " ", "Nitrogen treatment in pot", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val","  ", "In this study","Field"))

PW_group_lv <- levels(PW_pcoa_bray$points$group)
PW_group_color <- hcl(seq(15,325,length.out=17), c=100, l=65, alpha=1)

PW_PCoA_bray_p <- 
  ggplot(PW_pcoa_bray$points, aes(PCoA1, PCoA2, color=group, fill=group, group=group))+
  geom_encircle(expand = 0, s_shape=1, alpha =0.2)+
  geom_point(size=3, aes(shape=group)) +
  labs(x=paste0("PCoA 1 (",
                round(PW_pcoa_bray$eig[1]/sum(PW_pcoa_bray$eig), 3)*100,
                "%)"),
       y=paste0("PCoA 2 (",
                round(PW_pcoa_bray$eig[2]/sum(PW_pcoa_bray$eig), 3)*100,
                "%)"),
       color="",
       fill="",
       shape="",
       alpha=""
  )+
  scale_fill_manual(values=c("white", PW_group_color[1:10],
                            "white","white", PW_group_color[11:16],
                            "white","white", PW_group_color[17]),
                    drop=FALSE)+
  scale_color_manual(values=c("white", PW_group_color[1:10],
                             "white","white", PW_group_color[11:16],
                             "white","white", PW_group_color[17]),
                     drop=FALSE)+
  scale_shape_manual(values=c(NA, 1:10, NA, NA, 11:16, NA, NA, 17), drop=FALSE)+
  scale_alpha_manual(values=c(0, rep(1, 10),
                              0, 0, rep(1, 6),
                              0, 0, 1),
                     drop=FALSE)+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  guides(fill=guide_legend(ncol=1), color=guide_legend(ncol=1), shape=guide_legend(ncol=1))+
  coord_fixed()+
  theme_light()
PW_PCoA_bray_p


PW_permanova <- adonis(t(PW_norm.data) ~ group, PW_pcoa_bray$points, permutations=999, method="bray")
PW_permanova

library(pairwiseAdonis)
PW_pwadonis <- pairwise.adonis(PW_bray, PW_pcoa_bray$points[,"group"], p.adjust.m="fdr")
PW_pwadonis$sig <- ifelse(PW_pwadonis$p.adjusted > 0.05, "", 
       ifelse(PW_pwadonis$p.adjusted > 0.01, "*", 
              ifelse(PW_pwadonis$p.adjusted > 0.001, "**", "***" )))

PW_pwadonis[str_detect(PW_pwadonis$pairs, "Field"),]
colnames(PW_pwadonis) <- c("Compair", "Df", "Sum of Sqs", "F. Model", "R2", "*P*", "*P~adj~*", "significant")
require("ftExtra")
#===================
# Phylogenetic tree
#===================
# Load fast file
fasta_path <- "metagenome_analysis_results/DADA2_results/OTUs_sequence.fasta"
seqs <- readDNAStringSet(fasta_path)
seqs <- OrientNucleotides(seqs)

# Alignment
aligned_seqs <- AlignSeqs(seqs)

fasta_path <- "metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs_alignes_sequence.fasta"
writeXStringSet(aligned_seqs, file = fasta_path)

#===================================================================
# bNTI and beta RC
#==================================================================
rownames(OTU_counts_table) <- OTU_counts_table$`#OTU ID`

#------------ All
phylogeny <- read.tree("raxml_results/all3.raxml.bestTree")
tip.label <- phylogeny$tip.label
rownames(OTU_counts_table) <- OTU_counts_table$`#OTU ID`
match.phylo.OTU <- match.phylo.data(phylogeny, OTU_counts_table[tip.label, -1])
# calculate empirical betaMNTD
beta.mntd.weighted <- as.matrix(comdistnt(t(match.phylo.OTU$data), cophenetic(match.phylo.OTU$phy), abundance.weighted=T))
dir.create("metagenome_analysis_results/RCI and bNTI")
write.csv(beta.mntd.weighted, "metagenome_analysis_results/RCI and bNTI/betaMNTD_weighted.csv",quote=F);

identical(colnames(match.phylo.OTU$data),colnames(beta.mntd.weighted)) # just a check, should be TRUE

# calculate randomized betaMNTD
reps <- 999; # number of randomizations

rand.weighted.bMNTD.comp <- array(c(-999),dim=c(ncol(match.phylo.OTU$data),ncol(match.phylo.OTU$data), reps))

pb <- progress_bar$new(total = reps)
for (i in 1:reps) {
  rand.weighted.bMNTD.comp[,,i] <- as.matrix(comdistnt(t(match.phylo.OTU$data),taxaShuffle(cophenetic(match.phylo.OTU$phy)),abundance.weighted=T,exclude.conspecifics = F))
  pb$tick()
}

weighted.bNTI <- matrix(c(NA),nrow=ncol(match.phylo.OTU$data),ncol=ncol(match.phylo.OTU$data))

for (columns in 1:(ncol(match.phylo.OTU$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.OTU$data)) {
    rand.vals <- rand.weighted.bMNTD.comp[rows,columns,]
    weighted.bNTI[rows,columns] <- (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
    rm(rand.vals)
  }
}
weighted.bNTI_mat <- as.matrix(as.dist(weighted.bNTI))
colnames(weighted.bNTI_mat) <- sample_ord
rownames(weighted.bNTI_mat) <- sample_ord

for(i in 1:length(sample_ord)){
  weighted.bNTI_mat[i,i] <- NA
}

weighted.bNTI_ldf <- melt(weighted.bNTI_mat)

colnames(weighted.bNTI_ldf) <- c("Sample1" ,"Sample2", "weighted_bNTI")
weighted.bNTI_ldf$Sample1 <- factor(weighted.bNTI_ldf$Sample1, levels = sample_ord)
weighted.bNTI_ldf$Sample2 <- factor(weighted.bNTI_ldf$Sample2, levels = sample_ord)

weighted.bNTI_ldf$weighted_bNTI2 <- ifelse(weighted.bNTI_ldf$weighted_bNTI > 2, "betaNTI > 2",
                                                      ifelse(weighted.bNTI_ldf$weighted_bNTI > -2, "2 > betaNTI > -2", "-2 > betaNTI"))

weighted.bNTI_ldf$weighted_bNTI2 <- factor(weighted.bNTI_ldf$weighted_bNTI2 , levels=c("betaNTI > 2", "2 > betaNTI > -2", "-2 > betaNTI"))

betaNTI_feid <- ggplot(weighted.bNTI_ldf, aes(x=Sample1,y=Sample2))+
  geom_tile(aes(fill = weighted_bNTI))+
  scale_fill_gradient2(low = "steelblue", mid="grey80", high = "tomato", breaks=-3:3*2, na.value = "black")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=c(1:6*5, (7:16*4+6))+0.5, yend=c(1:6*5, (7:16*4+6))+0.5, size=0.2)+
  annotate(geom="segment", x=c(1:6*5, (7:16*4+6))+0.5, xend=c(1:6*5, (7:16*4+6))+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste(Rhizosphere~beta,NTI)))+

coord_fixed()

betaNTI_feid2 <- ggplot(weighted.bNTI_ldf, aes(x=Sample1,y=Sample2))+
  geom_tile(aes(fill = weighted_bNTI2))+
  scale_fill_manual(values=c("#FF866C","grey80","#1752FF"), na.value="black",
                    labels = c(expression(paste(beta,NTI>2)),
                               expression(paste(2>beta,NTI)>-2),
                               expression(paste(-2>beta,NTI)),
                               expression(self~pairwise)
                    ))+
  annotate(geom="segment", x=Inf, xend=-Inf, y=c(1:6*5, (7:16*4+6))+0.5, yend=c(1:6*5, (7:16*4+6))+0.5, size=0.2)+
  annotate(geom="segment", x=c(1:6*5, (7:16*4+6))+0.5, xend=c(1:6*5, (7:16*4+6))+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  coord_fixed()+
  labs(fill=expression(paste(Rhizosphere~beta,NTI)))

ggarrange(betaNTI_feid, betaNTI_feid2)

weighted.bNTI_ldf$Sample1
stochastic_betaD_similarity <- subset(weighted.bNTI_ldf, Sample1 %in% sample_ord[71:85])
stochastic_betaD_similarity <- stochastic_betaD_similarity[!(stochastic_betaD_similarity$Sample1 == stochastic_betaD_similarity$Sample2), ]
stochastic_betaD_similarity <- stochastic_betaD_similarity[stochastic_betaD_similarity$weighted_bNTI2 == "-2 > betaNTI", ]

stochastic_betaD_similarity$Sample1_char <- "Field"
stochastic_betaD_similarity$Sample2_char <- gsub("^R", "", stochastic_betaD_similarity$Sample2)
stochastic_betaD_similarity$Sample2_char <- gsub("^[0-9]-[0-9]$", "Field", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("[0-9]$", "", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("Con", "Control", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("NH4Cl", "NH₄Cl", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^A", "", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^sn", "Asn", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^sp", "Asp", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("C$", "", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^1$", "1st", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^2$", "2nd", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^3$", "3rd", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^4$", "4th", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^5$", "5th", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^6$", "6th", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^7$", "7th", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^8$", "8th", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("^9$", "9th", stochastic_betaD_similarity$Sample2_char)
stochastic_betaD_similarity$Sample2_char <- gsub("10$", "10th", stochastic_betaD_similarity$Sample2_char)

stochastic_betaD_similarity2 <- stochastic_betaD_similarity %>%
  group_by(Sample1_char, Sample2_char)%>%
  count(weighted_bNTI2)

stochastic_betaD_similarity2$total_paired_n <- 15*c(rep(4,10), rep(5, 3), 14, rep(5, 3))

stochastic_betaD_similarity2$similarity <- stochastic_betaD_similarity2$n / stochastic_betaD_similarity2$total_paired_n *100
stochastic_betaD_similarity2$Sample2_char <- factor(stochastic_betaD_similarity2$Sample2_char, levels = c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th", 
                                                                                             "8th", "9th", "10th", "Control", "NH₄Cl", "Glu",
                                                                                             "Asn", "Asp", "Val", "Field"))
stochastic_betaD_similarity2%>%
  ggplot(aes(x=Sample2_char, y=similarity))+
  geom_bar(stat='identity')+
  ylab("Phylogenetic weighted betadiviersity similarity (%)")

SAM$group2 <- ifelse(SAM$group %in% c("NH4Cl", "Field"), "Feild and NH₄Cl treatment", "Others")

require(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = round(PW),
                              colData = SAM,
                              design = ~ group2) 
dds <- DESeq(dds)
res <- results(dds, pAdjustMethod = "fdr")
res <- data.frame(res)
res$sig <- "Nothing"
res$sig[res$padj < 0.01&res$log2FoldChange > 0.5] <- "Less in field and NH₄Cl treatment"
res$sig[res$padj < 0.01&res$log2FoldChange < -0.5] <- "More in field and NH₄Cl treatment"
res$description <- PW_d[rownames(res), "description"]
important_PW <- c("AST-PWY", "ARGORNPROST-PWY" ,"DHGLUCONATE-PYR-CAT-PWY", "PWY0-1338", "PYRIDOXSYN-PWY", "PWY-5910", "PWY-7391", "PWY-922")
res$description <- gsub("geranylgeranyldiphosphate biosynthesis I [(]", "geranylgeranyldiphosphate biosynthesis I\n(", res$description)
res[important_PW,]

PW_volcano <-
  ggplot(data.frame(res), aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(color=sig), alpha=0.5)+
  geom_text_repel(data=res[important_PW,], aes(label=description), min.segment.length = 0.01, nudge_x = 3, nudge_y=3)+
  scale_color_manual(values=c("tomato", "steelblue","grey80"))+
  labs(color="")+
  ylab(expression(-log[10](italic(P[adj]))))+
  xlab(expression(log[2](FC)))+
  theme_light()+
  theme(axis.line = element_line(color="black"))

PW_volcano

sig_res <- subset(res, padj < 0.01&abs(log2FoldChange)>0.5)

Less_in_field_NH4Cl <- sig_res[sig_res$log2FoldChange>0,]
More_in_field_NH4Cl <- sig_res[sig_res$log2FoldChange<0,]

PW_d <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/pathways_out/path_abun_unstrat_descrip.tsv.gz"), sep="\t", row.names = 1)
PW_d
Less_in_field_NH4Cl <- Less_in_field_NH4Cl[order(Less_in_field_NH4Cl$log2FoldChange, decreasing = T),]
More_in_field_NH4Cl <- More_in_field_NH4Cl[order(More_in_field_NH4Cl$log2FoldChange),]

Less_in_field_NH4Cl$description <- PW_d[rownames(Less_in_field_NH4Cl),]$description
More_in_field_NH4Cl$description <- PW_d[rownames(More_in_field_NH4Cl),]$description
Less_in_field_NH4Cl
More_in_field_NH4Cl

MetaCYC_ontology <- read.csv("MetaCYC_ontology.tsv", sep="\t", row.names = 1)

Less_in_field_NH4Cl_ontology <- MetaCYC_ontology[rownames(Less_in_field_NH4Cl),]
More_in_field_NH4Cl_ontology <- MetaCYC_ontology[rownames(More_in_field_NH4Cl),]

Less_in_field_NH4Cl$ontology <- Less_in_field_NH4Cl_ontology$Ontology_parents_of_class_lv.5
More_in_field_NH4Cl$ontology <- More_in_field_NH4Cl_ontology$Ontology_parents_of_class_lv.5

Less_in_field_NH4Cl$group <- "Less in field NH₄Cl"
More_in_field_NH4Cl$group <- "More in field NH₄Cl"

PW_ontology <- rbind(Less_in_field_NH4Cl, More_in_field_NH4Cl)

PW_onto_count <- c()
PW_group <- c()
PW_onto <- c()
for(i in 1:length(unique(PW_ontology$ontology))){
  PW_onto <- c(PW_onto, unique(PW_ontology$ontology)[i])
  PW_onto_count <- c(PW_onto_count,  
                   nrow(
                     PW_ontology[which(PW_ontology$group == "Less in field NH₄Cl" & PW_ontology$ontology==unique(PW_ontology$ontology)[i]), ]
                   )
)
  PW_group <- c(PW_group, "Less in field NH₄Cl")
  
  PW_onto <- c(PW_onto, unique(PW_ontology$ontology)[i])
  PW_onto_count <- c(PW_onto_count,  
                     nrow(
                       PW_ontology[which(PW_ontology$group == "More in field NH₄Cl" & PW_ontology$ontology==unique(PW_ontology$ontology)[i]), ]
                     )
  )
  PW_group <- c(PW_group, "More in field NH₄Cl")

}

PW_ontology2 <- data.frame(PW_group,PW_onto, PW_onto_count)
PW_ontology2$PW_onto[is.na(PW_ontology2$PW_onto)]<-""
PW_ontology2 <- aggregate(PW_onto_count ~ PW_group +PW_onto,  PW_ontology2, sum)
PW_ontology2$PW_onto <- gsub("^$", "Not assigned", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("CARBO-BIOSYNTHESIS", "Carbohydrate Biosynthesis", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("AROMATIC-COMPOUNDS-DEGRADATION", "Aromatic-Compound-Ddegradation", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("IND-AMINO-ACID-SYN", "Individual Amino Acids Biosynthesis", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("MANDELATE-DEG", "Mandelate Degradation", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("C1-COMPOUNDS", "C1 Compound Utilization and Assimilation", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("PYR-NUC-SYN", "Pyrimidine Nucleotide Biosynthesis", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("GALLATE-DEG", "Gallate Degradation", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("AMINE-DEG", "Amine and Polyamine Degradation", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("SECONDARY METABOLITE DEGRADATION", "Secondary Metabolite Degradation", PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- gsub("-", " ", PW_ontology2$PW_onto)

unique(PW_ontology2$PW_onto)
PW_ontology2$PW_onto <- factor(PW_ontology2$PW_onto, levels=unique(PW_ontology2$PW_onto)[
  order(PW_ontology2$PW_onto_count[1:54*2-1] - PW_ontology2$PW_onto_count[1:54*2], decreasing = T)]
)

PW_ontology2$PW_onto <- factor(PW_ontology2$PW_onto, levels = levels(PW_ontology2$PW_onto)[c(1:16,18:53,17)])

PW_onthology <-
ggplot(data = PW_ontology2, aes(x=PW_onto, y= ifelse(PW_group=="More in field NH₄Cl", 1,-1)*PW_onto_count, fill = PW_group)) +
  geom_bar(stat = "identity")+
  scale_fill_discrete(labels=c("Less in field and NH₄Cl", "More in field and NH₄Cl"))+
  labs(fill="")+
  ylab("Count of MetaCyc ontology class")+
  xlab("MetaCyc ontology class")+
  scale_y_continuous(breaks = seq(-4,6,by=2), labels=abs(seq(-4,6,by=2)))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
        axis.line = element_line(color="black"),
        panel.background = element_blank(),
        legend.position = "top")
#=======

family_relative_abundance2 <- dcast(family_relative_abundance[,c("family", "Sample", "Abundance")], family ~ Sample)
family_relative_abundance2[,-1]<- family_relative_abundance2[,-1]*1e9
rownames(family_relative_abundance2) <- family_relative_abundance2[,1]
dds_f <- DESeqDataSetFromMatrix(countData = round(family_relative_abundance2[, rownames(SAM)]),
                                colData = SAM,
                                design = ~ group2) 
dds_f <- DESeq(dds_f)
res_f <- results(dds_f, pAdjustMethod = "fdr")
res_f <- data.frame(res_f)
res_f$baseMean <- res_f$baseMean/1e9
res_f$padj < 0.01 & abs(res_f$log2FoldChange)>1

res_f <- data.frame(res_f)
res_f$sig <- "Nothing"
res_f$sig[res_f$padj < 0.01&res_f$log2FoldChange>1] <- "Less in field and NH₄Cl treatment"
res_f$sig[res_f$padj < 0.01&res_f$log2FoldChange < -1] <- "More in field and NH₄Cl treatment"

require(ggrepel)

important_family <- res_f[which(rownames(res_f) %in% c("Pseudomonadaceae", "Beijerinckiaceae", "Devosiaceae", "Xanthomonadaceae", "Rhizobiaceae", "Sphingobacteriaceae", "Sphingomonadaceae", "Comamonadaceae")), ]
important_family$family <- rownames(important_family)
family_volcano <- 
  ggplot(res_f, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(color=sig), alpha=0.5)+
  geom_text_repel(data=important_family, aes(label=family), min.segment.length = 0.01, nudge_x = 2, nudge_y=2)+
  scale_color_manual(values=c("tomato", "steelblue","grey80"))+
  ylab(expression(-log[10](italic(P[adj]))))+
  xlab(expression(log[2](FC)))+
  theme_light()+
  theme(axis.line = element_line(color="black"))

ggarrange(PW_volcano, family_volcano, common.legend = T)
