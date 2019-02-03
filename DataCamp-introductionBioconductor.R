#CHAPTER 1 - What is bioconductor?

# Load the BiocInstaller package
library(BiocInstaller)

# Explicit syntax to check the Bioconductor version
BiocInstaller::biocVersion()

# When BiocInstaller is loaded use biocVersion alone
biocVersion()

# Load the BSgenome package
library(BSgenome)

# Check the version of the BSgenome package
packageVersion("BSgenome")

# Investigate about the a_genome using show()
show(a_genome)

# Investigate some other accesors

organism(a_genome)
provider(a_genome)
seqinfo(a_genome)

# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

# Select chromosome M, alias chrM
getSeq(yeastGenome, "chrM")

# Count characters of the chrM sequence
nchar(yeastGenome$chrM)
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the first 30 bases of each chromosome
getSeq(yeastGenome,start = 1, end = 30)

#########################################################
#CHAPTER 2 - Biostrings

# Load packages
library(Biostrings)

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

# Check the alphabetFrequency of the zikaVirus
alphabetFrequency(zikaVirus)

# Check alphabet of the zikaVirus using baseOnly = TRUE
alphabet(zikaVirus, baseOnly = TRUE)
# Unlist the set and select the first 21 letters as dna_seq, then print it
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# 1.1 Transcribe dna_seq as rna_seq, then print it
rna_seq <- RNAString(dna_seq) 
rna_seq

# 1.2 Translate rna_seq as aa_seq, then print it
aa_seq <- translate(rna_seq)
aa_seq

# 2.1 Translate dna_seq as aa_seq_2, then print it
aa_seq_2 <- translate(dna_seq)
aa_seq == aa_seq_2
# Create zikv with one collated sequence using `zikaVirus`
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Check the width of zikaVirus
width(zikaVirus)

# Subset zikv to only the first 30 bases
subZikv <- subseq(zikv, end = 30)
subZikv
length(subZikv)
# The reverse of zikv is
reverse(zikv)

# The complement of zikv is
complement(zikv)

# The reverse complement of zikv is
reverseComplement(zikv)

# The translation of zikv is
translate(zikv)
# Find palindromes in zikv
findPalindromes(zikv)
# print the rnaframesZikaSet 
rnaframesZikaSet

# translate all 6 reading frames 
AAzika6F <- translate(rnaframesZikaSet)
AAzika6F

# Count the matches allowing 15 mistmatches
vcountPattern(pattern = ns5, subject = AAzika6F, max.mismatch = 15)

# Select the frame that contains the match

selectedSet <- AAzika6F[3]

#Convert this frame into a single sequence
selectedSeq <- unlist(selectedSet)

selectedSet
selectedSeq
# Use vmatchPattern with the set
vmatchPattern(pattern = ns5, subject = selectedSet, max.mismatch = 15)

# Use matchPattern with the single sequence
matchPattern(ns5, selectedSeq, max.mismatch = 15)

#####################################################################
#CHAPTER 3 - IRanges and GenomicRanges

# load package IRanges
library(IRanges)

# start vector 1 to 5 and end 100 
IRnum1 <- IRanges(start = c(1,2,3,4,5), end = 100)

# end 100 and width 89 and 10
IRnum2 <- IRanges(end = 100, width = c(89, 10))

# logical argument start = Rle(c(F, T, T, T, F, T, T, T))
IRlog1 <- IRanges(start = Rle(c(F,T,T,T,F,T,T,T)))

# Printing objects in a list
print(list(IRnum1 = IRnum1, IRnum2 = IRnum2, IRlog1 = IRlog1))

library(GenomicRanges)

print(seq_intervals)
myGR = GRanges(seq_intervals)

# Load Package Genomic Ranges
library(GenomicRanges)

# Print the GRanges object myGR
myGR

# Check the metadata, if any
mcols(myGR)
# load human reference genome hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# assign hg38 to hg, then print it
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg

# filter 1: extract all the genes in chromosome X as hg_chrXg, then print it
hg_chrXg <- genes(hg, filter = list(tx_chrom = c("chrX")))
hg_chrXg

# filter 2: extract all positive stranded genes in chromosome X as hg_chrXgp, then sort it
hg_chrXgp <- genes(hg, filter = list(tx_chrom = c("chrX"), tx_strand = "+"))
sort(hg_chrXgp)
length(hg_chrXgp)
length(hg_chrXg)

# Store the overlapping range in rangefound
rangefound <- subsetByOverlaps(hg_chrX, ABCD1)

# Check names of rangefound
names(rangefound)

# Check the geneOfInterest 
ABCD1

# Check rangefound
rangefound

# load the human transcripts DB to hg
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene

# prefilter chromosome X "chrX" using seqlevels()
seqlevels(hg) <- c("chrX")

# get all transcripts by gene and print it
hg_chrXt <- transcriptsBy(hg, by = "gene")
hg_chrXt

# select gene `215` from the transcripts
hg_chrXt$`215`

####################################################################
#CHAPTER 3 - ShortRead

# load ShortRead
library(ShortRead)

# print fqsample
fqsample

# class of fqsample
class(fqsample)

# class sread fqsample
class(sread(fqsample))

# id fqsample
id(fqsample)

# load ShortRead
library(ShortRead)

# set a seed for sampling
set.seed(1234)

# Use FastqSampler with f and select 100 reads
fs <- FastqSampler(con = f, n = 100)

# new sample yield
my_sample <- yield(fs)

# print my_sample
my_sample

# load ShortRead
library(ShortRead)

# Check quality
quality(fqsample)

# Check encoding
encoding(quality(fqsample))

# Check baseQuality
qaSummary[["baseQuality"]]

# glimpse nucByCycle
glimpse(nucByCycle)

# make an awesome plot!
nucByCycle %>% 
  # gather the nucleotide letters in alphabet and get a new count column
  gather(key = alphabet, value = count , -cycle) %>% 
  ggplot(aes(x = cycle, y =  count, color = alphabet)) +
  geom_line(size = 0.5 ) +
  labs(y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

# Load package ShortRead
library(ShortRead)

# Check class of fqsample
class(fqsample)

# filter reads into selectedReads using myStartFilter
selectedReads <- fqsample[myStartFilter(fqsample)]

# Check class of selectedReads
class(selectedReads)

# Check detail of selectedReads
detail(selectedReads)

# Load package Rqc
library(Rqc)

# Average per cycle quality plot
rqcCycleAverageQualityPlot(qa)

# Average per cycle quality plot with white background
rqcCycleAverageQualityPlot(qa) + theme_minimal()

# Read quality plot with white background
class(qa)
rqcReadQualityPlot(qa) + theme_minimal()
print(myGR)
