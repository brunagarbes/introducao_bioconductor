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


