##Adapted from: https://benjjneb.github.io/dada2/tutorial.html
##and https://astrobiomike.github.io/amplicon/dada2_workflow_ex#dada2

##Remove primers sequences using 'cutadapt' (FWD: 17NTs; REV: 21NTs):
##FWD <- "CCTACGGGNGGCWGCAG"
## Reverse complement of Forward primer: "CTGCWGCCNCCCGTAGG"
##REV <- "GACTACHVGGGTATCTAATCC"
## Reverse complement of Reverse primer: "GGATTAGATACCCBDGTAGTC"

for sample in $(cat samples)
do

    echo "On sample: $sample"
    
    cutadapt -a ^CCTACGGGNGGCWGCAG...GGATTAGATACCCBDGTAGTC \
    -A ^GACTACHVGGGTATCTAATCC...CTGCWGCCNCCCGTAGG \
    -o ${sample}_R1_trimmed.fq -p ${sample}_R2_trimmed.fq \
    ${sample}_R1_001.fastq ${sample}_R2_001.fastq \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done

##Have a look at what fraction of reads were retained in each sample (column 2) and what fraction of bps were retained in each sample (column 3):

paste samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")

##We would expect to lose around 13-14% of bps just for cutting off the primers, and the remainder of lost bps would be from the relatively low percent of those reads totally removed (~92-97% across the samples), which could happen for several reasons.

#With primers removed, we’re now ready to switch R and start using DADA2!

R

library(dada2)
library("vegan")
library("ggplot2")
library("tidyr")
library(Biostrings)


setwd("/home/path/to/working/directory")

list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("samples", what="character")
## one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_R1_trimmed.fq")
## and one with the reverse
reverse_reads <- paste0(samples, "_R2_trimmed.fq")

## and variables holding file names for the forward and reverse
## filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_R1_filtered.fq")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq")

##Quality trimming/filtering
##Start by visualizing the trimmed reads
pdf("qualPlots_fwd.pdf")
plotQualityProfile(forward_reads[7:10])
dev.off()

pdf("qualPlots_rev.pdf")
plotQualityProfile(reverse_reads[17:20])
dev.off()

##The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed. Generally speaking, you want the observed (black dots) to track well with the estimated (black line).

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                reverse_reads, filtered_reverse_reads, maxEE=c(2,4),
                rm.phix=TRUE, minLen=175, truncLen=c(280,240), multithread=TRUE, compress = TRUE)

class(filtered_out)
dim(filtered_out)

##Visualize error plot of trimmed reads

pdf("qualPlots_filt_fwd.pdf")
plotQualityProfile(filtered_forward_reads[7:10])
dev.off()

pdf("qualPlots_filt_rev.pdf")
plotQualityProfile(filtered_reverse_reads[17:20])
dev.off()

##Dereplication
## When DADA2 dereplicates sequences, it also generates a new quality-score profile of each unique sequence based on the average quality scores of each base of all of the sequences that were replicates of it.

derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

##Inferring ASVs

##Dada2 infers true biological sequences by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence is more likely to be of biological origin or more likely to be spurious.

##Here we use the pseudo-pooling option - pooling allows information to be shared across samples, which makes it easier to resolve rare variants that were seen just once or twice in one sample but many times across samples.By default, all ASVs detected in at least two samples in the first sample processing step are input as priors to the second step, but that condition can be changed with the PSEUDO_PREVALENCE and PSEUDO_ABUNDANCE dada options.

dada_forward <- dada(derep_forward, err=err_forward_reads, pool = "pseudo", multithread=TRUE) 
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool = "pseudo", multithread=TRUE) 

##Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                    derep_reverse, verbose = TRUE)


## this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 170 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

##Generating a count table (aka OTU matrix)
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix

##Chimera identification
##DADA2 identifies likely chimeras by aligning each sequence with those that were recovered in greater abundance and then seeing if there are any lower-abundance sequences that can be made exactly by mixing left and right portions of two of the more-abundant ones. These are then removed.

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) 
## though we only lost some sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that:

sum(seqtab.nochim)/sum(seqtab) # 0.917476 # good, we barely lost any in terms of abundance

##Overview of counts throughout
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
               filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
               dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
               nonchim=rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

save.image("rev_opt_pool_orig.RData")

##Assigning taxonomy
library(DECIPHER)

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying

##Using the naive Bayesian classifier method and RDP database:
taxa<- assignTaxonomy(seqtab.nochim, "/home/path/to/rdp_tax_files/rdp_train_set_16.fa.gz", multithread=TRUE, minBoot=50)
taxa <- addSpecies(taxa, "/home/path/to/rdp_tax_files/rdp_species_assignment_16.fa.gz")

save.image("dada_output.RData") #for downstream analyses.








































































































































































































































































































































