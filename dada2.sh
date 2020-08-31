
##To save workspace and all functions created:
#Option1:
save(file = “d:/filename.RData”)
#Option2
save.image(“d:/filename.RData”)

library(rpart)
library(ips) #Raxml function for tree construction
library(ShortRead)
library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(knitr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(ggpubr)
library(dplyr)
library(btools)
library(MicrobiomeX)
library(aMiAD)
library(picante)
library(entropart)
library(microbiomeSeq)
library(DESeq2)
library(metagenomeSeq)
library(apeglm)
library(ashr)
library(btools)
library(GUniFrac) # UniFracs
library(ade4) # s.class
library(phytools) # Load tree - read.newick
library(cluster) # Clustering
library(clusterSim) # Clustering
library(ggplot2) # Graphs
library(car) # Anova
library(NMF) # Heatmap
library(qvalue) # False discovery rate
library(cowplot) # plot_grid
library(grid) #manual plot annotation
library(ggdendro) # Dendrogram
library(BiodiversityR) # Alpha diversity
library(reshape2) # melt
#library(remotes) # for intalling github packages (useful for Tax4Fun installation)
library(Tax4Fun) # Metagenome prediction
library(ecodist) # Mantel test
library(plot3D) # 3D graphs
library(phangorn) #Phylogenetic tree construction


##===========Dada2 Pre-processing=======================
###This follows the tutorial on https://benjjneb.github.io/dada2/tutorial.html and https://benjjneb.github.io/dada2/ITS_workflow.html (June 2019).

setwd("/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/")
path<-"/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2"
list.files(path)

## Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

## Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##Inspect read quality profiles
##We start by visualizing the quality profiles of the forward reads:

plotQualityProfile(fnFs[1:2])

##In gray-scale is a heat map of the frequency of each quality score at each base position. The median quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position (this is more useful for other sequencing technologies, as Illumina reads are typically all the same lenghth, hence the flat red line).

##The forward reads are good quality. We generally advise trimming the last few nucleotides to avoid less well-controlled errors that can arise there. These quality profiles do not suggest that any additional trimming is needed. We will truncate the forward reads at position 240 (trimming the last 10 nucleotides).

##Now we visualize the quality profile of the reverse reads:

plotQualityProfile(fnRs[1:2])

##Remove primer sequences (FWD: 17NTs; REV: 21NTs): Not really applicable to 16S. Primers accounted for with trunlen option in the FilterandTrim command.
##FWD <- "CCTACGGGNGGCWGCAG"
##REV <- "GACTACHVGGGTATCTAATCC"

##The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.

##Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

##We’ll use suggested modified* (standard) filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2 for the fwd and 4* for the reverse strand. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. The trimleft option contains the length of the fwd and rev primers

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,240), maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE, trimLeft=c(17, 21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

##Learn the Error Rates
##The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates. The learnErrors method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

##It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Dereplicate identical reads: how many unique sequences contained in the total reads in a sample?

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

##Name the derep-class objects by the sample names

names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Sample Inference: how many true ASVs are contained in the unique sequence reads?

##We are now ready to apply the core sample inference algorithm to the dereplicated data.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

##Inspecting the returned dada-class object of the 1st sample - it shows us the number of true ASVs are contained in the unique sequence reads:

dadaFs[[1]]

##Pooled samples option:
#dadaFsPool <- dada(derepFs, err = errF, pool=TRUE, multithread = TRUE)

#dadaRsPool <- dada(derepRs, err = errR, pool=TRUE, multithread = TRUE)

##Inspecting the returned dada-class object:

#dadaFsPool[[1]]

##Merge paired reads
##We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#mergers <- mergePairs(dadaFsPool, derepFs, dadaRsPool, derepRs, verbose=TRUE)

##Inspect the merger data.frame from the first sample

head(mergers[[1]])

##The mergers object is a list of data.frames from each sample. Each data.frame contains the merged $sequence, its $abundance, and the indices of the $forward and $reverse sequence variants that were merged. Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.

#Sequence table construction

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

##The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants.

##Chimeras removal
##Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

##Distribution of non-chimeric sequence lengths (Histogram)
hist(nchar(getSequences(seqtab.nochim)), main="Distribution of sequence lengths")

##Track reads through pipeline steps

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##Save table
write.table(track, file="dada2_output_track_reads.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

##Taxonomy Assignment

##Using the naive Bayesian classifier method

taxa <- assignTaxonomy(seqtab.nochim, "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/rdp_tax_files/rdp_train_set_16.fa.gz", multithread=TRUE)

#For >
set.seed(2000) #ASVs = 10,640
taxa_80 <- assignTaxonomy(seqtab.nochim, "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/rdp_tax_files/rdp_train_set_16.fa.gz", multithread=TRUE, minBoot = 80)

###Add species level
taxa <- addSpecies(taxa, "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/rdp_tax_files/rdp_species_assignment_16.fa.gz")

##Let’s inspect the taxonomic assignments:

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

##Save table
write.table(taxa, file="taxonomy_dada2_output_naive_bayes_maxee_2_4.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

##Construct phylogenetic tree

##Phylogenetic relatedness is commonly used to inform downstream analyses, especially the calculation of phylogeny-aware distances between microbial communities. The DADA2 sequence inference method is reference-free, so we must construct the phylogenetic tree relating the inferred sequence variants de novo. We begin by performing a multiple-alignment using the DECIPHER R package

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

##Export saved alignment file to cluster for tree construction.

##The phangorn R package is then used to construct a phylogenetic tree. Here we first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point.

#phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
#dm <- dist.ml(phang.align)
#treeNJ <- NJ(dm) # Note, tip order != sequence order
#fit = pml(treeNJ, data=phang.align)
#
### negative edges length changed to 0!
#
#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
#detach("package:phangorn", unload=TRUE)
write.tree(fitGTR$tree, file = "agt_swt.tre", append = FALSE, digits = 10, tree.names = FALSE))


##==============Phyloseq handoff==============
#Create phyloseq object

tree<-read.newick(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/rooted_agt_swt_maxee_2_4.tre")

##tree<-read.tree(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/rooted_agt_swt_maxee_2_4.tre",text = NULL, tree.names = NULL, skip = 0,comment.char = "", keep.multi = FALSE)

samdf<-read.csv("/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/metadata/agt_swt_metadata_study.csv")
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out
agt_swt <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa), phy_tree(tree))

#Save Phyloseq file
saveRDS(agt_swt, file="agt_swt_phyloseqfile.rds")


##Save phyloseq data that can be re-imported
save(agt_swt, file ="agt_swt_phyloseqData.RData")

##To load data in R:
load(“agt_swt_phyloseqData.RData”) 

##It is more convenient to use short names for our ASVs (e.g. ASV21) rather than the full DNA sequence when working with some of the tables and visualizations from phyloseq, but we want to keep the full DNA sequences for other purposes like merging with other datasets or indexing into reference databases like the Earth Microbiome Project. For that reason we’ll store the DNA sequences of our ASVs in the refseq slot of the phyloseq object, and then rename our taxa to a short string. That way, the short new taxa names will appear in tables and plots, and we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps)

dna <- Biostrings::DNAStringSet(taxa_names(agt_swt))
names(dna) <- taxa_names(agt_swt)
agt_swt <- merge_phyloseq(agt_swt, dna)
taxa_names(agt_swt) <- paste0("ASV", seq(ntaxa(agt_swt)))
agt_swt


##Save phyloseq data that can be re-imported
save(agt_swt, file ="agt_swt_phyloseqData.RData")

##To load data in R:
load(“agt_swt_phyloseqData.RData”) print(PS)

##Plot rarefaction data

#Rarefy data based on remaining non-chimeric reads from dada2 output tracking reads: agt_swt - 53,500 (removed 3 samples- BZEHG(12) and FQOSV(0) and ELHZE(25,277) and a total of 572 OTUs).

agt_swt.rarefied = rarefy_even_depth(agt_swt, rngseed=1, sample.size=53500, replace=F)

##Diversity Measurements

##Plot rarecurve
rarecurve(otu_table(agt_swt), step=10000, cex=0.5)
#High res save:
png("fig1_rarecurve_agt_swt.png", units="mm", width=200, height=150, res=300)
rarecurve(otu_table(agt_swt), step=10000, cex=0.8)+ theme(axis.text =element_text(size = 10), axis.title =element_text(size = 12))
dev.off()

##Test .ps
postscript("fig1_rarecurve_agt_swt.ps", width=10, height=7)

ggsave("fig1_rarecurve_agt_swt.png", plot = last_plot(), device = ps, path = NULL,
  scale = 1, width = 200, height = 150, units = c("mm"),
  dpi = 300, limitsize = TRUE) 




##Subset Data
agt.rarefied <-subset_samples(agt_swt.rarefied, site != "Soweto")
swt.rarefied <-subset_samples(agt_swt.rarefied, site != "Agincourt")
agt_lnob.rarefied <- subset_samples(agt_swt.rarefied, site != "Soweto")
agt_lnob.rarefied <- subset_samples(agt_lnob.rarefied, bmi_info != "Overweight")
swt_lnob.rarefied <- subset_samples(agt_swt.rarefied, site != "Agincourt")
swt_lnob.rarefied <- subset_samples(swt_lnob.rarefied, bmi_info != "Overweight")
agt_swt_lnob.rarefied  <- subset_samples(agt_swt.rarefied, bmi_info != "Overweight")
agt_swt_ln.rarefied <- subset_samples(agt_swt.rarefied, bmi_info != "Overweight")
agt_swt_ln.rarefied <- subset_samples(agt_swt_ln.rarefied, bmi_info != "Obese")
agt_swt_ob.rarefied <- subset_samples(agt_swt.rarefied, bmi_info != "Overweight")
agt_swt_ob.rarefied <- subset_samples(agt_swt_ob.rarefied, bmi_info != "Lean")

##Alpha Diversity (Phyloseq)

#rich = estimate_richness(ps.rarefied)
rich_agt = estimate_richness(agt_lnob.rarefied) #agt lean and obese samples only
rich_swt = estimate_richness(swt_lnob.rarefied) #swt lean and obese samples only
rich_agt_swt_lnob = estimate_richness(agt_swt_lnob.rarefied) #agt and swt lean and obese samples only
rich_agt_swt_site = estimate_richness(agt_swt.rarefied) #agt and swt lean, overweight and obese samples - all rarefied samples.
rich_agt_swt_ln = estimate_richness(agt_swt_ln.rarefied)
rich_agt_swt_ob = estimate_richness(agt_swt_ob.rarefied)

   
##Test whether the observed number of OTUs(ASVs) differ significantly between categories using Shannon diversity measure. We make a non-parametric test, the Wilcoxon rank-sum test (Mann-Whitney):

alpha_agt<-pairwise.wilcox.test(rich_agt$Shannon, sample_data(agt_lnob.rarefied)$bmi_info)
alpha_swt<-pairwise.wilcox.test(rich_swt$Shannon, sample_data(swt_lnob.rarefied)$bmi_info)
alpha_agt_swt_lnob<-pairwise.wilcox.test(rich_agt_swt_lnob$Shannon, sample_data(agt_swt_lnob.rarefied)$bmi_info)
alpha_agt_swt_site<-pairwise.wilcox.test(rich_agt_swt_site$Shannon, sample_data(agt_swt.rarefied)$site)
alpha_agt_swt_ln<-pairwise.wilcox.test(rich_agt_swt_ln$Shannon, sample_data(agt_swt_ln.rarefied)$site)
alpha_agt_swt_ob<-pairwise.wilcox.test(rich_agt_swt_ob$Shannon, sample_data(agt_swt_ob.rarefied)$site)



##Plot Richness
o)
##alpha_agt_swt_lnob<-pairwise.wilcox.test(rich_agt_swt_lnob$Shannon, sample_data(agt_swt_lnob.rarefied)$bmi_info)
alpha_agt_swt_site<-pairwise.wilcox.t
plot_rich_agt_all<-plot_richness(agt.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + geom_boxplot() #Plots measure for lean, overweight and obese categories
plot_rich_swt_all<-plot_richness(swt.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + geom_boxplot()
plot_rich_agt_swt_all<-plot_richness(agt_swt.rarefied, x="site", measures= "Shannon", color ="bmi_info") + geom_boxplot() 

##Save combined plot
require(cowplot)
plot_grid(plot_rich_agt_all, plot_rich_swt_all, plot_rich_agt_swt_all, labels = c('A', 'B', 'C'), align = 'h')

dir.create(path = "./alpha_diversity/")
ggsave("shannon_agt_swt_site_ovwt_all.pdf", plot=last_plot(), path ="./alpha_diversity/")

#High Resolution Save example
ggsave("richness_agt_all.png", units="in", width=4, height=4, dpi=300)


plot_rich_agt_lnob<-plot_richness(agt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + geom_boxplot() #Plots measure for only lean and obese categories
plot_rich_swt_lnob<-plot_richness(swt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + geom_boxplot()
plot_rich_agt_swt_lnob<-plot_richness(agt_swt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + geom_boxplot()
plot_rich_agt_swt_all_site<-plot_richness(agt_swt.rarefied, x="site", measures= "Shannon", color ="site") + geom_boxplot() #Plots measure for site comparison - all data: Agt vs Swt

##plot_rich_swt_lnob<-plot_richness(swt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + theme(legend.position = "none") + geom_boxplot() #To exclude legends from plot

##plot_richness(swt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + theme(legend.position = "none") + geom_boxplot() + xlab("BMI group") #Modify x-axis label

##plot_richness(swt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust=0.5), plot.title = element_text(face = 'plain', hjust =0)) + geom_boxplot() + labs(title = "Soweto", hjust =0, x ="BMI group") ##plain title, aligned to the left.

plot_rich_agt_lnob<-plot_richness(agt_lnob.rarefied, x="bmi_info", measures= "Shannon") + geom_boxplot() #Plots measure for only lean and obese categories
plot_rich_swt_lnob<-plot_richness(swt_lnob.rarefied, x="bmi_info", measures= "Shannon")  + geom_boxplot()
plot_rich_agt_swt_lnob<-plot_richness(agt_swt_lnob.rarefied, x="bmi_info", measures= "Shannon") + geom_boxplot()
plot_rich_agt_swt_all_site<-plot_richness(agt_swt.rarefied, x="site", measures= "Shannon") + geom_boxplot() #Plots measure for site comparison - all data: Agt vs Swt
plot_rich_as_ob<-plot_richness(agt_swt_ob.rarefied, x="site", measures= "Shannon", color ="site") + geom_boxplot()
plot_rich_as_ln<-plot_richness(agt_swt_ln.rarefied, x="site", measures= "Shannon", color ="site") + geom_boxplot()

##Modified Plots
grob <- grobTree(textGrob("p = 0.77", x=0.45,  y=0.92, hjust=0, gp=gpar(fontsize=12)))
plot_rich_agt_lnob<-plot_richness(agt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust=0.5, size =13, margin = margin(t = 0, r = 0, b = -17, l = 0)), axis.text.y = element_text(size = 13), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + geom_boxplot() + labs(title = "Bushbuckridge", hjust =0, x ="BMI group") + annotation_custom(grob)  

grob <- grobTree(textGrob("p = 0.41", x=0.45,  y=0.92, hjust=0, gp=gpar(fontsize=12)))
plot_rich_swt_lnob<-plot_richness(swt_lnob.rarefied, x="bmi_info", measures= "Shannon", color ="bmi_info") + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust=0.5, size =13, margin = margin(t = 0, r = 0, b = -17, l = 0)), axis.text.y = element_text(size = 13), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + geom_boxplot() + labs(title = "Soweto", hjust =0, x ="BMI group") + annotation_custom(grob) 

grob <- grobTree(textGrob("p = 0.01", x=0.45,  y=0.92, hjust=0, gp=gpar(fontsize=12)))
plot_rich_agt_swt_all_site <- plot_richness(agt_swt.rarefied, x="site", measures= "Shannon", color = "site") + geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, size =13, hjust=0.35, margin = margin(t = 0, r = 0, b = -20, l = 0)), axis.text.y = element_text(size = 13), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + geom_boxplot() + labs(title = "Site Comparisons", hjust=0, x="Study site") + scale_x_discrete(labels = c("Bushbuckridge","Soweto")) + annotation_custom(grob)


agt_swt_lnob_combo<-plot_grid(plot_rich_as_ln, plot_rich_as_ob, labels = c('A', 'B'), align = 'h')
ggsave("richness_plot_agt_swt_lnob_combo.pdf", plot=agt_swt_lnob_combo, path ="./alpha_diversity/")

##Save files:
##One method to save as high res grid:
##Method used: #Fig2A
png("fig2_combined_richness3.png", units="mm", width=250, height=300, res=300)
plot_grid(plot_rich_agt_lnob, plot_rich_swt_lnob, plot_rich_agt_swt_all_site, labels = c('A', 'B', 'C'), align = 'h')
dev.off()

##pdf
pdf("fig2_combined_richness3.pdf", width=10, height=12)
plot_grid(plot_rich_agt_lnob, plot_rich_swt_lnob, plot_rich_agt_swt_all_site, labels = c('A', 'B', 'C'), align = 'h')
dev.off()


##Save plots as grid
require(cowplot)
theme_set(theme_cowplot(font_size=12)) # reduce default font size
combined_richness=plot_grid(plot_rich_agt_lnob, plot_rich_swt_lnob, plot_rich_agt_swt_all_site, labels = c('A', 'B', 'C'), align = 'h')

ggsave("richness_plot_agt_swt_site.pdf", plot=combined, path ="./alpha_diversity/")
#High Resolution Save 
ggsave("combined_richness.png", units="in", width=4, height=4, dpi=300, plot=combined_richness)
##Beta Diversity

##PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(agt_lnob.rarefied, method="unifrac", weighted=F)
ordination = ordinate(agt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_lnob.rarefied, ordination, color="site", shape = "bmi_info) + theme(aspect.ratio=1) + labs(color="BMI group")

#done for agt_lnob.rarefied (agt_lnob_rare.plot_wt, agt_lnob_rare.plot_unwt), swt_lnob.rarefied (swt_lnob_rare.plot_unwt, swt_lnob_rare.plot_wt), agt_swt.rarefied(agt_swt_rare.plot_wt, agt_swt_rare.plot_unwt), agt_swt_lnob.rarefied (agt_swt_lnob_rare.plot_unwt, agt_swt_lnob_rare.plot_wt)

plot_ordination(agt_swt_lnob.rarefied, ordination, color="site", shape = "bmi_info") + geom_point(size=0.05)

dir.create(path ="./beta_diversity")
ggsave("./beta_diversity/swt_lnob_unweightedUniFrac.pdf", plot=last_plot())

##Save plots as grid
combined_beta_unweighted=plot_grid(agt_lnob_rare.plot_unwt, swt_lnob_rare.plot_unwt, agt_swt_lnob_rare.plot_unwt, labels = c('A', 'B', 'C'), align = 'h')




##Beta Diversity (https://micca.readthedocs.io/en/latest/phyloseq.html)
##Permanova
##Test whether the categories differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:

dir.create("./beta_diversity")

##Beta Diversity - weighted/unweighted Unifrac (Phyloseq package)

wunifrac_dist = phyloseq::distance(agt_swt_lnob.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
grob <- grobTree(textGrob("p = 0.41", x=0.70,  y=0.92, hjust=0, gp=gpar(fontsize=12)))
weighted_agt_swt_lnob<-plot_ordination(agt_swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1) + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Weighted UniFrac") + annotation_custom(grob)
permanova_weighted_agt_swt_lnob<-adonis(wunifrac_dist ~ sample_data(agt_swt_lnob.rarefied)$bmi_info)

####LabelTest
weighted_agt_swt_lnob_label<-plot_ordination(agt_swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1) + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = expression(~bold(B)~ "  Weighted UniFrac")) + annotation_custom(grob)

##High res save:
#
#png("weighted_agt_swt_lnob.png", units="mm", width=250, height=100, res=300)
#
#plot_ordination(agt_swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1)
#ggsave("./beta_diversity/wt_uniFrac_agt_swt_lnob.pdf", plot=last_plot())
#permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt_lnob.rarefied)$bmi_info) #Test whether the bmi categories differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis with weighted or wnweighted  uniFrac distance

wunifrac_dist = phyloseq::distance(agt_swt_lnob.rarefied, method="unifrac", weighted=F)
ordination = ordinate(agt_swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
grob <- grobTree(textGrob("p = 0.01", x=0.70,  y=0.92, hjust=0, gp=gpar(fontsize=12)))
unweighted_agt_swt_lnob<-plot_ordination(agt_swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Unweighted UniFrac") + annotation_custom(grob)
permanova_unweighted_agt_swt_lnob<-adonis(wunifrac_dist ~ sample_data(agt_swt_lnob.rarefied)$bmi_info)

####LabelTest
unweighted_agt_swt_lnob_label<-plot_ordination(agt_swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1) + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = expression(~bold(A)~ "  Unweighted UniFrac")) + annotation_custom(grob)



#High res save:

png("unweighted_agt_swt_lnob.png", units="mm", width=250, height=100, res=300)

#Or save grid:
png("SF1ed_all_site_lnob_comparison.png", units="mm", width=250, height=100, res=300)
plot_grid(unweighted_agt_swt_lnob, weighted_agt_swt_lnob, labels = c('A', 'B'), align = 'h')
dev.off()

#pdf save label
pdf("SF1ed_label_all_site_lnob_comparison.pdf", width=10, height=8)
plot_grid(unweighted_agt_swt_lnob_label, weighted_agt_swt_lnob_label, align = 'h')
dev.off()

wunifrac_dist = phyloseq::distance(agt_swt_lnob.rarefied, method="unifrac", weighted=F) #Unweighted
ordination = ordinate(agt_swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1)
ggsave("./beta_diversity/unwt_uniFrac_agt_swt_lnob.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt_lnob.rarefied)$bmi_info)


wunifrac_dist = phyloseq::distance(agt_swt_lnob.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_lnob.rarefied, ordination, color="site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/wt_uniFrac_agt_swt_lnob_site.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt_lnob.rarefied)$site) #Test whether the bmi categories differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis with weighted or wnweighted  uniFrac distance


wunifrac_dist = phyloseq::distance(agt_swt_lnob.rarefied, method="unifrac", weighted=T) #Unweighted
ordination = ordinate(agt_swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_lnob.rarefied, ordination, color="site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/unwt_uniFrac_agt_swt_lnob_site.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt_lnob.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_swt.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_swt.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt.rarefied, ordination, color= "site", shape = "bmi_info") + theme(aspect.ratio=1)
ggsave("./beta_diversity/wt_uniFrac_agt_swt_all_site_shapes.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_swt.rarefied, method="unifrac", weighted=F) #Unweighted
ordination = ordinate(agt_swt.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt.rarefied, ordination, color= "site", shape = "bmi_info") + theme(aspect.ratio=1)
ggsave("./beta_diversity/unwt_uniFrac_agt_swt_all_site_shapes.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt.rarefied)$site)


wunifrac_dist = phyloseq::distance(agt_swt.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_swt.rarefied, method="PCoA", distance=wunifrac_dist)
grob <- grobTree(textGrob("p < 0.001", x=0.70,  y=0.96, hjust=0, gp=gpar(fontsize=12)))
weighted_agt_swt_all<-plot_ordination(agt_swt.rarefied, ordination, color= "site") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "Site", title = "Weighted UniFrac") + scale_color_discrete(name="Site", breaks=c("Agincourt", "Soweto"), labels=c("Bushbuckridge", "Soweto")) + annotation_custom(grob)
permanova_weighted_agt_swt_all_site<-adonis(wunifrac_dist ~ sample_data(agt_swt.rarefied)$site)

####LabelTest
weighted_agt_swt_all_label<-plot_ordination(agt_swt.rarefied, ordination, color= "site") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "Site", title = expression(~bold(B)~ "  Weighted UniFrac")) + scale_color_discrete(name="Site", breaks=c("Agincourt", "Soweto"), labels=c("Bushbuckridge", "Soweto")) + annotation_custom(grob)



wunifrac_dist = phyloseq::distance(agt_swt.rarefied, method="unifrac", weighted=F)
ordination = ordinate(agt_swt.rarefied, method="PCoA", distance=wunifrac_dist)
grob <- grobTree(textGrob("p < 0.001", x=0.70,  y=0.96, hjust=0, gp=gpar(fontsize=12)))
unweighted_agt_swt_all<-plot_ordination(agt_swt.rarefied, ordination, color= "site") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "Site",fill = c("Bushbuckridge", "Soweto"), title = "Unweighted UniFrac") + scale_color_discrete(name="Site", breaks=c("Agincourt", "Soweto"), labels=c("Bushbuckridge", "Soweto")) + annotation_custom(grob)
permanova_unweighted_agt_swt_all_site<-adonis(wunifrac_dist ~ sample_data(agt_swt.rarefied)$site)


####LabelTest
unweighted_agt_swt_all<-plot_ordination(agt_swt.rarefied, ordination, color= "site") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "Site",fill = c("Bushbuckridge", "Soweto"), title = expression(~bold(A)~ "   Unweighted UniFrac")) + scale_color_discrete(name="Site", breaks=c("Agincourt", "Soweto"), labels=c("Bushbuckridge", "Soweto")) + annotation_custom(grob)



wunifrac_dist = phyloseq::distance(agt_lnob.rarefied, method="unifrac", weighted=F) 
ordination = ordinate(agt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
grob <- grobTree(textGrob("p = <0.001", x=0.63,  y=0.88, hjust=0, gp=gpar(fontsize=12)))
unweighted_agt_lnob<-plot_ordination(agt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Bushbuckridge") + annotation_custom(grob)
permanova_unweighted_agt_lnob<-adonis(wunifrac_dist ~ sample_data(agt_lnob.rarefied)$bmi_info)

####LabelTest
unweighted_agt_lnob_label<-plot_ordination(agt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = expression(~bold(C)~ "  Bushbuckridge")) + annotation_custom(grob)



wunifrac_dist = phyloseq::distance(agt_lnob.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
#grob <- grobTree(textGrob("p = <0.001", x=0.63,  y=0.88, hjust=0, gp=gpar(fontsize=12)))
#unweighted_agt_lnob<-plot_ordination(agt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Bushbuckridge") + annotation_custom(grob)
permanova_weighted_agt_lnob<-adonis(wunifrac_dist ~ sample_data(agt_lnob.rarefied)$bmi_info)


wunifrac_dist = phyloseq::distance(swt_lnob.rarefied, method="unifrac", weighted=F)
ordination = ordinate(swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
grob <- grobTree(textGrob("p = 0.81", x=0.63,  y=0.88, hjust=0, gp=gpar(fontsize=12)))
unweighted_swt_lnob<-plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Soweto") + annotation_custom(grob)
permanova_unweighted_swt_lnob<-adonis(wunifrac_dist ~ sample_data(swt_lnob.rarefied)$bmi_info)

####LabelTest
unweighted_swt_lnob_label<-plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = expression(~bold(D)~ "  Soweto")) + annotation_custom(grob)



wunifrac_dist = phyloseq::distance(swt_lnob.rarefied, method="unifrac", weighted=T) 
ordination = ordinate(swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
#grob <- grobTree(textGrob("p = 0.81", x=0.63,  y=0.88, hjust=0, gp=gpar(fontsize=12)))
#unweighted_swt_lnob<-plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Soweto") + annotation_custom(grob)
permanova_weighted_swt_lnob<-adonis(wunifrac_dist ~ sample_data(swt_lnob.rarefied)$bmi_info)


png("fig3_beta_diversity_all_site_unweighted_lnob_comparison.png", units="mm", width=250, height=200, res=300)
plot_grid(unweighted_agt_swt_all, weighted_agt_swt_all, unweighted_agt_lnob, unweighted_swt_lnob, labels = c('A', 'B', 'C', 'D'), align = 'h')
dev.off()

#pdf save label
pdf("fig3_beta_diversity_all_site_unweighted_lnob_comparison_label.pdf", width=10, height=12)
plot_grid(unweighted_agt_swt_all_label, weighted_agt_swt_all_label, unweighted_agt_lnob_label, unweighted_swt_lnob_label, align = 'h')
dev.off()



ggsave("./beta_diversity/wt_uniFrac_agt_swt_all_site.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_swt.rarefied, method="unifrac", weighted=F) #Unweighted
ordination = ordinate(agt_swt.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt.rarefied, ordination, color= "site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/unwt_uniFrac_agt_swt_all_site.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_swt.rarefied)$site)



wunifrac_dist = phyloseq::distance(agt_swt_ln.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_swt_ln.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_ln.rarefied, ordination, color= "site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/wt_uniFrac_agt_swt_ln_site.pdf", plot=last_plot())
permanova_weighted_agt_swt_ln<-adonis(wunifrac_dist ~ sample_data(agt_swt_ln.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_swt_ln.rarefied, method="unifrac", weighted=F)
ordination = ordinate(agt_swt_ln.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_ln.rarefied, ordination, color= "site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/unwt_uniFrac_agt_swt_ln_site.pdf", plot=last_plot())
permanova_unweighted_agt_swt_ln<-adonis(wunifrac_dist ~ sample_data(agt_swt_ln.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_swt_ob.rarefied, method="unifrac", weighted=T)
ordination = ordinate(agt_swt_ob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_ob.rarefied, ordination, color= "site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/wt_uniFrac_agt_swt_ob_site.pdf", plot=last_plot())
permanova_weighted_agt_swt_ob<-adonis(wunifrac_dist ~ sample_data(agt_swt_ob.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_swt_ob.rarefied, method="unifrac", weighted=F)
ordination = ordinate(agt_swt_ob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_swt_ob.rarefied, ordination, color= "site") + theme(aspect.ratio=1)
ggsave("./beta_diversity/unwt_uniFrac_agt_swt_ob_site.pdf", plot=last_plot())
permanova_unweighted_agt_swt_ob<-adonis(wunifrac_dist ~ sample_data(agt_swt_ob.rarefied)$site)

wunifrac_dist = phyloseq::distance(agt_lnob.rarefied, method="unifrac", weighted=T) #Fig2B
ordination = ordinate(agt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
weighted_agt_lnob<-plot_ordination(agt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Bushbuckridge") 

wunifrac_dist = phyloseq::distance(swt_lnob.rarefied, method="unifrac", weighted=T) #Fig2B
ordination = ordinate(swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
weighted_swt_lnob<-plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Soweto")

##Method used: #Fig2B
png("weighted.png", units="mm", width=250, height=100, res=300)
plot_grid(weighted_agt_lnob, weighted_swt_lnob, labels = c('A', 'B'), align = 'h')
dev.off()

ggsave("./beta_diversity/wt_uniFrac_agt_lnob.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_lnob.rarefied)$bmi_info)

wunifrac_dist = phyloseq::distance(agt_lnob.rarefied, method="unifrac", weighted=F) #Unweighted - Fig2B
ordination = ordinate(agt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(agt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1) + labs(color = "BMI group")
unweighted_agt_lnob<-plot_ordination(agt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Bushbuckridge")

#ggsave("./beta_diversity/unwt_uniFrac_agt_lnob.pdf", plot=last_plot())
#permanova_weighted<-adonis(wunifrac_dist ~ sample_data(agt_lnob.rarefied)$bmi_info)

wunifrac_dist = phyloseq::distance(swt_lnob.rarefied, method="unifrac", weighted=F) #Unweighted - Fig2B
ordination = ordinate(swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
unweighted_swt_lnob<-plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Soweto")

##Method used: #Fig2B
png("unweighted.png", units="mm", width=250, height=100, res=300)
plot_grid(unweighted_agt_lnob, unweighted_swt_lnob, labels = c('A', 'B'), align = 'h')
dev.off()

#wunifrac_dist = phyloseq::distance(swt_lnob.rarefied, method="unifrac", weighted=T) #Fig2B
#ordination = ordinate(swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
#plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1) + labs(color = "BMI group")
#ggsave("./beta_diversity/wt_uniFrac_swt_lnob.pdf", plot=last_plot())
#permanova_weighted<-adonis(wunifrac_dist ~ sample_data(swt_lnob.rarefied)$bmi_info)

#wunifrac_dist = phyloseq::distance(swt_lnob.rarefied, method="unifrac", weighted=F) #Unweighted - Fig2B
#ordination = ordinate(swt_lnob.rarefied, method="PCoA", distance=wunifrac_dist)
#unweighted_swt_lnob<-plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1, axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), plot.title = element_text(face = 'plain', hjust =0, size = 15)) + labs(color = "BMI group", title = "Soweto")
#
#
#plot_ordination(swt_lnob.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1) + labs(color = "BMI group")
#ggsave("./beta_diversity/unwt_uniFrac_swt_lnob.pdf", plot=last_plot())
permanova_weighted<-adonis(wunifrac_dist ~ sample_data(swt_lnob.rarefied)$bmi_info)

#ordination = ordinate(ps.rarefied, method="PCoA", distance=wunifrac_dist)
#plot_ordination(ps.rarefied, ordination, color="bmi_info") + theme(aspect.ratio=1)

##Plot grid for fig3
png("fig3_beta_diversity.png", units="mm", width=250, height=200, res=300)
plot_grid(unweighted_agt_swt_all, weighted_agt_swt_all,  labels = c('A', 'B'), align = 'h')
dev.off()

=======================================================================
##DESeq2 Differential Abundance Testing

##First subset samples, then remove (filter out) taxa not seen more than 2 times in at least 5% of the samples. This protects against an OTU with small mean & trivially large C.V. This step was primarily used for deseq2 analyses.

#Subset samples

agt_swt_sub<-subset_samples(agt_swt, SampleID != c("BZEHG","FQOSV", "ELHZE"))

filt_agt_swt_sub= filter_taxa(agt_swt_sub, function(x) sum(x > 2) > (0.05*length(x)), TRUE) ##This resulted in 1,395 high abundance ASVs (taxa)

##!Remember to subset other comparisons from the filt_agt_swt_sub data!

##Subset Data
filt_agt_swt_lnob <-subset_samples(filt_agt_swt_sub, bmi_info != "Overweight") 
filt_agt_lnob <- subset_samples(filt_agt_swt_lnob, site != "Soweto")
filt_swt_lnob <- subset_samples(filt_agt_swt_lnob, site != "Agincourt")
filt_agt_swt_ln <- subset_samples(filt_agt_swt_lnob, bmi_info != "Obese")
filt_agt_swt_ob <- subset_samples(filt_agt_swt_lnob, bmi_info != "Lean")




##Phyloseq to DESEq2 object - Condition of interest should always be the last term - covariates should come before.

diagdds = phyloseq_to_deseq2(filt_agt_lnob, ~ batch + bmi_info) #Agt lean vs obese comparisons whilst controlling for potential batch effects
diagdds = phyloseq_to_deseq2(filt_swt_lnob, ~ bmi_info) #Swt lean vs obese comparisons - sequenced in one batch
diagdds = phyloseq_to_deseq2(filt_agt_swt_lnob, ~ batch + bmi_info) #Lean vs obese comparisons of pilot cohort excluding overweight samples whilst controlling for batch effects
#diagdds = phyloseq_to_deseq2(filt_agt_swt_lnob, ~ batch + site + bmi_info) #Lean vs obese comparisons of pilot cohort excluding overweight samples whilst controlling for site and batch effects
diagdds = phyloseq_to_deseq2(filt_agt_swt_sub, ~ batch + bmi_info + site) #Agt vs Soweto comparisons of all data including overweight samples
diagdds = phyloseq_to_deseq2(filt_agt_swt_sub_ob, ~ site) #Comparing compositional differences between sites in only obese samples
diagdds = phyloseq_to_deseq2(filt_agt_swt_ln, ~ batch + site) #Comparing compositional differences between sites in only lean samples whilst controlling for potential batch effects


###Controlling for htn
#diagdds = phyloseq_to_deseq2(filt_agt_swt_lnob, ~ batch + htn + site + bmi_info) #Lean vs obese comparisons of pilot cohort excluding overweight samples whilst controlling for batch effects
#diagdds = phyloseq_to_deseq2(filt_agt_swt_sub, ~ batch + bmi_info + htn + site) #Agt vs Soweto comparisons of all data including overweight samples
#
###Htn as the outcome variable whilst controlling for BMI as a continuous variable using log10 BMI values
#diagdds = phyloseq_to_deseq2(filt_agt_lnob, ~ batch + log_bmi + htn)
#diagdds = phyloseq_to_deseq2(filt_swt_lnob, ~ log_bmi + htn)
#diagdds = phyloseq_to_deseq2(filt_agt_swt_sub, ~ batch + log_bmi + htn)  #control for batch and log bmi.
#diagdds = phyloseq_to_deseq2(filt_agt_swt_sub, ~ batch + log_bmi + site + htn) #control for batch, log bmi and site.
#
####Controlling for both site and bmi_info:
###diagdds = phyloseq_to_deseq2(filt_agt_swt_sub, ~ batch + bmi_info + site)
###
####diagdds = phyloseq_to_deseq2(swted, ~ bmi_info + htn)# account (control) for differences due to bmi_info whilst estimating the effect due to htn. Condition of interest should always be the last term - covariates should come before.
###
####For agt_swt site analyses whilst controlling for batch effects:res = results(diagdds)
###diagdds = phyloseq_to_deseq2(filt_agt_swt, ~ batch + site)
###
####Controlling for both site and bmi_info:
###diagdds = phyloseq_to_deseq2(filt_agt_swt_05, ~ batch + bmi_info_2 + site)
####For Swt and Agt Lean vs Obese Comparisons:
###diagdds = phyloseq_to_deseq2(agt_obln, ~ bmi_info) #No batch control because one batch for all samples in Swt; and for Agt, batch categorizatio is pretty much equivalent to lean/obese categories.
###
####For Agt_Swt site comparison of lean and obese
###diagdds = phyloseq_to_deseq2(agt_swt_ob, ~ site) #all the obese samples were done in one batch; controlling for batch in the lean samples is equivalent to site compparisons, so no batch comparisons were done.
##


# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
##Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.


res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_agt_swt_lnob*)[rownames(sigtab),], "matrix"))
head(sigtab)

##First, cleaning up the table a little for legibility.
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

##Writing Deseq Report
dir.create(path="./deseq2")
write.csv(as.data.frame(sigtab), file="./deseq2/filt_agt_batch_controlled.csv")
write.csv(as.data.frame(sigtab), file="./deseq2/filt_swt_lnob.csv")
write.csv(as.data.frame(sigtab), file="./deseq2/filt_agt_swt_lnob.csv")
write.csv(as.data.frame(sigtab), file="./deseq2/filt_agt_swt_sub.csv")
write.csv(as.data.frame(sigtab), file="./deseq2/filt_agt_swt_ob.csv")
write.csv(as.data.frame(sigtab), file="./deseq2/filt_agt_swt_ln.csv")


##Plot Results
sigtabgen<-sigtab
theme_set(theme_bw())

# Genus order <-this was used
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

####Text size adjustment in plot: ##For DESeq figures
agt_lnob_plot<- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1.3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=11), axis.text.y=element_text(size = 11, face = "italic"), legend.text=element_text(size=11, face ="italic"), axis.title=element_text(size=12), legend.title=element_text(size=12)) 

swt_lnob_plot<- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1.3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=11), axis.text.y=element_text(size = 11, face = "italic"), legend.text=element_text(size=11, face = "italic"), axis.title=element_text(size=12), legend.title=element_text(size=12))

agt_swt_lnob_plot<- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1.3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=11), axis.text.y=element_text(size = 11, face = "italic"), legend.text=element_text(size=11, face = "italic"), axis.title=element_text(size=12), legend.title=element_text(size=12))

agt_swt_ln_plot<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1.3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=11), axis.text.y=element_text(size = 11, face = "italic"), legend.text=element_text(size=11, face = "italic"), axis.title=element_text(size=12), legend.title=element_text(size=12))

agt_swt_ob_plot<- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1.3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=11), axis.text.y=element_text(size = 11, face = "italic"), legend.text=element_text(size=11, face = "italic"), axis.title=element_text(size=12), legend.title=element_text(size=12))

agt_swt_bmi_batch_ctrled_plot <- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1.3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=11), axis.text.y=element_text(size = 11, face = "italic"), legend.text=element_text(size=11, face = "italic"), axis.title=element_text(size=12), legend.title=element_text(size=12))



ggsave("filt_agt_swt_sub_site_batch_bmi_controlled.pdf", plot=last_plot())

png("ln_agt_vs_swt_left_ob_agt_vs_swt_right.png", units="mm", width=250, height=150, res=300) #Used height =150
plot_grid(agt_lnob_plot, swt_lnob_plot, labels = c('A', 'B'), align = 'h')
dev.off()

png("ln_vs_ob_agt_left_swt_right_1.png", units="mm", width=250, height=150, res=300) #Used height =150
plot_grid(agt_lnob_plot, swt_lnob_plot, labels = c('A', 'B'), align = 'h')
dev.off()

png("ln_agt_vs_swt_left_ob_agt_vs_swt_right_1.png", units="mm", width=250, height=150, res=300) #Used height =150
plot_grid(agt_swt_ln_plot, agt_swt_ob_plot, labels = c('A', 'B'), align = 'h')
dev.off()

theme_set(theme_bw())
png("ln_agt_vs_swt_left_A__ob_agt_vs_swt_right_B_agt_swt_bmi_batch_ctrled_all_C.png", units="mm", width=250, height=300, res=300) 
plot_grid(agt_swt_ln_plot, agt_swt_ob_plot,agt_swt_bmi_batch_ctrled_plot, labels = c('A', 'B', 'C'), align = 'h') #Figure5
dev.off()

png("lnob_agt_A_swt_B_agt_vs_swt_C.png", units="mm", width=250, height=300, res=300)             
plot_grid(agt_lnob_plot, swt_lnob_plot,agt_swt_lnob_plot, labels = c('A', 'B', 'C'), align = 'h') #Figure6
dev.off()


##Edited figures for phyla color uniformity across plots
#Set colours manually
#Actinobacteria - red
#Bacteroidetes - forestgreen
#Elusimicrobia - blue
#Firmicutes - orange
#Lentisphaerae - purple
#Proteobacteria - black
#Tenericutes - pink
#Verrucomicrobia - green
#Spirochaetes - brown
#NA - gray

#scale_colour_manual(values=c('red','forestgreen','blue','orange','purple','black','pink','green'))

#For agt_swt_ln
grob <- grobTree(textGrob("BBR      SWT", x=0.20,  y=0.93, hjust=0, gp=gpar(fontsize=11)))
agt_swt_ln_plot_color <- agt_swt_ln_plot + scale_colour_manual(values=c('red','forestgreen','blue','orange','purple','black','pink','green'), na.value= 'gray') + annotation_custom(grob)


#For agt_swt_ob
grob <- grobTree(textGrob("BBR       SWT", x=0.20,  y=0.93, hjust=0, gp=gpar(fontsize=11)))
agt_swt_ob_plot_color <- agt_swt_ob_plot + scale_colour_manual(values=c('forestgreen','orange','purple','black','pink','green'), na.value= 'gray') + annotation_custom(grob)

#For agt_swt_bmi_batch_ctrled
grob <- grobTree(textGrob("BBR    SWT", x=0.17,  y=0.94, hjust=0, gp=gpar(fontsize=11)))
agt_swt_bmi_batch_ctrled_plot_color <- agt_swt_bmi_batch_ctrled_plot + scale_colour_manual(values=c('forestgreen','orange','black','pink'), na.value= 'gray') + annotation_custom(grob)

#Figure5
theme_set(theme_bw())
png("fig5_ln_agt_vs_swt_left_A_ob_agt_vs_swt_right_B_agt_swt_bmi_batch_ctrled_all_C.png", units="mm", width=250, height=300, res=300)
plot_grid(agt_swt_ln_plot_color, agt_swt_ob_plot_color,agt_swt_bmi_batch_ctrled_plot_color, labels = c('A', 'B', 'C'), align = 'h') #Figure5
dev.off()

#pdf save:
pdf("fig5_ln_agt_vs_swt_left_A_ob_agt_vs_swt_right_B_agt_swt_bmi_batch_ctrled_all_C.pdf", width=10, height=12)
plot_grid(agt_swt_ln_plot_color, agt_swt_ob_plot_color,agt_swt_bmi_batch_ctrled_plot_color, labels = c('A', 'B', 'C'), align = 'h') #Figure5
dev.off()


#For agt_lnob
grob <- grobTree(textGrob("LN    OB", x=0.31,  y=0.97, hjust=0, gp=gpar(fontsize=11)))
agt_lnob_plot_color <- agt_lnob_plot + scale_colour_manual(values=c('forestgreen','orange','black','green'),na.value= 'gray') + annotation_custom(grob)

#For swt_lnob
grob <- grobTree(textGrob("LN    OB", x=0.27,  y=0.97, hjust=0, gp=gpar(fontsize=11)))
swt_lnob_plot_color <- swt_lnob_plot + scale_colour_manual(values=c('forestgreen','orange','black','green')) + annotation_custom(grob) 

#For agt_swt_lnob
grob <- grobTree(textGrob("LN    OB", x=0.33,  y=0.95, hjust=0, gp=gpar(fontsize=11)))
agt_swt_lnob_plot_color <- agt_swt_lnob_plot + scale_colour_manual(values=c('forestgreen','blue','orange','purple','black','brown','green'), na.value= 'gray') + annotation_custom(grob)

#Figure 6
png("fig6_lnob_agt_A_swt_B_agt_vs_swt_C.png", units="mm", width=250, height=300, res=300)             
plot_grid(agt_lnob_plot_color, swt_lnob_plot_color,agt_swt_lnob_plot_color, labels = c('A', 'B', 'C'), align = 'h') 
dev.off()

#pdf save:
pdf("fig6_lnob_agt_A_swt_B_agt_vs_swt_C.pdf", width=10, height=12)
plot_grid(agt_lnob_plot_color, swt_lnob_plot_color,agt_swt_lnob_plot_color, labels = c('A', 'B', 'C'), align = 'h')
dev.off()








===================================================

Enterotyping
dir.create(./enterotype)
setwd("/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/july_analyses/enterotype")


##Agt_Swt
microbio.taxonomy= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/taxonomy_dada2_output_naive_bayes_maxee_2_4_ed.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)
microbio.tree = read.newick(file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/rooted_agt_swt_maxee_2_4.tre")
microbio.meta= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/metadata/agt_swt_metadata_study.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)
microbio.otus= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/otu_table_maxee_2_4_ed.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)
microbio.taxonomy= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/dada_cluster_Dec2019/taxa_80_taxonomy_dada2_output_naive_bayes_maxee_2_4_ed.txt", header =T, sep ="\t", row.names=1, check.names=FALSE) #Using taxa_80
##Agt
microbio.meta= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/metadata/agt_metadata_study.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)
microbio.otus= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/agt_otu_table_maxee_2_4_ed.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)

##Swt
microbio.meta= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/metadata/swt_metadata_study.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)
microbio.otus= read.table(file="/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/dada2/swt_otu_table_maxee_2_4_ed.txt", header =T, sep ="\t", row.names=1, check.names=FALSE)


# OTU rarefaction
# By default, it uses the minimum number. Verify with rowSums(microbio.rare)
#microbio.rare = Rarefy(microbio.otus)$otu.tab.rff

#However, I rarefied to a depth of 53,500:
microbio.rare = Rarefy(microbio.otus, depth = 53500)$otu.tab.rff

# Calculate OTU relative frequencies
microbio.relative = t(microbio.otus/rowSums(microbio.otus))

# To create phylotypes for each taxonomic level ----
# Sums all the OTUs with the same taxonomy

# The following code was taken from a tutorial found here:
# http://www.r-bloggers.com/from-otu-table-to-heatmap/
# https://learningomics.wordpress.com/2013/02/23/from-otu-table-to-heatma/
#
# Functions
# Function to separate taxonomies
extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}


# Function to sumarize the data at different taxonomic levels
otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa)){
    taxa = colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length as the number of columns in the OTU table")
    return;
  }
level.names = sapply(as.character(taxa),
                       function(x)
                         extract.name.level(x,level=level))
  t(apply(x, 1,
          function(y)
            tapply(y,level.names,sum)))
}

# Separate the abundance data from the taxonomy
taxa.names = microbio.taxonomy$Taxonomy #Edit RDP taxonomy file to concatenate levels with ";".
dat2 = t(microbio.otus)

# Remove samples with low sequence count
s_abundances = apply(dat2,2,sum)

# Separate the data that are above and below the threshold (1000 reads in this case)
bads = dat2[,s_abundances<53500] #53,500 threshold for Agt_Swt 
goods = dat2[,s_abundances>53500]

# Number of samples that are above and below the threshold
ncol(goods)
ncol(bads)

# Keep only 'good' samples
dat2 = goods

dat2 = scale(dat2, center=F, scale=colSums(dat2))
dat2 <-t(dat2)

#dat3<-t(dat2)
#rownames(dat3)<-taxa.names
#dat3melt<-melt(dat3)
#dat4<-cbind(as.data.frame(matrix(unlist(strsplit(as.character(dat3melt[,1]),split=';')), nrow=length(dat3melt[,1]),byrow=T)), dat3melt)
#aggregate(value~V6, data=dat4, FUN=mean)

# Separate objects for each taxonomic level
# RDP taxonomy has 7 levels, starting from 1 (Kingdom), 2 (Phylum), ...
d.phylum = otu2taxonomy(dat2,level=2,taxa=taxa.names)
d.class = otu2taxonomy(dat2,level=3,taxa=taxa.names)
d.order = otu2taxonomy(dat2,level=4,taxa=taxa.names)
d.family = otu2taxonomy(dat2,level=5,taxa=taxa.names)
d.genus = otu2taxonomy(dat2,level=6,taxa=taxa.names)
d.species = otu2taxonomy(dat2,level=7,taxa=taxa.names)

# Transpose the tables and export the files
 phylum2 <-t(d.phylum)
 class2 <-t(d.class)
 order2 <-t(d.order)
 family2 <-t(d.family)
 genus2 <-t(d.genus)
 species2 <-t(d.species)

 # Create a new folder and save taxonomy (phylotype) tables
 dir.create(path = "./phylotypes/")

write.table(phylum2, file="phylotypes/phyla.txt", col.names=NA,row.names=TRUE,
 sep="\t", quote=FALSE)
 write.table(class2, file="phylotypes/classes.txt", col.names=NA,row.names=TRUE,
 sep="\t", quote=FALSE)
 write.table(order2, file="phylotypes/orders.txt", col.names=NA,row.names=TRUE,
 sep="\t", quote=FALSE)
 write.table(family2, file="phylotypes/families.txt", col.names=NA,row.names=TRUE,
 sep="\t", quote=FALSE)
 write.table(genus2, file="phylotypes/genera.txt", col.names=NA,row.names=TRUE,
 sep="\t", quote=FALSE)
 write.table(species2, file="phylotypes/species.txt", col.names=NA,row.names=TRUE,
 sep="\t", quote=FALSE)


# Descriptive statistics by phylum and OTU ----

# By phylum
# Melt the phylum table
phylum = t(d.phylum)
phylum_melt = melt(phylum)

# Compute the mean and standard deviation
mean_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = mean)
sd_abund_phylum = aggregate(value ~ Var1, data = phylum_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_phylum = cbind(mean_abund_phylum, sd = sd_abund_phylum$value)
mean_abund_phylum = mean_abund_phylum[order(mean_abund_phylum[,2], decreasing = T),]

dir.create(path = "./mean_abundances/")

write.table(mean_abund_phylum, file="./mean_abundances/agt_swt_mean_abundance_sd_phylum.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE) #For Agt_Swt
write.table(mean_abund_phylum, file="./mean_abundances/agt_mean_abundance_sd_phylum.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE) #For Agt
write.table(mean_abund_phylum, file="./mean_abundances/swt_mean_abundance_sd_phylum.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE) #For Swt



# Complete OTU table
# Empty OTUs (that were excluded based on sampling depth) should be removed
otu_summary = data.frame(mean.abundance = round(rowMeans(microbio.relative[rowSums(microbio.relative) > 0,])*100, 2), Taxonomy = microbio.taxonomy[rowSums(microbio.relative) > 0, 2]) #might need to transpose microbio.relative

#Save Table
write.table(otu_summary, file="./agt_swt_complete_otu_summary_orig.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE) #For Agt_Swt
write.table(otu_summary, file="./agt_complete_otu_summary_orig.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE) #For Agt
write.table(otu_summary, file="./swt_complete_otu_summary_orig.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE) #For Swt



# By OTUs
# Melt the OTU table
otus_melt = melt(microbio.relative)

# Compute mean and standard deviation
mean_abund_otus = aggregate(value ~ Var1, data = otus_melt, FUN = mean)
sd_abund_otus = aggregate(value ~ Var1, data = otus_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_otus = cbind(mean_abund_otus, sd = sd_abund_otus$value)
mean_abund_otus = mean_abund_otus[order(mean_abund_otus[,2], decreasing = T),]
top_ten_otus = head(mean_abund_otus, 10)

write.table(mean_abund_otus, file="./mean_abundances/agt_swt_mean_abundance_sd_otus.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)
write.table(top_ten_otus, file="./mean_abundances/agt_swt_mean_abundance_otus_top10.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)

write.table(mean_abund_otus, file="./mean_abundances/agt_mean_abundance_sd_otus.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)
write.table(top_ten_otus, file="./mean_abundances/agt_mean_abundance_otus_top10.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)

write.table(mean_abund_otus, file="./mean_abundances/swt_mean_abundance_sd_otus.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)
write.table(top_ten_otus, file="./mean_abundances/swt_mean_abundance_otus_top10.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)


# Boxplot of phyla and top OTUs
# Boxplot of phyla

theme_set(theme_classic())
# Combine phyla with very low abundance
phyla_median = aggregate(value ~ Var1, data = phylum_melt, FUN = median)
top_phyla = phyla_median$Var1[phyla_median$value > 0]
bottom_phyla = phyla_median$Var1[phyla_median$value == 0]

top_bottom_phyla = rbind(phylum[top_phyla, ], "Other" = colSums(phylum[bottom_phyla, ])) 

phylum_melt = melt(top_bottom_phyla)

phylum_melt$value[phylum_melt$value < 0.00005] = 0.00005

#Agt_Swt Labels Redone
phyla_labels = c("Bacteroidetes","Firmicutes", "Actinobacteria", "Proteobacteria", "Verrucomicrobia","Tenericutes","Unclassified*","Synergistetes", "Fusobacteria", "Elusimicrobia","Other phyla", "Unclassified**")

#Agt_Swt Labels Original
phyla_labels = c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Verrucomicrobia","Actinobacteria","Lentisphaerae","Unclassified*","Other phyla", "Tenericutes","Synergistetes")

#Agt Labels Redone
phyla_labels = c("Bacteroidetes","Firmicutes","Actinobacteria", "Proteobacteria", "Verrucomicrobia","Tenericutes","Unclassified*","Synergistetes", "Fusobacteria", "Elusimicrobia", "Other phyla")

#Agt Labels Original
phyla_labels = c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria","Verrucomicrobia","Lentisphaerae","Unclassified*","Other phyla", "Tenericutes","Synergistetes")

#Swt Labels Original
phyla_labels = c("Bacteroidetes","Firmicutes",  "Proteobacteria","Actinobacteria", "Verrucomicrobia","Lentisphaerae","Other phyla","Unclassified*", "Unclassified****")

#Swt Labels Redone
#phyla_labels = c("Bacteroidetes","Firmicutes","Actinobacteria", "Proteobacteria", "Verrucomicrobia","Tenericutes","Unclassified*","Synergistetes","Other phyla","Fusobacteria","Euryarchaeota", "Unclassified**")



box_phylum_agt_swt = ggplot(phylum_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title=element_text(size=15)) +
  scale_x_discrete(labels = phyla_labels)

box_phylum_agt = ggplot(phylum_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title=element_text(size=15)) +
  scale_x_discrete(labels = phyla_labels)

box_phylum_swt = ggplot(phylum_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title=element_text(size=15)) +
  scale_x_discrete(labels = phyla_labels)



#Save Plot
ggsave("agt_swt_top_phyla_raw.pdf", plot=last_plot())
ggsave("agt_swt_top_phyla.pdf", plot=box_phylum)

# Boxplot of OTUs
top_otus_melt = melt(microbio.relative[top_ten_otus$Var1,])

OTU_labels = c("Succinivibrio", "F. Prausnitzii","Prevotella", "Prevotella","Bacteroides","Succinivibrio", "F. Prausnitzii","Bacteroides","Prevotella", "Bacteroides vulgatus")

box_otus = ggplot(top_otus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title=element_text(size=9)) +
  scale_x_discrete(labels = OTU_labels)

ggsave("agt_swt_top_otus.pdf", plot=box_otus)

# Figure 1
combined=ggdraw() +
  draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_otus,  0.45, 0.35, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)

ggsave("agt_swt_combined_topOTUs_topPhyla_plots.pdf", plot=combined)

#High res save: ##Genus melt below! Used genus and not OTUs
png("agt_swt_combined_topGenus_topPhyla_plot.png", units="mm", width=250, height=200, res=300)
ggdraw() +
  draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus,  0.45, 0.35, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)
dev.off()

#To save in workspace:
combined_agt_swt <- ggdraw() +
  draw_plot(box_phylum_agt_swt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_agt_swt,  0.45, 0.32, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)

combined_agt <- ggdraw() +
  draw_plot(box_phylum_agt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_agt,  0.45, 0.32, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("E", "F"), c(0, 0.50), c(1, 1), size = 15)

combined_swt <- ggdraw() +
  draw_plot(box_phylum_swt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_swt,  0.45, 0.32, 0.60, 0.73, scale = 0.8) +
  draw_plot_label(c("C", "D"), c(0, 0.50), c(1, 1), size = 15)

##High res save grid:
png("fig4_top_phyla_genus_A_agt_swt_C_swt_E_agt_2.png", units="mm", width=300, height=300, res=300) 
plot_grid(combined_agt_swt, combined_swt, combined_agt, align = 'h')
dev.off()


##Test
combined_agt_swt <- ggdraw() +
  draw_plot(box_phylum_agt_swt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_agt_swt,  0.42, 0.36, 0.68, 0.73, scale = 0.75) +
  draw_plot_label(c("A", "B"), c(0, 0.49), c(1, 1), size = 15)

combined_agt <- ggdraw() +
  draw_plot(box_phylum_agt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_agt,  0.43, 0.36, 0.68, 0.73, scale = 0.75) +
  draw_plot_label(c("E", "F"), c(0, 0.50), c(1, 1), size = 15)

combined_swt <- ggdraw() +
  draw_plot(box_phylum_swt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_swt,  0.41, 0.36, 0.68, 0.73, scale = 0.75) +
  draw_plot_label(c("C", "D"), c(0, 0.48), c(1, 1), size = 15)

box_genus_swt=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
axis.title=element_text(size=12)) + scale_x_discrete(labels = genus_labels)

box_genus_agt_swt=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
axis.title=element_text(size=12)) + scale_x_discrete(labels = genus_labels)

box_genus_agt=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
axis.title=element_text(size=12)) + scale_x_discrete(labels = genus_labels)



# By genus
# Melt the genus table
genus = t(d.genus)
genus_melt = melt(genus)

# Compute the mean and standard deviation
mean_abund_genus = aggregate(value ~ Var1, data = genus_melt, FUN = mean)
sd_abund_genus = aggregate(value ~ Var1, data = genus_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_genus = cbind(mean_abund_genus, sd = sd_abund_genus$value)
mean_abund_genus = mean_abund_genus[order(mean_abund_genus[,2], decreasing = T),]

top_ten_genus = head(mean_abund_genus, 10)
#top_genus_melt = melt(microbio.relative[top_ten_genus$Var1,])
top_genus_melt = melt(top_ten_genus) 

write.table(mean_abund_genus, file="./mean_abundances/agt_swt_mean_abundance_sd_genus.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)
write.table(top_ten_genus, file="./mean_abundances/agt_swt_mean_abundance_genus_top10.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)

write.table(mean_abund_genus, file="./mean_abundances/agt_mean_abundance_sd_genus.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)
write.table(top_ten_genus, file="./mean_abundances/agt_mean_abundance_genus_top10.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)

write.table(mean_abund_genus, file="./mean_abundances/swt_mean_abundance_sd_genus.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)
write.table(top_ten_genus, file="./mean_abundances/swt_mean_abundance_genus_top10.txt", col.names=TRUE,row.names=FALSE,sep="\t", quote=FALSE)

##Agt_Swt Genus Labels Redone
genus_labels = c("Alistipes", "Acidaminococcus", "Bacteroides", "Actinomyces", "Alloprevotella", "Acetanaerobacterium", "Abiotrophia", "Akkermansia", "Allisonella", "Butyricimonas")


##Agt_Swt Genus Labels
genus_labels = c("Prevotella", "Bacteroides", "Faecalibacterium","Ruminococcaceae**","Succinivibrio","Clostridiales***", "Oscillibacter", "Clostridium_XIVa", "Roseburia", "Ruminococcus")

##Agt Genus Labels Redone
genus_labels = c("Alistipes", "Bacteroides", "Acidaminococcus", "Actinomyces","Abiotrophia", "Acetanaerobacterium", "Alloprevotella", "Akkermansia", "Allisonella", "Butyricimonas")

##Agt Genus Labels
genus_labels = c("Prevotella", "Succinivibrio", "Bacteroides", "Faecalibacterium", "Ruminococcaceae**","Clostridiales***", "Oscillibacter", "Clostridium_XIVa", "Ruminococcus", "Roseburia")

##Swt Genus Labels Redone
#genus_labels = c("Acidaminococcus","Alistipes","Actinomyces", "Bacteroides","Alloprevotella","Acetanaerobacterium", "Akkermansia", "Achromobacter","Butyricimonas", "Allisonella")

##Swt Genus Labels

genus_labels = c("Prevotella", "Bacteroides", "Faecalibacterium", "Ruminococcaceae**","Clostridiales***", "Succinivibrio", "Alistipes", "Oscillibacter", "Clostridium_XIVa", "Roseburia")



box_genus=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
axis.title=element_text(size=14)) + scale_x_discrete(labels = genus_labels)

box_genus_agt_swt=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
axis.title=element_text(size=14)) + scale_x_discrete(labels = genus_labels)

box_genus_swt=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
axis.title=element_text(size=14)) + scale_x_discrete(labels = genus_labels)

box_genus_agt=ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() +
  labs(x = "", y = "Relative abundance") +
  theme(axis.text.x = element_text(face = 'italic', angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
axis.title=element_text(size=14)) + scale_x_discrete(labels = genus_labels)


ggsave("agt_swt_top_genus_raw.pdf", plot=last_plot())
ggsave("agt_swt_top_genus.pdf", plot=box_genus)

##Figure 1 alternative(genus)

combined_alt=ggdraw() +
  draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus,  0.35, 0.35, 0.55, 0.73, scale = 0.8) +
  draw_plot_label(c("A", "B"), c(0, 0.40), c(1, 1), size = 15)
 

combined_alt_swt=ggdraw() +
  draw_plot(box_phylum_swt, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus_swt,  0.35, 0.35, 0.55, 0.73, scale = 0.8) +
  draw_plot_label(c("C", "D"), c(0, 0.40), c(1, 1), size = 15)
 


#ggdraw() + 
#  draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
#  draw_plot(box_genus,  0.4, 0.5, 0.6, 0.7, scale = 0.6) +
#  draw_plot_label(c("A", "B"), c(0, 0.40), c(1, 1), size = 15)

##For Swt:
combined_alt=ggdraw() +
draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
  draw_plot(box_genus,  0.45, 0.5, 0.60, 0.6, scale = 0.7) +
  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)


#combined_alt=ggdraw() +
#draw_plot(box_phylum, 0, 0, 1, 1, scale = 1) +
#  draw_plot(box_genus,  0.45, 0.35, 0.60, 0.73, scale = 0.8) +
#  draw_plot_label(c("A", "B"), c(0, 0.50), c(1, 1), size = 15)
#

ggsave("agt_swt_combined_topGenus_topPhyla_plots.pdf", plot=combined_alt)                                                                                       






======================================================
### Calculating UniFrac distance matrices ----
##Setwd() as needed
#For Agt:
setwd("/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/july_analyses/enterotype/agt/enterogradient_agt")

#For Agt:
setwd("/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2019_Analyses/july_analyses/enterotype/swt/enterogradient_swt")

unifracs <- GUniFrac(microbio.rare, microbio.tree, alpha=c(0, 0.5, 1))$unifracs

dw <- unifracs[, , "d_1"]   # Weighted UniFrac
du <- unifracs[, , "d_UW"]  # Unweighted UniFrac
#d5 <- unifracs[, , "d_0.5"] # Generalized UniFrac with alpha 0.5

# Enterotype analysis ----
# Clustering analysis and evaluation of the enterotyping
# Modified from Arumugam et al., 2011 (doi:10.1038/nature09944)
# Code available at enterotyping.embl.de

# K-medoids function
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

# Weighted UniFrac for use with enterotypes
data.dist = as.dist(dw)

# PCoA analysis and computation of the proportion of the variance of each axis
e.pcoa = cmdscale(data.dist, k=5, eig = T)
e.PC1 = round(e.pcoa$eig[1]/sum(e.pcoa$eig), 4)* 100
e.PC2 = round(e.pcoa$eig[2]/sum(e.pcoa$eig), 4)* 100
e.PC3 = round(e.pcoa$eig[3]/sum(e.pcoa$eig), 4)* 100


# Optimum number of clusters using the Silhouette index instead of Calinski-Harabasz index
nclusters=NULL

for (k in 1:20) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k] = mean(silhouette(data.cluster_temp, data.dist)[,3])
    #nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist, centrotypes = "medoids")
  }
}

# Grouping with k clusters
data.cluster=pam.clustering(data.dist, k=2)

# SI index for the above clustering
# -1 <= S(i) <= 1
# A sample which is closer to its own cluster than to any other cluster has a high S(i),
# S(i) close to 0 implies that the given sample lies somewhere between two clusters.
# Large negative S(i) values indicate that the sample was assigned to the wrong cluster.
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
cat(obs.silhouette) #For agt_swt:  k=2, S(i) = 0.2462824 (previous = 0.2529376); redone in Dec2019 S(i) = 0.251652 for dw;S(i)=0.1534415 (previous=0.150321) for du. For Agt: k=2, S(i) = 0.2560603. For Swt: k=2, S(i) = 0.2277168


# plots SI index ~ number of clusters
ggplot(data.frame(n = factor(2:20), nclusters[-1]), aes(n, nclusters[-1])) + geom_bar(stat = "identity") + labs(x = "Number of clusters", y = "Average silhouette index")

ggsave("agt_silhouettes.pdf", plot =last_plot())

##High res save:
png("agt_swt_silhouettes_weighted.png", units="mm", width=250, height=200, res=300)
ggplot(data.frame(n = factor(2:20), nclusters[-1]), aes(n, nclusters[-1])) + geom_bar(stat = "identity") + labs(x = "Number of clusters", y = "Average silhouette index")
dev.off()


png("agt_swt_silhouettes_weighted_left_unweighted_right.png", units="mm", width=250, height=200, res=300)
plot_grid(weighted_silhouette_index_agt_swt, unweighted_silhouette_index_agt_swt, labels = c('A', 'B'), align = 'h')
dev.off()

png("agt_swt_silhouettes_weighted.png", units="mm", width=150, height=100, res=300)
ggplot(data.frame(n = factor(2:20), nclusters[-1]), aes(n, nclusters[-1])) + geom_bar(stat = "identity") + labs(x = "Number of clusters", y = "Average silhouette index") + theme(axis.text.x = element_text(size = 12))
dev.off()

# PCoA graph of the enterotyping with k clusters
clusters = as.factor(data.cluster)
w = as.data.frame(cmdscale(dw, k = 2))
ggplot(w, aes(V1, V2)) + geom_point(aes(color = clusters)) +
  stat_ellipse(type = "t", mapping = aes(color = clusters), level = 0.75) +
  labs(x = paste("PC1",e.PC1, "%"),
       y = paste("PC2",e.PC2, "%"),
title = "PCoA Weighted UniFrac")

ggsave("agt_swt_enterotype_cluster_k2.pdf", plot=last_plot()) #For k=2

###For High res save:
w = as.data.frame(cmdscale(dw, k = 2))
ggplot(w, aes(V1, V2)) + geom_point(aes(color = clusters)) + stat_ellipse(type = "t", mapping = aes(color = clusters), level = 0.75) +
  labs(x = paste("PC1",e.PC1, "%"),
       y = paste("PC2",e.PC2, "%"),
title = "PCoA Weighted UniFrac") + theme(axis.text.x = element_text(size = 12)

png("SF3_agt_swt_PCoA_weighted_uniFrac.png", units="mm", width=150, height=100, res=300)
ggplot(w, aes(V1, V2)) + geom_point(aes(color = clusters)) +
  stat_ellipse(type = "t", mapping = aes(color = clusters), level = 0.75) +
  labs(x = paste("PC1",e.PC1, "%"),
       y = paste("PC2",e.PC2, "%"),
title = "PCoA Weighted UniFrac")
dev.off()


#To identify which samples are clustered together

d<-data.frame(dw) #format distance matrix as df
cl<-data.frame(d[1], clusters ) #Merge clusters to df to annotate samples with cluster information.

#Save as .csv file
write.csv(cl, "agt_swt_cluster_assignment_weighted_uniFrac.csv")

================================================================================
# Enterogradient analysis ----
# Selection of OTUs with median relative abundance >= 0.01%
microbio.relative_tax = cbind(microbio.taxonomy[,2], microbio.relative)

# Compute the median of each OTU and filter
median_abund = apply(microbio.relative, MARGIN = 1, FUN = median)
abundant_otus = microbio.relative[median_abund >= 0.0001, ] #Agt_Swt: 124 ASVs; Agt: 151 ASVs; Swt: 
nrow(abundant_otus)

# Descriptive statistics of the subset of most abundant OTUs
mean(colSums(abundant_otus))
sd(colSums(abundant_otus))
range(range(colSums(abundant_otus)))
summary(colSums(abundant_otus))
rownames(abundant_otus)

# Correlation analysis between the first 3 axis of weighted UniFrac PCoA and the most abundant OTUs
# For each axis
# 1. Spearman's correlation tests are done between each axis and OTU
# 2. A table is generated with rho and p-value of each result of point 1
# 3. p-value correction is made using qvalue
# 4. Significant correlations are selected if q<0.05
# 5. A table is generated with the taxonomy, rho, p-value and q-value of the significant OTUs

# PC1
cor_entero_PC1 = apply(X = abundant_otus, MARGIN = 1, 
                       FUN = function(x) cor.test(x, e.pcoa$points[, 1], m = "s"))

cor_PC1 = sapply(X = cor_entero_PC1,
                 FUN = function(x) cbind(rho = as.vector(x$estimate),
                                         p =as.vector(x$p.value)))

cor_PC1 = as.data.frame(t(cor_PC1))
colnames(cor_PC1) = c("rho", "p")
q_value = as.vector(qvalue(cor_PC1$p)$qvalue) #First try this before next line.
##q_value = as.vector(qvalue(cor_PC1$p,pi0=1)$qvalue) #Use this line in case of an error in previous line
cor_PC1$q = q_value
sig_cor_PC1 = cor_PC1[cor_PC1$q <= 0.05,] #Agt_Swt: 54 ASVs; Agt: 68 ASVs; Swt: 21 ASVs
result_cor_PC1 = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(sig_cor_PC1),2], sig_cor_PC1)
result_cor_PC1 = result_cor_PC1[abs(result_cor_PC1$rho) >= 0.3,] #Agt_Swt: 27 ASVs; Agt: 29 ASVs; Swt:21 ASVs
sig_otus_PC1 = abundant_otus[cor_PC1[,3] <= 0.05,]
nrow(result_cor_PC1) #Agt_Swt: 27 ASVs; Agt: 29 ASVs; Swt:21 ASVs


#PC2
cor_entero_PC2 = apply(X = abundant_otus, MARGIN = 1, 
                       FUN = function(x) cor.test(x, e.pcoa$points[, 2], m = "s"))
cor_PC2 = sapply(X = cor_entero_PC2, 
                 FUN = function(x) cbind(rho = as.vector(x$estimate), 
                                         p =as.vector(x$p.value)))

cor_PC2 = as.data.frame(t(cor_PC2))
colnames(cor_PC2) = c("rho", "p")
q_value = as.vector(qvalue(cor_PC2$p)$qvalue)
cor_PC2$q = q_value
sig_cor_PC2 = cor_PC2[cor_PC2$q <= 0.05,] #Agt_Swt: 82 ASVs; Agt: 68 ASVs; Swt:42 ASVs.
result_cor_PC2 = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(sig_cor_PC2),2], sig_cor_PC2)
result_cor_PC2 = result_cor_PC2[abs(result_cor_PC2$rho) >= 0.3,] #Agt_Swt: 36 ASVs; Agt: 19 ASvs; Swt:30 ASVs.
sig_otus_PC2 = abundant_otus[cor_PC2[,3] <= 0.05,] 
nrow(result_cor_PC2) #Agt_Swt: 36 ASVs; Agt: 19 ASVs; Swt: 30 ASVs.

# PC3
cor_entero_PC3 = apply(X = abundant_otus, MARGIN = 1, 
                       FUN = function(x) cor.test(x, e.pcoa$points[, 3],
                                                  m = "s"))
cor_PC3 = sapply(X = cor_entero_PC3, 
                 FUN = function(x) cbind(rho = as.vector(x$estimate),
                                         p =as.vector(x$p.value)))

cor_PC3 = as.data.frame(t(cor_PC3))
colnames(cor_PC3) = c("rho", "p")
q_value = as.vector(qvalue(cor_PC3$p)$qvalue)
##q_value = as.vector(qvalue(cor_PC1$p,pi0=1)$qvalue)
cor_PC3$q = q_value
sig_cor_PC3 = cor_PC3[cor_PC3$q <= 0.05,] #Agt_Swt: 54 ASVs; Agt: 57 ASVs; Swt:
result_cor_PC3 = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(sig_cor_PC3),2], sig_cor_PC3)
result_cor_PC3 = result_cor_PC3[abs(result_cor_PC3$rho) >= 0.3,]
sig_otus_PC3 = abundant_otus[cor_PC3[,3] <= 0.05,] #Agt_Swt: 5 ASVs; Agt: 10 ASVs; Swt:
nrow(result_cor_PC3) #Agt_Swt: 5 ASVs; Agt: 10 ASVs; Swt: 

# Table with the OTUs correlated with at least one PCoA axis 
all_axes_otus_names = sort(unique(c(rownames(result_cor_PC1), 
                                    rownames(result_cor_PC2), 
                                    rownames(result_cor_PC3))))

all_axes_otus = abundant_otus[all_axes_otus_names,]
cor_all_PC1 = apply(X = all_axes_otus, MARGIN = 1, 
                    FUN = function(x) cor(x, e.pcoa$points[, 1], m = "s"))

cor_all_PC2 = apply(X = all_axes_otus, MARGIN = 1, 
                    FUN = function(x) cor(x, e.pcoa$points[, 2], m = "s"))

cor_all_PC3 = apply(X = all_axes_otus, MARGIN = 1, 
                    FUN = function(x) cor(x, e.pcoa$points[, 3], m = "s"))

result_cor_all = data.frame(Tax = microbio.taxonomy[row.names(microbio.taxonomy) %in% row.names(all_axes_otus),2], 
                            round(cor_all_PC1, 2), 
                            round(cor_all_PC2, 2), 
                            round(cor_all_PC3, 2))

#Save Table
write.table(result_cor_all, file="./agt_swt_ASVs_correlated_with_PCOA_axes_weighted_unifrac_enterogradient_analysis.txt", col.names=TRUE,row.names=TRUE,sep="\t", quote=FALSE)

================================================================================
# Select phylotypes characteristic of Western and non-Western microbiota ----

bacteroides = d.genus[,"Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides"]
rownames(bacteroides) = NULL

prevotella = d.genus[,"Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella"]
rownames(prevotella) = NULL

ruminococcus = d.genus[,"Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcus"]
rownames(ruminococcus) = NULL

treponema = d.genus[,"Bacteria;Spirochaetes;Spirochaetia;Spirochaetales;Spirochaetaceae;Treponema"]
rownames(treponema) = NULL

bifidobacterium = d.genus[,"Bacteria;Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium"]
rownames(bifidobacterium) = NULL

cytophagales = d.order[,"Bacteria;Bacteroidetes;Cytophagia;Cytophagales"]
rownames(cytophagales) = NULL

barnesiella = d.genus[,"Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Barnesiella"]
rownames(barnesiella) = NULL


select_entero_table = data.frame(ID = rownames(d.genus),
                                 Prevotella = prevotella,
                                 Bacteroides = bacteroides,
                                 Ruminococcus = ruminococcus,
                                 Treponema = treponema,
                                 Bifidobacterium = bifidobacterium,
                                 Cytophagales = cytophagales,
                                 Barnesiella = barnesiella,
                                 PC1 = e.pcoa$points[, 1],
                                 PC2 = e.pcoa$points[, 2],
                                 PC3 = e.pcoa$points[, 3])

# Test that the mean phylotype abundance is different from that found in the
# meta-analysis of benchmark datasets from curatedMetagenomicData
# (run the accompanying "curated_metagenomes.R" script from the Escobar paper to get these values)
##Used pre-calculated mu values for the Escobar paper

t_prevotella=t.test(prevotella, mu = 0.245)
t_bifidobacterium = t.test(bifidobacterium, mu = 0.074)
t_bacteroides=t.test(bacteroides, mu = 0.230)
t_treponema=t.test(treponema, mu = 0.021)
t_cytophagales=t.test(cytophagales, mu = 6e-05)
t_barnsiella=t.test(barnesiella, mu = 0.0126)

# Test that the mean phylotype abundance is different from 0
t_bacteroides_0=t.test(bacteroides, mu = 0, alternative = "g")
t_prevotella_0=t.test(prevotella, mu = 0, alternative = "g")
t_bifidobacterium_0 = t.test(bifidobacterium, mu = 0, alternative = "g")
t_treponema_0=t.test(treponema, mu = 0, alternative = "g")
t_cytophagales_0=t.test(cytophagales, mu = 0, alternative = "g")
t_barnsiella_0=t.test(barnesiella, mu = 0, alternative = "g")


## Western - non-Western histograms (Figure 2)
hist_bacteroides = ggplot(select_entero_table, aes(Bacteroides)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Bacteroides abundance", y = "Number of subjects")

 hist_bacteroides = ggplot(select_entero_table, aes(Bacteroides)) +
+   geom_histogram(boundary = 0, bins = 100) + theme_gray() +
+   labs(x = expression(~italic(Bacteroides)~ "abundance"), y = "Number of individuals")


hist_prevotella = ggplot(select_entero_table, aes(Prevotella)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Prevotella abundance", y = "Number of subjects")

hist_prevotella = ggplot(select_entero_table, aes(Prevotella)) +
+   geom_histogram(boundary = 0, bins = 100) + theme_gray() +
+   labs(x = expression(~italic(Prevotella)~ "abundance"), y = "Number of individuals")


hist_bifidobacterium = ggplot(select_entero_table, aes(Bifidobacterium)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Bifidobacterium abundance", y = "Number of subjects")

hist_bifidobacterium = ggplot(select_entero_table, aes(Bifidobacterium)) +
+   geom_histogram(boundary = 0, bins = 100) + theme_gray() +
+   labs(x = expression(~italic(Bifidobacterium)~ "abundance"), y = "Number of individuals")


hist_treponema = ggplot(select_entero_table, aes(Treponema)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Treponema abundance", y = "Number of subjects")
grob <- grobTree(textGrob("#165", x=0.06,  y=0.98, hjust=0, gp=gpar(fontsize=11))) #153 subjects with zero abundance values
hist_treponema_modified <-  hist_treponema +  coord_cartesian(ylim = c(0, 10)) + annotation_custom(grob) #requires library(grid) #modified y-axis

hist_treponema = ggplot(select_entero_table, aes(Treponema)) +
+   geom_histogram(boundary = 0, bins = 100) + theme_gray() +
+   labs(x = expression(~italic(Treponema)~ "abundance"), y = "Number of individuals")
> grob <- grobTree(textGrob("#165", x=0.06,  y=0.98, hjust=0, gp=gpar(fontsize=11))) #153 individuals with zero abundance values
> hist_treponema_modified <-  hist_treponema +  coord_cartesian(ylim = c(0, 10)) + annotation_custom(grob) #requires library(grid) #modified y-axis




hist_cytophagales = ggplot(select_entero_table, aes(Cytophagales)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Cytophagales abundance", y = "Number of subjects")

hist_barnesiella = ggplot(select_entero_table, aes(Barnesiella)) +
  geom_histogram(boundary = 0, bins = 100) + theme_gray() +
  labs(x = "Barnesiella abundance", y = "Number of subjects")
hist_barnesiella = ggplot(select_entero_table, aes(Barnesiella)) +
+   geom_histogram(boundary = 0, bins = 100) + theme_gray() +
+   labs(x = expression(~italic(Barnesiella)~ "abundance"), y = "Number of individuals")



hist_grid_agt_swt=plot_grid(hist_bacteroides, hist_bifidobacterium,
          hist_barnesiella, hist_prevotella, hist_treponema,
          labels="AUTO", nrow = 2, ncol = 3)

ggsave("swt_hist_grid_western_non_western_big.pdf", plot =last_plot())

#High res save:
png("fig7_agt_swt_hist_grid_western_non_western.png", units="mm", width=250, height=200, res=300)
plot_grid(hist_bacteroides, hist_bifidobacterium,
          hist_barnesiella, hist_prevotella, hist_treponema,
          labels="AUTO", nrow = 2, ncol = 3)
dev.off()

#Modified:
png("fig7_agt_swt_hist_grid_western_non_western.png", units="mm", width=250, height=200, res=300)
plot_grid(hist_bacteroides, hist_bifidobacterium,
          hist_barnesiella, hist_prevotella, hist_treponema_modified,
          labels="AUTO", nrow = 2, ncol = 3)
dev.off()


# Prevotella-Bacteroides ratio in the PCoA of weighted UniFrac (Figure SR1)
prev_bact_qplot=qplot(PC1, PC2, data = select_entero_table,
      colour = bacteroides/(bacteroides+prevotella) ) +
  scale_colour_gradientn(colours=(rainbow(10)), name = "") +
  labs(x = "PCo1 36.30%", y = "PCo2 16.98%")

ggsave("agt_qplot_prevotella_bacteroides_ratio.pdf", plot =prev_bact_qplot)

prev_bact_cortest=cor.test(prevotella, bacteroides, method = "s")

#Save workspace
save.image(file = "agt_enterogradient_july2019.RData")











##Evaluating htn across sites controlling for batch effects
otu_table(ps_trans) <- otu_table(obj_trans, taxa_are_rows = TRUE)


##DESEq2
diagdds = phyloseq_to_deseq2(filt_agt_lnob, ~ batch + bmi_info) #Agt lean vs obese comparisons whilst controlling for potential batch effects
diagdds = phyloseq_to_deseq2(filt_swt_lnob, ~ bmi_info) #Swt lean vs obese comparisons - sequenced in one batch
diagdds = phyloseq_to_deseq2(filt_agt_swt_lnob, ~ batch + bmi_info) #Lean vs obese comparisons of pilot cohort excluding overweight samples whilst controlling for batch effects
diagdds = phyloseq_to_deseq2(filt_agt_swt, ~ batch + bmi_info + site) #Agt vs Soweto comparisons of all data including overweight samples
diagdds = phyloseq_to_deseq2(filt_agt_swt_ob, ~ bmi_info) #Comparing compositional differences between sites in only obese samples
diagdds = phyloseq_to_deseq2(filt_agt_swt_ln, ~ batch + bmi_info) #Comparing compositional differences between sites in only lean samples whilst controlling for potential batch effects


##Controlling for htn
diagdds = phyloseq_to_deseq2(filt_agt_swt_lnob, ~ batch + htn + bmi_info) #Lean vs obese comparisons of pilot cohort excluding overweight samples whilst controlling for batch effects
diagdds = phyloseq_to_deseq2(filt_agt_swt, ~ batch + bmi_info + htn + site) #Agt vs Soweto comparisons of all data including overweight samples

##Htn as the outcome variable whilst controlling for BMI as a continuous variable using log10 BMI values
diagdds = phyloseq_to_deseq2(filt_agt_lnob, ~ batch + log_bmi + htn)
diagdds = phyloseq_to_deseq2(filt_swt_lnob, ~ log_bmi + htn)
diagdds = phyloseq_to_deseq2(filt_agt_swt, ~ batch + log_bmi + site)


#diagdds = phyloseq_to_deseq2(swted, ~ bmi_info + htn)# account (control) for differences due to bmi_info whilst estimating the effect due to htn. Condition of interest should always be the last term - covariates should come before.

#For agt_swt site analyses whilst controlling for batch effects:
diagdds = phyloseq_to_deseq2(filt_agt_swt_05, ~ batch + site)

#Controlling for both site and bmi_info:
diagdds = phyloseq_to_deseq2(filt_agt_swt_05, ~ batch + bmi_info_2 + site)
#For Swt and Agt Lean vs Obese Comparisons:
diagdds = phyloseq_to_deseq2(agt_obln, ~ bmi_info) #No batch control because one batch for all samples in Swt; and for Agt, batch categorizatio is pretty much equivalent to lean/obese categories.

#For Agt_Swt site comparison of lean and obese
diagdds = phyloseq_to_deseq2(agt_swt_ob, ~ site) #all the obese samples were done in one batch; controlling for batch in the lean samples is equivalent to site compparisons, so no batch comparisons were done.

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
##Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(swted)[rownames(sigtab),], "matrix"))
head(sigtab)

##First, cleaning up the table a little for legibility.
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

##Writing Deseq Report
write.csv(as.data.frame(sigtab), file="condition_treated_results.csv")
#write.csv(as.data.frame(sigtab), file="agt_swt_site_batch_bmi_controlled_filt_05.csv")

##Plot Results
sigtabgen<-sigtab
theme_set(theme_bw())

# Genus order <-this was used
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

####Text size adjustment in plot:
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=1) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), axis.text.y=element_text(size = 7.5))



















extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,'\t')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}

otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa)){
    taxa = colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length as the number of columns in the OTU table")
    return;
  }
level.names = sapply(as.character(taxa),
                       function(x)
                         extract.name.level(x,level=level))
  t(apply(x, 1,
          function(y)
            tapply(y,level.names,sum)))
}






# By genus
# Melt the genus table
genus = t(d.genus)
genus_melt = melt(genus)

# Compute the mean and standard deviation
mean_abund_genus = aggregate(value ~ Var1, data = genus_melt, FUN = mean)
sd_abund_genus = aggregate(value ~ Var1, data = genus_melt, FUN = sd)

# Ordered mean and SD table
mean_abund_genus = cbind(mean_abund_genus, sd = sd_abund_genus$value)
mean_abund_genus = mean_abund_genus[order(mean_abund_genus[,2], decreasing = T),]
top_ten_genus = head(mean_abund_genus, 10)
#top_genus_melt = melt(microbio.relative[top_ten_genus$Var1,])
top_genus_melt = melt(top_ten_genus$Var1)

ggplot(top_genus_melt, aes(x = reorder(Var1, -value, FUN=median), y = value)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative abundance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 8), 
axis.title=element_text(size=9)) + scale_x_discrete(labels = genus_labels)

genus_labels = c("Prevotella", "Bacteroides", "Succinivibrio", "Faecalibacterium", "Ruminococcaceae","Clostridiales", "Oscillibacter", "Clostridium_XIVa", "Roseburia", "Ruminococcus")











#Diversity Measurements
##alpha diversity table - microbiome package

tab<-diversities(agt_swt.rarefied, index="all")

##Visualize tab:
kable(head(tab))

###Adding diversity data to meta table
agt_swt.rarefied.meta<-meta(agt_swt.rarefied)
agt_swt.rarefied.meta$Shannon <-tab$shannon
agt_swt.rarefied.meta$Inverse_Simpson <-tab$inverse_simpson
bmi<-levels(agt_swt.rarefied.meta$bmi_info)

#Save file
write_phyloseq(agt_swt.rarefied.meta, type = "all", path = getwd())
write_phyloseq(tab, type = "all", path = getwd())

##Subset samples to be site-specific and exclude overweight samples when necessary

agt<-subset_samples(agt_swt_sub, site != "Soweto")

##Alpha Diversity 

##Estimate Richness
rich = estimate_richness(ps.rarefied)
rich_agt = estimate_richness(agt_lnob.rarefied) #agt lean and obese samples only

##Test whether the observed number of OTUs differs significantly between categories. We make a non-parametric test, the Wilcoxon rank-sum test (Mann-Whitney):

pairwise.wilcox.test(rich$Observed, sample_data(agt_lnob.rarefied)$bmi_info)

##By default, the function pairwise.wilcox.test() reports the pairwise adjusted (Holm) p-values.

##Plot Richness (Phyloseq)
##Between sites:
plot_richness(agt_swt.rarefied, x="site", measures=c("Observed", "Shannon"), color ="bmi_info") + geom_boxplot()

##ggsave()

##Categories within sites:
plot_richness(agt.rarefied, x="bmi_info", measures=c("Observed", "Shannon", "Chao1"), color ="bmi_info") + geom_boxplot()

plot_richness(agt.rarefied, x="bmi_info", measures=c("Shannon"), color ="bmi_info") + geom_boxplot()

plot_richness(swt.rarefied, x="bmi_info", measures=c("Observed", "Shannon", "Chao1"), color ="bmi_info") + geom_boxplot()

p <- plot_anova_diversity(agt_lnob.rarefied, method = c("richness", "simpson", "shannon"), grouping_column = "bmi_info", pValueCutoff = 0.05)
