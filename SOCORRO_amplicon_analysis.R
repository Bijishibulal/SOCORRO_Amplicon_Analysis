if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
BiocManager::install("DECIPHER")
apackageVersion("dada2")
packageVersion("DECIPHER")


list.files()

samples <- scan("samples", what="character")

forward_reads <- paste0(samples, "_1.fastq.gz")
reverse_reads <- paste0(samples, "_2.fastq.gz")
forward_reads
reverse_reads

filtered_forward_reads <- paste0(samples, "_1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_2_filtered.fq.gz")

plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, 
                              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                              compress=TRUE, multithread= FALSE, truncLen=c(240,170))

list.files()
class(filtered_out)
dim(filtered_out)
filtered_out


plotQualityProfile(filtered_forward_reads)
plotQualityProfile(filtered_reverse_reads)

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)


dada_forward <- dada(filtered_forward_reads, err=err_forward_reads, pool="pseudo", multithread=TRUE)
dada_reverse <- dada(filtered_reverse_reads, err=err_reverse_reads, pool="pseudo", multithread=TRUE)


dada_forward[[1]]
dada_reverse[[1]]

merged_amplicons <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse,
                               filtered_reverse_reads, trimOverhang=TRUE, verbose = TRUE)


merged_amplicons <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse,
                               filtered_reverse_reads, trimOverhang=TRUE, justConcatenate=TRUE, verbose = TRUE)

head(merged_amplicons)

seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)
num_chim_removed <- 1 - (sum(seqtab.nochim)/sum(seqtab))
num_chim_removed


getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1], filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN), dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN), nonchim=rowSums(seqtab.nochim), final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
summary_tab

write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)


download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")
load("SILVA_SSU_r138_2019.RData")
dna <- DNAStringSet(getSequences(seqtab.nochim))
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

###renamed the sample names from .filtered.gz to only sample name in "asv_counts.tsv"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)


library(decontam)
packageVersion("decontam")

colnames(asv_tab)
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 5))
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)
table(contam_df$contaminant)
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
asv_tax[row.names(asv_tax) %in% contam_asvs, ]

library(phyloseq)
library(DESeq2)
library(vegan)
library(ggplot2)
library(dendextend)
library(tidyr)
library(viridis)
library(reshape)
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")

count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info_tab <- read.table(file = "clipboard", sep = "\t", header=TRUE, row.names =1, check.names = F)
sample_info_tab
sample_info_tab$Color <- as.character(sample_info_tab$Color)

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~Source) 
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
##enlarge the plot area + using full screen
plot(euc_clust) 
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$Color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.")
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
sample_names(vst_count_phy)
sample_names(sample_info_tab_phy)
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues
eigen_vals
write.table(eigen_vals, "PCoA_Eigen_vals.tsv", sep = "\t", quote=F, col.names=NA)

sample_info_tab
sample_info_tab$Color
plot_ordination(vst_physeq, vst_pcoa, color="Source") + 
  geom_point(size=4) + labs(col="Color") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)])) + 
  theme_bw() + theme(legend.position="none")
rarecurve(t(count_tab), step=100, lwd=2, ylab="ASVs", label=F, col = sample_info_tab$Color)
abline(v=(min(rowSums(t(count_tab)))))
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
plot_richness(ASV_physeq, x="Source", color="Source", measures=c("Chao1", "Shannon" , "Simpson")) + scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)])) + theme(legend.title = element_blank())

phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="phylum")) 
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="phylum"))[,2])
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="class"))
phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])
rownames(class_counts_tab) <- as.vector(class_tax_tab$class) 
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ]
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)
major_taxa_counts_tab
write.table(major_taxa_counts_tab, "major_taxa_counts_tab.tsv", quote=FALSE, sep="\t", col.names=NA)
identical(colSums(major_taxa_counts_tab), colSums(count_tab))
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
dim(major_taxa_proportions_tab)
major_taxa_proportions_tab
write.table(major_taxa_proportions_tab, "major_taxa_proportions_tab.tsv", quote=FALSE, sep="\t", col.names=NA)
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 1, ])
dim(temp_filt_major_taxa_proportions_tab) 
temp_filt_major_taxa_proportions_tab
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)
filt_major_taxa_proportions_tab
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab
filt_major_taxa_proportions_tab_for_plot
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa
filt_major_taxa_proportions_tab_for_plot.g <- gather(filt_major_taxa_proportions_tab_for_plot, Sample, Proportion, -Major_Taxa)
filt_major_taxa_proportions_tab_for_plot.g
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)

sample_info_tab
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab),  "color"=sample_info_tab$Color, "char"= sample_info_tab$Source, stringsAsFactors = F)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)
filt_major_taxa_proportions_tab_for_plot.g2

library(ggplot2)
library(RColorBrewer)
library(colorRamps)
filt_major_taxa_proportions_tab_for_plot.g2
colourCount <- length(unique(filt_major_taxa_proportions_tab_for_plot.g2$Major_Taxa))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples") +  scale_fill_manual(values = colorRampPalette(brewer.pal(16, "Set1"))(colourCount))

filt_major_taxa_proportions_tab_for_plot.g2
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(source), shape=factor(source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")

filt_major_taxa_proportions_tab_for_plot.g3 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)
filt_major_taxa_proportions_tab_for_plot.g3
sediment_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$Source == "Sediment"] 
sediment_sample_IDs

filt_major_taxa_proportions_sediments_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g3[filt_major_taxa_proportions_tab_for_plot.g3$Sample %in% sediment_sample_IDs, ] 
ggplot(filt_major_taxa_proportions_sediments_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(source), shape=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_sediments_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_sediments_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Sediments only")

corrosion_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$Source == "Corrosion"] 
corrosion_sample_IDs

filt_major_taxa_proportions_corrosion_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g3[filt_major_taxa_proportions_tab_for_plot.g3$Sample %in% corrosion_sample_IDs, ] 
ggplot(filt_major_taxa_proportions_corrosion_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(source), shape=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_corrosion_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_corrosion_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Corrosion only")

Seawater_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$Source == "Seawater"]
Seawater_sample_IDs
filt_major_taxa_proportions_sw_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g3[filt_major_taxa_proportions_tab_for_plot.g3$Sample %in% Seawater_sample_IDs, ]
ggplot(filt_major_taxa_proportions_sw_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,70)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(source), shape=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_sw_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_sw_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Seawater only")


Anova_result <- anova(betadisper(euc_dist, sample_info_tab$Source))
Anova_result
write.table(Anova_result, "Anova_result.xls", quote=FALSE, sep="\t", col.names=NA)


seqtab.nochim
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
install.packages("devtools")
BiocManager::install("Biostrings")
library(devtools)
install_github("liamrevell/phytools")
install_github("KlausVigo/phangorn")
library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA") 
dm <- dist.ml(phang.align)
treeUPGMA <- upgma(dm)
plot(treeUPGMA)

treeNJ <- NJ(dm) # Note, tip order != sequence order
plot(treeNJ)

fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitGTR)
saveRDS(fitGTR, "phangorn.tree.RDS")


sample_info_tab1 <- read.table(file = "clipboard", sep = "\t", header=TRUE, row.names =1, check.names = F)
sample_info_tab$Color <- as.character(sample_info_tab$Color)
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sample_info_tab1), tax_table(taxa), phy_tree(fitGTR$tree))


#Diversity
##1. Performed a MDS with euclidean distance (mathematically equivalent to a PCA)

ord <- ordinate(physeq, 'MDS', 'euclidean')
plot_ordination(physeq, ord, type='Source', color='season', title='PCA of the samples') + theme_minimal()

##2. Performed with Bray-Curtis distance

ord <- ordinate(physeq, 'NMDS', 'bray')


stressplot(ord)## Assess goodness of ordination fit (stress plot)
plot_ordination(physeq, ord, type='Samples', color='Source', title='Bray-Curtis distance/dissimilarity') + theme_minimal() + geom_point(size=3) +
  scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)]))

mds.fig <- ordiplot(ord, type = "none")

points(mds.fig, "sites", pch = 19, col = "navy", select = sample_info_tab1$season == 
         "winter")

points(mds.fig, "sites", pch = 19, col = "red", select = sample_info_tab1$season == 
         "summer")

points(mds.fig, "sites", pch = 19, col = "wheat4", select = sample_info_tab1$Source == 
         "Corrosion")

ordiellipse(ord, sample_info_tab1$Source, conf = 0.95, label = TRUE)




top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
physeq_top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
physeq_top20 <- prune_taxa(top20, physeq_top20)

Zixibacteria <- subset_taxa(physeq, Phylum %in% c('Zixibacteria'))
plot_tree(Zixibacteria, ladderize='left', size='abundance', color='Source', label.tips='Phylum') + scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)]))

Sulfurovaceae <- subset_taxa(physeq, Family %in% c('Sulfurovaceae'))
plot_tree(Sulfurovaceae, ladderize='left', size='abundance', color='Source', label.tips='Genus') + scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)]))


Shewanellaceae <- subset_taxa(physeq, Family %in% c('Shewanellaceae'))
plot_tree(Shewanellaceae, ladderize='left', size='abundance', color='Source', label.tips='Genus') + scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)]))


Desulfonatronum <- subset_taxa(physeq, Genus %in% c('Desulfonatronum', 'Desulfocarbo', 'Desulfomonile', 'Desulfomonile', 'Desulfobacula', 'Desulfobacca', 'Desulfocapsa', 'Desulfofrigus'))

plot_tree(Desulfatitalea, ladderize='left', size='abundance', color='Source', label.tips='Genus') + scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)]))


Candidatus <- subset_taxa(physeq, Genus %in% c('Candidatus Electrothrix', 'Candidatus Candidatus Magnetomorum'))
plot_tree(Candidatus , ladderize='left', size='abundance', color='char', label.tips='Genus') + scale_color_manual(values=unique(sample_info_tab1$Color[order(sample_info_tab1$Source)]))



Fe_oxid <- subset_taxa(physeq, Genus %in% c('Mariprofundus', 'Gallionella', 'Sideroxydans', 'Leptothrix', 'Thiobacillus', 'Dechloromonas'))
plot_tree(Fe_oxid, ladderize='left', size='abundance', color='Source', label.tips='Genus') + scale_color_manual(values=unique(sample_info_tab$Color[order(sample_info_tab$Source)]))


