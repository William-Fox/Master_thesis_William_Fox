library(tidyverse)
library(reshape2)
library(GenomicRanges)
library(readr)
library(UpSetR)
library(kableExtra)
library(dplyr)
library(GenomicFeatures)
library(stringr)
library(ballgown)
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(RColorBrewer)
library(gplots)
library(corrplot)
# awk -F'\t' -v OFS='\t' '{$10=$9; gsub(/.*transcript_id "/,"",$10); gsub(/".*/,"",$10); print}' CC_ALL_Copci2021.annotated.gtf > CC_ALL_Copci2021_ids
setwd("/Volumes/student_users/williamfox/R_final_final")
#load in gffcompare output with trnascripts_ids as a seperate column
Copci2021_gtf <- read.csv("CC_ALL_Cpci2021.annotated_known_all_structural_all_ids.gtf", sep = "\t", header = FALSE, quote = "")
colnames(Copci2021_gtf) <-  c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "Name" )

# load in coding predictions
hypothetical_prediction_df <- read.delim("hypothetical_Prediction", header =FALSE)
novel_prediction_df <- read.delim("new_novel_prediction", header = FALSE)
colnames(hypothetical_prediction_df) <- c("Program", "Prediction", "Sequence")
colnames(novel_prediction_df) <- c("Program", "Prediction", "Sequence")

# Reformat prediction
hypothetical_prediction_df_wide <- dcast(hypothetical_prediction_df, Sequence ~ Program, value.var = "Prediction")
novel_prediction_df_wide = dcast(novel_prediction_df, Sequence ~ Program, value.var = "Prediction")

# replace "na" from RNAcode column
hypothetical_prediction_df_wide$RNAcode <- hypothetical_prediction_df_wide$RNAcode %>% replace_na('None')
novel_prediction_df_wide$RNAcode <- novel_prediction_df_wide$RNAcode %>% replace_na('None')

## identify novel identified as noncoding
novel_noncoding_all <- novel_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "Noncoding")
novel_noncoding_none <- novel_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "None")
novel_ncRNA_prediction <- rbind(novel_noncoding_all, novel_noncoding_none)

##identify novel identified as coding
novel_coding <- novel_prediction_df_wide %>% filter(CPC2 == "Coding" | PLEK == "Coding" | RNAcode == "Coding")

#### identify hypothetical identified as noncoding
hypothetical_noncoding_all <- hypothetical_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "Noncoding")
hypothetical_noncoding_none <- hypothetical_prediction_df_wide %>% filter(CPC2 == "Noncoding" & PLEK == "Noncoding" & RNAcode == "None")
hypothetical_ncRNA_prediction <- rbind(hypothetical_noncoding_all, hypothetical_noncoding_none)

##identify hypothetical identified as coding
hypothetical_coding <- hypothetical_prediction_df_wide %>% filter(CPC2 == "Coding" | PLEK == "Coding" | RNAcode == "Coding")

#combining coding predictions
coding_putative <- rbind(hypothetical_coding, novel_coding)
noncoding_putative<- rbind(hypothetical_ncRNA_prediction, novel_ncRNA_prediction)
all_putative_analysis <- rbind(noncoding_putative,coding_putative)

##UpsetR plot of coding analysis
expressionInput <- c(RNAcode_none =nrow(subset(all_putative_analysis, RNAcode == "None"  & PLEK == "Coding" & CPC2 == "Coding")),
                     PLEK_noncoding = nrow(subset(all_putative_analysis, RNAcode == "Coding"  & PLEK == "Noncoding" & CPC2 == "Coding")), 
                     CPC2_noncoding = nrow(subset(all_putative_analysis, RNAcode == "Coding"  & PLEK == "Coding" & CPC2 == "Noncoding")),
                     RNAcode_noncoding = nrow(subset(all_putative_analysis, RNAcode == "Noncoding"  & PLEK == "Coding" & CPC2 == "Coding")), 
                     `PLEK_noncoding&RNAcode_none` = nrow(subset(all_putative_analysis, RNAcode == "None"  & PLEK == "Noncoding" & CPC2 == "Coding")), 
                     `PLEK_noncoding&RNAcode_noncoding` = nrow(subset(all_putative_analysis, RNAcode == "Noncoding"  & PLEK == "Noncoding" & CPC2 == "Coding")), 
                     `PLEK_noncoding&CPC2_noncoding` = nrow(subset(all_putative_analysis, CPC2 == "Noncoding"  & PLEK == "Noncoding" & RNAcode == "Coding")), 
                     `RNAcode_none&PLEK_noncoding&CPC2_noncoding` = nrow(subset(all_putative_analysis, RNAcode == "None"  & PLEK == "Noncoding" & CPC2 == "Noncoding")),
                     `CPC2_noncoding&RNAcode_none` = nrow(subset(all_putative_analysis, RNAcode == "None"  & PLEK == "Coding" & CPC2 == "Noncoding")), 
                     `CPC2_noncoding&RNAcode_noncoding` = nrow(subset(all_putative_analysis, RNAcode == "Noncoding"  & PLEK == "Coding" & CPC2 == "Noncoding")), 
                     `PLEK_noncoding&RNAcode_noncoding&CPC2_noncoding` = nrow(subset(all_putative_analysis, RNAcode == "Noncoding"  & PLEK == "Noncoding" & CPC2 == "Noncoding")))

pdf("Results/upsetR.pdf")
upset(fromExpression(expressionInput), order.by = "freq")   
dev.off()


#make into a format comparable to the gtf
noncoding_putative<- rbind(hypothetical_ncRNA_prediction, novel_ncRNA_prediction)
noncoding_putative <- separate(noncoding_putative, col = Sequence, into = c("seqname", "start", "end"), sep = "[:-]")
noncoding_putative$start <- as.numeric(noncoding_putative$start) + 1
noncoding_putative$end <- as.numeric(noncoding_putative$end)
noncoding_putative$seqname <- gsub(">","",noncoding_putative$seqname)
noncoding_putative[,7] <- "ncRNA"

#combine coding results to gtf
Copci2021_noncoding_gtf <- left_join(Copci2021_gtf, noncoding_putative, by= c("seqname", "start", "end"), multiple= "first")

# class codes into a seperate column
Copci2021_noncoding_gtf$class_code[Copci2021_noncoding_gtf$feature == "transcript"] <- gsub('.*class_code "([=a-z]).*','\\1', Copci2021_noncoding_gtf$attribute[Copci2021_noncoding_gtf$feature == "transcript"])

#load hypothetical interpro 
hypothetical_interpro <- read.csv("hypothetical_interpro.tabular", sep = "\t", header=FALSE)
colnames(hypothetical_interpro) [1] = "Name"
hypothetical_interpro <- hypothetical_interpro[grep("Pfam", hypothetical_interpro$V4),]

#add to gtf
Copci2021_noncoding_interpro <- left_join(Copci2021_noncoding_gtf, hypothetical_interpro, by="Name", multiple = "first")
Copci2021_noncoding_interpro <- Copci2021_noncoding_interpro[,c(1:10,14:15,18:20, 26:28)]



#add novel_interpro
novel_interpro <- read.csv("novel_interpro.tabular", sep = "\t", header=FALSE)
colnames(novel_interpro) [1] = "Name"
novel_interpro <- novel_interpro[grep("Pfam", novel_interpro$V5),]
novel_interpro <- separate(novel_interpro, col = Name, into = c("seqname", "start", "end"), sep = "[:-]")
novel_interpro$start <- as.numeric(novel_interpro$start)
novel_interpro$end <- as.numeric(novel_interpro$end)
Copci2021_noncoding_interpro <- left_join(Copci2021_noncoding_interpro, novel_interpro, by= c("seqname", "start", "end"), multiple = "first")

#table of Go terms for each transcript
go_terms <- Copci2021_noncoding_interpro[, c("Name", "V14.x", "V15")]
go_terms$V14.x <- replace(go_terms$V14.x, go_terms$V14.x == "-", NA)
go_terms$V15 <- replace(go_terms$V15, go_terms$V15 == "-", NA)
go_terms <- go_terms %>%
  filter(!is.na(V14.x) | !is.na(V15))
go_terms <- go_terms %>%
  mutate(V14.x = coalesce(V14.x, V15),
         V15 = ifelse(is.na(V14.x), V15, NA))
go_terms <- go_terms[, c("Name", "V14.x")]


#clean up table
Copci2021_noncoding_interpro <- Copci2021_noncoding_interpro[,c(1:18,22:24,30:32)]

# replace ncRNA if a feature has pfam results
Copci2021_noncoding_interpro$V7.x[Copci2021_noncoding_interpro$V4.x == "Pfam" | Copci2021_noncoding_interpro$V5.y == "Pfam"] <- "coding"

# make a table of just transcripts and key info for investigation
Copci2021_annotated_transcripts <- Copci2021_noncoding_interpro[grep("transcript", Copci2021_noncoding_interpro$feature), ]
Copci2021_annotated_transcripts$Name <- ifelse(
  is.na(Copci2021_annotated_transcripts$Name) | Copci2021_annotated_transcripts$Name == "",
  str_extract(Copci2021_annotated_transcripts$attribute, "(?<=transcript_id=)[^;]+"),
  Copci2021_annotated_transcripts$Name
)
Copci2021_annotated_transcripts$Feature_length <- ( Copci2021_annotated_transcripts$end - Copci2021_annotated_transcripts$start ) +1
Copci2021_annotated_transcripts$V7.x[grepl("unassigned_transcript", Copci2021_annotated_transcripts$attribute)] <- "ncRNA"
Copci2021_annotated_transcripts$class_code[grepl("unassigned_transcript", Copci2021_annotated_transcripts$attribute)] <- "tRNA"


Key_features <- Copci2021_annotated_transcripts[,c(1,3:5,10:12)]
Key_features$V7.x <- ifelse(is.na(Key_features$V7.x), "coding", Key_features$V7.x)
Key_features <- Key_features[grepl("[ixump=t]", Key_features$class_code), ]

Key_features$Feature_length <- ( Key_features$end - Key_features$start ) +1
Key_features_coding <- Key_features[grepl("coding", Key_features$V7.x), ]
Key_features_coding <- Key_features_coding[grepl("=", Key_features_coding$class_code), ]
Key_features_coding <- unique(Key_features_coding)
Key_features_ncRNA <- Key_features[grepl("ncRNA", Key_features$V7.x), ]
Key_features_ncRNA <- unique(Key_features_ncRNA)

Key_features_coding_all <- Key_features[grepl("coding", Key_features$V7.x), ]

ncRNA_lengths <- Key_features$Feature_length[Key_features$V7.x == "ncRNA"]
coding_lengths <- Key_features$Feature_length[Key_features$V7.x == "coding" & Key_features$class_code == "="]
feature_lengths <- data.frame(length=c(ncRNA_lengths, coding_lengths), type=c(rep("ncRNA", length(ncRNA_lengths)), rep("mRNA", length(coding_lengths))))
fake_coding_lengths <- data.frame(Key_features$Feature_length[Key_features$class_code == "="])
colnames(fake_coding_lengths) <- "length"
fake_coding_lengths$length <- as.numeric(fake_coding_lengths$length)
print(summary(fake_coding_lengths))


# Create histogram using ggplot2, only = coding and all ncRNA
ggplot(feature_lengths, aes(x=length, fill=type)) +
  geom_histogram(binwidth=100) +
  xlab("Feature length") +
  ylab("Count") +
  ggtitle("Histogram of feature lengths by annotation type") +
  scale_x_continuous(labels = function(x) x/1000, name = "Feature length (Kb)", minor_breaks = seq(0, 8000, by = 100), limits=c(0, 8000))

Key_features$V7.x <- sub("coding", "mRNA", Key_features$V7.x)
# Create histogram using ggplot2, all coding and all ncRNA
ggplot(Key_features, aes(x=Feature_length, fill = V7.x)) +
  geom_histogram(binwidth = 100) +
  scale_x_continuous(labels = function(x) x/1000, name = "Feature length (Kb)", minor_breaks = seq(0, 8000, by = 100), limits=c(0, 8000)) +
  labs(fill = "RNA type") +
  theme(panel.grid.major = element_line(color= "white", size=1), 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14)) +
  scale_y_continuous( breaks = seq(0,1500,100))
ggsave("ncRNA_coding_hist.pdf")

ggplot(Key_features, aes(y=Feature_length, fill = V7.x, x= V7.x)) +
  geom_boxplot() +
  scale_y_continuous(name = "Feature length (Kb)", labels = function(x) x/1000, 
                     breaks = seq(0,20000,1000)) +
  theme(panel.grid.major = element_line(color= "white", size=1), 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14)) +
  theme(legend.position = "none") +
  labs(fill = "RNA type", x = "RNA type") 
ggsave("ncRNA_coding_boxplot.pdf")

#mann whit stat test of lengths
coding_ncRNA_mann_whit <-wilcox.test(ncRNA_lengths, coding_lengths)
print(coding_ncRNA_mann_whit)

ncRNA_coding_summary <- Key_features %>% 
  group_by(V7.x) %>% 
  summarise(
    Feature_length_median = median(Feature_length),
    Feature_length_min = min(Feature_length),
    Feature_length_max = max(Feature_length),
    Feature_length_q1 = quantile(Feature_length, 0.25),
    Feature_length_q3 = quantile(Feature_length, 0.75),
    Feature_length_iqr = IQR(Feature_length)
  )

##legnth of class codes
class_lengths <- Key_features[, c("class_code", "Feature_length")]
class_lengths <- class_lengths[!grepl("transcript_id", class_lengths$class_code), ]
class_lengths <- class_lengths[!grepl("=", class_lengths$class_code), ]
class_lengths <- class_lengths[!grepl("tRNA", class_lengths$class_code), ]
mean_lengths <- aggregate(Feature_length ~ class_code, data = class_lengths, FUN = summary)
write.csv(mean_lengths, "Results/class_length_summary_stats.csv")

ggplot(class_lengths, aes(x =Feature_length, fill = class_code)) +
  geom_histogram(binwidth = 100)

ggsave("Results/class_code_feature_length.pdf")

## count key features
count_var_protein <-  Key_features_coding$class_code == '='
Key_features_ncRNA$class_code <- ifelse(grepl("transcript_id", Key_features_ncRNA$class_code), "structural only", Key_features_ncRNA$class_code)


Key_features_df <- data.frame(table(Key_features_ncRNA$class_code))
Key_features_df <- rbind(Key_features_df, data.frame(table(Key_features_coding$class_code)))
Key_features_df$Var1  <- c("class code = ncRNA", "class code i ncRNA", "class code m ncRNA", "class code p ncRNA", "structural only ncRNA", "tRNA", "class code u ncRNA", "class code x ncRNA", "Protein coding genes")
write.csv(Key_features_df, "Results/novel_ncRNA_count.csv")
#make into a table

key_feature_coding_names <- Key_features_coding$Name
key_feature_ncRNA_names <- Key_features_ncRNA$Name

# make table into gtf style 
Copci2021_noncoding_interpro$combined_attributes <- apply(Copci2021_noncoding_interpro[, c(9, 11:24)], 1, function(x) {
  paste(ifelse(is.na(x) | x == "", "", x), collapse="; ")
})

Copci2021_noncoding_interpro  <- subset(Copci2021_noncoding_interpro, select=c(1:8, 25))
write.table(Copci2021_noncoding_interpro, "Results/Copci2021_final.gtf", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#DESeq2 analysis of feature counts data using newly constructed
##lets try DEseq2
library(viridis)
library(DESeq2)
library(ggdendro)
library(gridExtra)
library(grid)
library(gtable)
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
gene2go <- readMappings("gene_go.tbl")
setwd("/Volumes/student_users/williamfox/R_final_final/")
CC_data= read.csv('CC_IDs.csv')

fcData <- read.table("counts.txt", sep='\t', header=TRUE)
fcData <- fcData[,c(1, 7:35)]
metaData <- CC_data
colnames(metaData) <- c("ID", "stage", "loc", "sample")
metaData$ID <- gsub("^CC_|\\.sam$", "", metaData$ID)
metaData$sample <- gsub("^CC_|\\.sam$", "", metaData$sample)

rownames(fcData) <- fcData[, 1]

# Remove gene_ids column from fcData as it's now redundant
fcData <- fcData[, -1]

# Make sure the column names of fcData (samples) match the IDs of metaData
all(colnames(fcData) == metaData$ID) # should return TRUE

# Set the rownames of metaData to match the IDs
rownames(metaData) <- metaData$ID

# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = fcData, 
                              colData = metaData, 
                              design = ~ sample)

# Run DESeq
dds <- DESeq(dds)

# Create a list to store results for each sample
results_list <- list()

for (sample in unique(metaData$sample)) {
  # Set the reference level for this comparison
  dds$sample <- relevel(dds$sample, ref = sample)
  
  # Re-run DESeq
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  
  # Remove rows with NA in padj or log2FoldChange
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  
  # Filter for significant genes
  sig_res_up <- res[res$padj < 0.05 & res$log2FoldChange > 4, ]
  sig_res_down <- res[res$padj < 0.05 & res$log2FoldChange < -4, ]
  
  # Store in the list
  results_list[[sample]] <- list("Up" = sig_res_up, "Down" = sig_res_down)
  combined <- rbind(cbind(sig_res_up, "Direction" = "Up"), cbind(sig_res_down, "Direction" = "Down"))
  
  # Write to a CSV file
  write.csv(combined, file.path("Results/new", paste0("deseq_results_", sample, ".csv")), row.names = TRUE)
}






## calculate vst and reformat data frame
vst_transformer <- vst(dds, blind = FALSE)
vst_data <- assay(vst_transformer)
vst_data <- as.data.frame(vst_data)
colnames(vst_data) <- metaData$ID
vst_data_sub <- vst_data

long_vst_data_sub <- vst_data_sub %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(cols = -gene_id, names_to = "ID", values_to = "expression")

long_vst_data_sub <- left_join(long_vst_data_sub, metaData, by = "ID")

avg_vst_data_sub <- long_vst_data_sub %>%
  group_by(sample, gene_id) %>%
  summarise(expression = mean(expression, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = expression)

rownames(avg_vst_data_sub) <- avg_vst_data_sub$gene_id

## coding
avg_vst_data_sub_coding <-avg_vst_data_sub[avg_vst_data_sub$gene_id %in% key_feature_coding_names, ]

avg_vst_data_sub_coding$gene_id <- NULL
avg_vst_data_sub_coding <- avg_vst_data_sub_coding[, c("VM", "H", "P1", "P2", "YFB_S", "YFB_L", "YFB_C", "FB_S" ,"FB_CL")]

#make and save heatmap
pdf("Results/New_heatmap_coding.pdf")
heatmap.2(as.matrix(avg_vst_data_sub_coding), trace = "none", 
          col = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(255),
          scale = "row", 
          dendrogram = "none",
          labRow = FALSE,
          hclustfun = function(x) hclust(x, method = "complete"), 
          distfun = function(x) dist(x, method = "euclidean"))
dev.off()

avg_vst_data_sub_nc <-avg_vst_data_sub[avg_vst_data_sub$gene_id %in% key_feature_ncRNA_names, ]

avg_vst_data_sub_nc$gene_id <- NULL
avg_vst_data_sub_nc <- avg_vst_data_sub_nc[, c("VM", "H", "P1", "P2", "YFB_S", "YFB_L", "YFB_C", "FB_S" ,"FB_CL")]

#make and save heatmap
pdf("Results/New_heatmap_nc.pdf")
heatmap.2(as.matrix(avg_vst_data_sub_nc), trace = "none", 
          col = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(255),
          dendrogram = "none",
          scale = "row", 
          labRow = FALSE,
          hclustfun = function(x) hclust(x, method = "complete"),
          distfun = function(x) dist(x, method = "euclidean"))
dev.off()


##get the lists
upregulated_VM <- rownames(results_list[["VM"]][["Up"]])
downregulated_VM <- rownames(results_list[["VM"]][["Down"]])

# Get all gene names
all_genes <- rownames(fcData)

# Initialize gene list with all genes, marking them as not in the set of genes of interest
gene_list <- rep(0, length(all_genes))
names(gene_list) <- all_genes

# Loop over the results list
for (sample in names(results_list)) {
  # Get upregulated and downregulated gene names
  upregulated_genes <- rownames(results_list[[sample]][["Up"]])
  downregulated_genes <- rownames(results_list[[sample]][["Down"]])
  
  # Mark these genes as being in the set of genes of interest
  gene_list[upregulated_genes] <- 1
  gene_list[downregulated_genes] <- 1
}

# Specify root nodes
root_nodes <- c("MF", "BP", "CC")

# Function to select genes of interest
sel_funct <- function(all_genes) {
  return(all_genes == 1)
}

# Function for GO enrichment analysis
go_enrichment <- function(gene_list, sample_name, up_down) {
  for (ontology in root_nodes) {
    # Prepare topGO object
    topgo_data <- new("topGOdata",
                      ontology = ontology,
                      allGenes = gene_list,
                      geneSel = sel_funct,
                      annotationFun = annFUN.gene2GO,
                      gene2GO = gene2go)
    
    # Run topGO analysis
    resultFisher <- runTest(topgo_data, algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(topgo_data, algorithm="classic", statistic="ks")
    resultKS.elim <- runTest(topgo_data, algorithm ="elim",statistic ="ks")
    
    # Run Fisher's weight01 test
    weight01.fisher <- runTest(topgo_data, statistic = "fisher")
    
    # Create a table of the results
    allRes <- GenTable(topgo_data, classicFisher = resultFisher,
                       classicKS = resultKS, elimKS = resultKS.elim,
                       weight01Fisher = weight01.fisher,
                       orderBy = "weight01Fisher",
                       ranksOf = "classicFisher", topNodes = 50)
    
    
    # Write results to file
    output_filename <- paste("/Volumes/student_users/williamfox/R_final_final/go_enrichment_samples/go_enrichment_", sample_name, "_", up_down, "_", ontology, ".csv", sep = "")
    write.csv(allRes, output_filename)
  }
}

# Perform GO enrichment for each sample's list of up and downregulated genes
for (sample in names(results_list)) {
  upregulated_genes <- rownames(results_list[[sample]][["Up"]])
  downregulated_genes <- rownames(results_list[[sample]][["Down"]])
  
  # Create gene list for upregulated genes
  gene_list <- rep(0, length(all_genes))
  names(gene_list) <- all_genes
  gene_list[upregulated_genes] <- 1
  
  # Perform GO enrichment analysis
  go_enrichment(gene_list, sample, "Up")
  
  # Create gene list for downregulated genes
  gene_list <- rep(0, length(all_genes))
  names(gene_list) <- all_genes
  gene_list[downregulated_genes] <- 1
  
  # Perform GO enrichment analysis
  go_enrichment(gene_list, sample, "Down")
}

## Graph of regulated genes
library(dplyr)
library(purrr)

# Function to process each DESeqResults object
process_DESeqResults <- function(item, name, group) {
  as.data.frame(item) %>% 
    rownames_to_column("Row") %>% 
    mutate(Name = name, Group = group)
}

# Initialize empty dataframes
up_df <- data.frame()
down_df <- data.frame()

# Loop over each element in the up_list and down_list
for (name in names(results_list)) {
  up_df <- bind_rows(up_df, process_DESeqResults(results_list[[name]]$Up, name, "up"))
  down_df <- bind_rows(down_df, process_DESeqResults(results_list[[name]]$Down, name, "down"))
}

# Combine dataframes
combined_df <- bind_rows(up_df, down_df)
combined_df$Row <- gsub("ALL", "CC_ALL", combined_df$Row)
combined_df <- combined_df[,c(1, 8:9)]
colnames(combined_df) <- c("transcript", "sample", "direction")
combined_df <- left_join(combined_df, Key_features, by = c("transcript" = "Name"), relationship = "many-to-many")
combined_df <- drop_na(combined_df, V7.x)
diff_expr_ncRNA <- combined_df %>% filter(V7.x == "ncRNA")
table(diff_expr_ncRNA$class_code)
combined_df <- combined_df[,c(1:3, 8)]
combined_df$sample <- factor(combined_df$sample, levels = c("VM", "H", "P1", "P2", "YFB_S", "YFB_L", "YFB_C", "FB_S" ,"FB_CL"))
combined_df$Type_Direction <- interaction(combined_df$V7.x, combined_df$direction, sep = ".", lex.order = TRUE)
combined_df$Type_Direction <- factor(combined_df$Type_Direction, levels = c("ncRNA.down", "ncRNA.up", "mRNA.down", "mRNA.up"))
table_df <- as.data.frame.table(table(combined_df$sample, combined_df$Type_Direction), responseName = "count")
table_df <- dcast(table_df, Var1 ~ Var2)
write.csv(table_df, "diff_expr_count.csv", quote = FALSE)


##make a plot
ggplot(combined_df, aes(x = sample, fill = interaction(V7.x, direction))) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(x = "Sample", y = "Count", fill = "V7.x-Direction Combination") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(combined_df, aes(x = sample, fill = direction)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(x = "Sample", y = "Proportion", fill = "Direction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(combined_df, aes(x = sample, fill = Type_Direction)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("mRNA.down" = "#078188", "mRNA.up" = "#16E8F3", "ncRNA.up" = "#FC746C", "ncRNA.down" = "#B40C04"),
                    labels = c("Down-regulated ncRNA", "Up-regulated ncRNA", "Down-regulated mRNA", "Up-regulated mRNA")) +
  labs(x = "Growth stage and location", y = "Proportion", fill = "transcript type and \nregulation direction") +
  theme(panel.grid.major = element_line(color= "white", size=1), 
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Results/transcript_regulation.pdf")

## correlation of mRNA - ncRNA pairs mentioned in this study
## make the list
corr_list <- read.table("correlation-list.txt")
corr_list <-  df_agg_sample[df_agg_sample$transcript_id %in% corr_list$V1, ]
selected_cols <- c("transcript_id", "FB_CL", "FB_S", "H", "P1", "P2", "VM", "YFB_C", "YFB_L", "YFB_S")
corr_list <- corr_list[, selected_cols]

#measure ccf (cross-correlation)
# List to store results
ccf_results <- list()

# Get the gene names from corr_list
genes <- corr_list$transcript_id

# Loop over each pair of genes
for(i in 1:(length(genes)-1)) {
  for(j in (i+1):length(genes)) {
    # Extract the time series for the two genes
    gene1 <- as.numeric(corr_list[i, -1]) # -1 to exclude the transcript_id column
    gene2 <- as.numeric(corr_list[j, -1]) # -1 to exclude the transcript_id column
    
    # Calculate the cross-correlation
    ccf_res <- ccf(gene1, gene2, plot = FALSE) # plot=FALSE to not generate plots
    
    # Store the result
    ccf_results[[paste(genes[i], genes[j], sep="___")]] <- ccf_res
  }
}

# Initialize an empty data frame to store results
ccf_df <- data.frame(
  gene1 = character(),
  gene2 = character(),
  max_acf = numeric(),
  lag_at_max_acf = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each entry in ccf_results
for (gene_pair in names(ccf_results)) {
  # Extract gene names by splitting at "___"
  genes <- strsplit(gene_pair, "___")[[1]]
  
  # Extract ccf results
  ccf_res <- ccf_results[[gene_pair]]
  
  # Find the maximum absolute correlation and its corresponding lag
  max_acf <- max(abs(ccf_res$acf))
  lag_at_max_acf <- ccf_res$lag[which.max(abs(ccf_res$acf))]
  
  # Create a data frame for this pair
  pair_df <- data.frame(
    gene1 = genes[1],
    gene2 = genes[2],
    max_acf = max_acf,
    lag_at_max_acf = lag_at_max_acf,
    stringsAsFactors = FALSE
  )
  
  # Append to the results data frame
  ccf_df <- rbind(ccf_df, pair_df)
}

correlation_matrix <- cor(t(corr_list[,-1]), use = "pairwise.complete.obs")
rownames(correlation_matrix) <- colnames(correlation_matrix) <- corr_list$transcript_id


##Zymo short selected coding prediction
zymo_prediction <- read.delim("Zymo_prediction.txt")
zymo_gtf <- read.csv("Zymo_class_codes.gtf", sep = "\t", header = FALSE, quote = "")
colnames(zymo_gtf) <-  c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


colnames(zymo_prediction) <- c("Program", "Prediction", "Sequence")

# Reformat prediction
zymo_prediction_wide <- dcast(zymo_prediction, Sequence ~ Program, value.var = "Prediction")
zymo_prediction_wide$CPC2 <- ifelse(is.na(zymo_prediction_wide$CPC2), "none", zymo_prediction_wide$CPC2)
zymo_prediction_wide$PLEK <- ifelse(is.na(zymo_prediction_wide$PLEK), "none", zymo_prediction_wide$PLEK)



##identify novel identified as coding
zymo_coding <- zymo_prediction_wide %>% filter(CPC2 == "Coding" | PLEK == "Coding")
zymo_coding <- separate(zymo_coding, col = Sequence, into = c("seqname", "start", "end"), sep = "[:-]")
zymo_coding$start <- as.numeric(zymo_coding$start) + 1
zymo_coding$end <- as.numeric(zymo_coding$end)
zymo_coding$seqname <- gsub(">","",zymo_coding$seqname)
zymo_coding[,6] <- "coding"

#combine coding results to gtf

zymo_gtf <- left_join(zymo_gtf, zymo_coding, by= c("seqname", "start", "end"), multiple= "first")
zymo_gtf$V6 <- ifelse(is.na(zymo_gtf$V6), "ncRNA", zymo_gtf$V6)
# class codes into a seperate column and save gtf
zymo_gtf$classcode <- str_extract(zymo_gtf$attribute, "(?<=class_code )[A-Za-z]")
write.table(zymo_gtf, "Results/Zymo_final.gtf", sep= '\t', quote = FALSE)

#look at class codes by RNA type
zymo_nc <- zymo_gtf[zymo_gtf$feature == "transcript" & zymo_gtf$V6 == "ncRNA", ]
table(zymo_nc$classcode)
zymo_c <- zymo_gtf[zymo_gtf$feature == "transcript" & zymo_gtf$V6 == "coding", ]
table(zymo_c$classcode)


