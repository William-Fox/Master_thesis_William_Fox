#packages
library(ggplot2)

##load gtf
Copci2021_gtf <- read.csv("Copci2021_ids.gtf", sep = "\t", header = TRUE, quote = "")
colnames(Copci2021_gtf) <-  c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "Name" )

#load intron and intergenic gtf
intron_gtf <-  read.delim("input_genome.final.intron.gff", header = FALSE)
intergenic_gtf <-  read.delim("input_genome.final.ig.gff", header = FALSE)


##hypothetical proteins
Copci2021_hypotheticals <- Copci2021_gtf[grep("hypothetical", Copci2021_gtf$attribute),] 
Copci2021_hypotheticals <- Copci2021_hypotheticals[grep("transcript", Copci2021_hypotheticals$feature),]
Copci2021_hypotheticals_length <- length(Copci2021_hypotheticals$seqname)


##tRNA
Copci2021_trna <- Copci2021_gtf[grep('transcript_biotype "tRNA"', Copci2021_gtf$attribute),] 
Copci2021_trna <- Copci2021_trna[grep("transcript", Copci2021_trna$feature),]
Copci2021_trna_length <- length(Copci2021_trna$seqname)

unassigned_transcripts <- Copci2021_gtf[grep("unassigned_transcript", Copci2021_gtf$attribute),]
unassigned_transcripts <- unassigned_transcripts[grep("transcript", unassigned_transcripts$feature),]
## Known proteins
Copci2021_known_proteins <- Copci2021_gtf[grep("hypothetical", Copci2021_gtf$attribute, invert = TRUE),]
Copci2021_known_proteins <- Copci2021_known_proteins[grep("transcript_biotype tRNA", Copci2021_known_proteins$attribute, invert = TRUE),]
Copci2021_known_proteins <- Copci2021_known_proteins[grep("unassigned_transcript", Copci2021_known_proteins$attribute, invert = TRUE),]
Copci2021_known_proteins <- Copci2021_known_proteins[grep("transcript", Copci2021_known_proteins$feature),]
Copci2021_known_proteins$length <- as.numeric(Copci2021_known_proteins$end[Copci2021_known_proteins$feature == "transcript"] - Copci2021_known_proteins$start[Copci2021_known_proteins$feature == "transcript"]) +1
Copci2021_known_proteins_length <- length(Copci2021_known_proteins$seqname)

Copci2021_known_short <- Copci2021_known_proteins[Copci2021_known_proteins$length <= 200, ]



## calculate length of features
exon_length <- abs(as.numeric(Copci2021_gtf$end[Copci2021_gtf$feature == "exon"] - Copci2021_gtf$start[Copci2021_gtf$feature == "exon"]) + 1)

known_length <- abs(as.numeric(Copci2021_known_proteins$end[Copci2021_known_proteins$feature == "transcript"] - Copci2021_known_proteins$start[Copci2021_known_proteins$feature == "transcript"]) +1)

hypothetical_length <- abs(as.numeric(Copci2021_hypotheticals$end[Copci2021_hypotheticals$feature == "transcript"] - Copci2021_hypotheticals$start[Copci2021_hypotheticals$feature == "transcript"]) +1)

trna_length <- abs(as.numeric(Copci2021_trna$end[Copci2021_trna$feature == "transcript"] - Copci2021_trna$start[Copci2021_trna$feature == "transcript"]) +1)

intron_length <- abs(as.numeric(intron_gtf$V5 - intron_gtf$V4) +1)

ig_length <- abs(as.numeric(intergenic_gtf$V5 - intergenic_gtf$V4) +1)



##data frame of lengths
Copci2021_lengths <- data.frame(lengths=c(exon_length, known_length, hypothetical_length, trna_length, intron_length, ig_length),
                                Feature=c(rep("Exon", length(exon_length)), rep("Characterised protein mRNA", length(known_length)),
                                          rep("Uncharacterised protein mRNA", length(hypothetical_length)), rep("tRNA", length(trna_length)), 
                                          rep("Intron", length(intron_length)), rep("Intergenic", length(ig_length))))

Copci2021_proteins_lengths <- data.frame(lengths=c(known_length, hypothetical_length),
                                         Feature=c(rep("Characterised protein mRNA", length(known_length)),rep("Uncharacterised protein mRNA", length(hypothetical_length))))
hypothetical_summary <-summary(hypothetical_length)
known_summary <- summary(known_length)
protein_length_summary <- cbind(hypothetical_summary,known_summary)



## remove outliers
Copci2021_lengths_clean <- Copci2021_lengths[Copci2021_lengths$lengths <= 8000, ]

##histogram of size frequency of features in original annotation
ggplot(Copci2021_lengths_clean, aes(x = lengths, fill = Feature)) +
  geom_histogram(binwidth = 100) +
  labs(x = "Length (bp)", y = "Count") +
  facet_wrap(~factor(Feature, levels = c('Characterised protein mRNA', 'Exon', 'Uncharacterised protein mRNA', 
                                         'Intron',  'tRNA','Intergenic')),
             ncol = 2, scales = "free_y", strip.position = "top") +
  scale_fill_brewer(palette = "Set2") +
  theme(panel.grid.major = element_line(color= "white", size=1), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position= "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(labels = function(x) x/1000, name = "Feature length (Kb)", minor_breaks = seq(0, 8000, by = 100))

ggsave("Results/original_feature_length_hist.pdf")  
