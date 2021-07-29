#!/usr/bin/Rscript

args <- commandArgs(TRUE)

if(!"tidyverse" %in% rownames(installed.packages())) {
  print("Installing tidyverse; please be patient")
  ###devtools::install_git("https://gitlab.com/mmb-umcu/mlplasmids.git", ref = "efaecalis")
  install.packages("tidyverse")
}


suppressMessages(library(tidyverse))

input_path <- args[1]
output_path <- args[2]

kraken_prediction <- read.csv(file = input_path, header = FALSE, sep = '\t')
colnames(kraken_prediction) <- c('C','spades_contig','Taxonomy','Length','Kmer_information')

kraken_prediction$Prediction <- ifelse(kraken_prediction$Taxonomy == '1351', 'Chromosome',
                                     ifelse(kraken_prediction$Taxonomy == '36549', 'Plasmid',
                                            ifelse(kraken_prediction$Taxonomy == '1', 'Unclassified', 'Other')))

kraken_prediction$Frequency_class <- 0

for(row in 1:nrow(kraken_prediction))
{
  df <- data.frame(Kmer_information = matrix(as.list(unlist(strsplit(as.character(kraken_prediction$Kmer_information[row]), '[[:space:]]')))))
  
  df$taxonomy <- str_split_fixed(string = df$Kmer_information, pattern = ':', n = 2)[,1]  
  df$bp <- as.numeric(str_split_fixed(string = df$Kmer_information, pattern = ':', n = 2)[,2])
  
  taxonomy_bp <- df %>%
    group_by(taxonomy) %>%
    summarise(sum(bp))
  
  taxonomy_bp$Frequency <- taxonomy_bp$`sum(bp)`/sum(taxonomy_bp$`sum(bp)`)
  
  taxonomy_bp <- taxonomy_bp[order(taxonomy_bp$Frequency, decreasing = TRUE),]
  
  kraken_prediction$Frequency_class[row] <- taxonomy_bp$Frequency[1]
  
}

ambiguous_contigs <- which(kraken_prediction$Prediction %in% c('Unclassified','Other'))

kraken_prediction$Prediction[ambiguous_contigs] <- 'Plasmid'
kraken_prediction$Frequency_class[ambiguous_contigs] <- 0.5

kraken_prediction$Frequency_other_class <- 1-kraken_prediction$Frequency_class

kraken_chr <- subset(kraken_prediction, kraken_prediction$Prediction == 'Chromosome')
kraken_pl <- subset(kraken_prediction, kraken_prediction$Prediction == 'Plasmid')


mlplasmids_chr_df <- data.frame(Prob_Chromosome = kraken_chr$Frequency_class, 
                                Prob_Plasmid = kraken_chr$Frequency_other_class,
                                Prediction = kraken_chr$Prediction,
                                Contig_name = kraken_chr$spades_contig,
                                Contig_length = kraken_chr$Length)

mlplasmids_pl_df <- data.frame(Prob_Chromosome = kraken_pl$Frequency_other_class, 
                                Prob_Plasmid = kraken_pl$Frequency_class,
                                Prediction = kraken_pl$Prediction,
                                Contig_name = kraken_pl$spades_contig,
                                Contig_length = kraken_pl$Length)

mlplasmids_df <- rbind(mlplasmids_chr_df, mlplasmids_pl_df)

mlplasmids_df <- mlplasmids_df[order(mlplasmids_df$Contig_length, decreasing = TRUE),]

write.table(x = mlplasmids_df, file = output_path, row.names = FALSE, sep = '\t')

