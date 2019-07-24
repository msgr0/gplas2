#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 

suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(Biostrings))

path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["links"]]
classifier <- snakemake@params[["classifier"]]
threshold <- snakemake@params[["threshold"]]

raw_nodes <- readDNAStringSet(filepath = path_nodes, format="fasta")

raw_contig_names <- names(raw_nodes)

raw_number <- str_split_fixed(string = raw_contig_names, pattern = '_', n = 2)[,1]
number <- gsub(pattern = 'S', replacement = '', x = raw_number)

raw_length <- str_split_fixed(string = raw_contig_names, pattern = ':', n = 4)[,3]
length <- gsub(pattern = '_dp', replacement = '', x = raw_length)

coverage <- str_split_fixed(string = raw_contig_names, pattern = ':', n = 5)[,5]

contig_info <- data.frame(number = number,
                          length = length,
                          coverage = coverage)


contig_info$length <- as.numeric(as.character(contig_info$length)) # Converting the column length into a numeric column
contig_info$coverage <- as.numeric(as.character(contig_info$coverage)) # Converting the column coverage into a coverage column

graph_pos_contigs <- contig_info
graph_pos_contigs$number <- paste(graph_pos_contigs$number, '+', sep = '')

graph_neg_contigs <- contig_info
graph_neg_contigs$number <- paste(graph_neg_contigs$number, '-', sep = '')

graph_contigs <- rbind(graph_pos_contigs, graph_neg_contigs)

write.table(x = graph_contigs, file = snakemake@output[["graph_contigs"]], row.names = FALSE)

small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

raw_links <- read.table(path_links)
colnames(raw_links) <- c('L','first_node','first_sign','second_node','second_sign','OM')

links <- NULL
for(row in 1:nrow(raw_links))
{
  connection <- raw_links[row,]

  reverse_first_sign <- ifelse(connection$first_sign == '+','-','+')
  reverse_second_sign <- ifelse(connection$second_sign == '+','-','+')

  info_reverse <- c('L',connection$second_node,reverse_second_sign,connection$first_node,reverse_first_sign,'OM')

  clean_first <- paste(connection$first_node, connection$first_sign, sep = '')
  clean_second <- paste(connection$second_node, connection$second_sign, sep = '')

  info_forward <- c(clean_first,'to',clean_second)

  clean_rev_first <- paste(connection$second_node, reverse_second_sign, sep = '')
  clean_rev_second <- paste(connection$first_node, reverse_first_sign, sep = '')

  info_reverse <- c(clean_rev_first,'to',clean_rev_second)

  full_connection <- rbind(info_forward, info_reverse)
  links <- rbind(links, full_connection)

}

links <- as.data.frame(links)

write.table(x = links, file = snakemake@output[["clean_links"]], row.names = FALSE)

unique_nodes <- unique(links$V1)

# Loop to create a dataframe with the number of links from each node

repeat_info <- NULL
for(node in unique_nodes)
{
  repeat_links <- subset(links, links$V1 == node)
  repeat_links <- repeat_links[! duplicated(repeat_links),]

  node_info <- data.frame(number = node,
                          connecting_nodes = length(repeat_links$V3))

  repeat_info <- rbind(repeat_info, node_info)
}

repeats <- subset(repeat_info, repeat_info$connecting_nodes > 1) # Transposases may be identified as hubs in the graph. They should have more than one link

repeats_graph <- repeats

write.table(x = repeats_graph, file = snakemake@output[["graph_repeats"]])




