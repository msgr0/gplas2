#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(cooccur))
suppressMessages(library(ggrepel))
suppressMessages(library(Biostrings))
suppressMessages(library(seqinr))
suppressMessages(library(MCL))

# Inputs

path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["clean_links"]]
path_prediction <- snakemake@input[["clean_prediction"]]
path_graph_contigs <- snakemake@input[["graph_contigs"]]
path_graph_repeats <- snakemake@input[["graph_repeats"]]
path_init_nodes <- snakemake@input[["initialize_nodes"]]
path_cov_variation <- snakemake@input[["coverage"]]
input_solutions <- snakemake@input[["solutions"]]
classifier <- snakemake@params[["classifier"]]
threshold <- snakemake@params[["threshold"]]
iterations <- snakemake@params[["iterations"]]

links <- read.table(file = path_links, header = TRUE)
graph_contigs <- read.table(file = path_graph_contigs, header = TRUE)

small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

repeats_graph <- read.table(file = path_graph_repeats)
repeats <- repeats_graph
repeats$number <- gsub(pattern = '\\+',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)
repeats$number <- gsub(pattern = '\\-',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)

initialize_nodes <- read.table(file = path_init_nodes, header = TRUE)
initialize_nodes <- initialize_nodes[,1]

clean_pred <- read.table(file = path_prediction, header = TRUE)


max_variation <- read.table(file = path_cov_variation, header = TRUE)
max_variation <- as.numeric(max_variation[1,1])
max_variation_small <- max_variation*5.0

max_nodes <- max(count.fields(input_solutions, sep = ","))

solutions <- read.table(input_solutions ,sep=",",fill=TRUE,col.names=1:max_nodes)


all_nodes <- NULL
  for(solution in 1:nrow(solutions))
  {
    iteration <- solutions[solution,]
    iteration <- t(iteration)
    nodes <- as.character(iteration[,1])
    nodes <- nodes[nodes != '']
    all_nodes <- append(x = all_nodes, values = nodes, after = length(all_nodes))
  }

unique_nodes <- unique(all_nodes)
unique_nodes <- unique_nodes[unique_nodes != ""]
unique_nodes <- as.character(na.omit(unique_nodes))

co_ocurrence <- data.frame(matrix(0, nrow = nrow(solutions), ncol = length(unique_nodes)))
colnames(co_ocurrence) <- unique_nodes

for(row in 1:nrow(co_ocurrence))
{
  iteration <- co_ocurrence[row,]
  particular_solution <- solutions[row,]
  particular_solution <- t(particular_solution)
  particular_solution <- particular_solution[,1]
  particular_solution <- particular_solution[particular_solution != ""]
  presence_absence <- ifelse(colnames(co_ocurrence) %in% particular_solution == TRUE, 1, 0)
  co_ocurrence[row,] <- presence_absence
}

co_ocurrence <- apply(co_ocurrence, 2, as.integer)

starting_nodes <- subset(unique_nodes, unique_nodes %in% unique(solutions[,1]))

number_iterations <- as.numeric(iterations)

count_left <- 1

count_right <- number_iterations

total_pairs <- NULL

for(iteration in 1:length(starting_nodes))
{

  node_sol <- co_ocurrence[c(count_left:count_right),]
  sumatory <- colSums(node_sol)
  df_test <- data.frame(Starting_node = solutions[count_left,1], 
                        Connecting_node = colnames(node_sol),
                        weight = sumatory)
  total_pairs <- rbind(total_pairs, df_test)
  
  count_left <- count_left + number_iterations
  count_right <- count_right + number_iterations
  
  
  }

total_pairs <- subset(total_pairs, total_pairs$Starting_node %in% starting_nodes)
total_pairs <- subset(total_pairs, total_pairs$Connecting_node %in% starting_nodes)
total_pairs <- subset(total_pairs, as.character(total_pairs$Connecting_node) != as.character(total_pairs$Starting_node))
total_pairs <- subset(total_pairs, total_pairs$weight > 1)

## Looking back in the solutions if there are circular graphs

circular_sequences <- NULL
for(solution in 1:nrow(solutions))
{
  iteration <- solutions[solution,]
  iteration <- t(iteration)
  nodes <- as.character(iteration[,1])
  nodes <- nodes[nodes != '']
  if(length(nodes) > 1 & length(nodes) < 20 & nodes[1] == nodes[length(nodes)])
  {
    circular_sequences <- rbind(circular_sequences, c(nodes[1],nodes[length(nodes)]))
  }
}

if(is.null(circular_sequences) == FALSE)
{
  
}

no_duplicated <- circular_sequences[!duplicated(circular_sequences),]
for(combination in 1:nrow(no_duplicated))
{
  combi <- no_duplicated[combination,]
  total_ocurrences <- subset(circular_sequences, circular_sequences[,2] == combi[2])
  if(nrow(total_ocurrences) == number_iterations)
  {
    df_test <- data.frame(Starting_node = combi[1],
                          Connecting_node = combi[2],
                          weight = nrow(total_ocurrences))
    
    total_pairs <- rbind(total_pairs, df_test)
  }
  
}

total_pairs$Starting_node <- gsub(pattern = '\\+', replacement = '', x = total_pairs$Starting_node)
total_pairs$Starting_node <- gsub(pattern = '\\-', replacement = '', x = total_pairs$Starting_node)

total_pairs$Connecting_node <- gsub(pattern = '\\+', replacement = '', x = total_pairs$Connecting_node)
total_pairs$Connecting_node <- gsub(pattern = '\\-', replacement = '', x = total_pairs$Connecting_node)


graph_pairs <- graph_from_data_frame(total_pairs, directed = FALSE)

#E(graph_pairs)$width <- E(graph_pairs)$weight

V(graph_pairs)$name <- names(as.table(V(graph_pairs)))

graph_viz <- ggraph(graph_pairs, layout = 'nicely') + 
  geom_edge_link() +
  geom_node_point(size = 1.0) +
  geom_node_text(aes(label = name), size = 7.0) +
  geom_edge_loop()

#mcl_input <- as_adj(graph_pairs,attr="weight")

mcl_input <- as_adj(graph_pairs)


results_mcl <- mcl(x = mcl_input, addLoops = FALSE, allow1 = TRUE)

results_subgraph <- data.frame(number = rownames(mcl_input),
                               Component = results_mcl$Cluster)


suppressMessages(ggsave(filename = snakemake@output[["plot_graph"]], plot = graph_viz))


pl_nodes <- subset(clean_pred, clean_pred$Prediction == 'Plasmid' & clean_pred$Prob_Plasmid > as.numeric(as.character(threshold))) # Selecting only contigs predicted as plasmid-derived 

raw_number <- str_split_fixed(string = pl_nodes$Contig_name, pattern = '_', n = 2)[,1]
pl_nodes$number <- gsub(pattern = 'S', replacement = '', x = raw_number)


pl_notassigned <- subset(pl_nodes,! pl_nodes$number %in% results_subgraph$number)

pl_repeats <- subset(pl_notassigned, pl_notassigned$number %in% repeats$number)

pl_unassigned <- subset(pl_notassigned,! pl_notassigned$number %in% repeats$number)

pl_assigned <- subset(pl_nodes, pl_nodes$number %in% results_subgraph$number)
full_info_assigned <- merge(pl_assigned, results_subgraph, by = 'number')

if(nrow(pl_unassigned) > 1)
{
  pl_unassigned$Component <- 'Unbinned'
  full_info_assigned <- rbind(full_info_assigned, pl_unassigned)
}

if(nrow(pl_repeats) > 1)
{
  pl_repeats$Component <- 'Repeat-like'
  full_info_assigned <- rbind(full_info_assigned, pl_repeats)
}

full_info_assigned$Contig_length <- NULL
full_info_assigned$Prob_Chromosome <- round(full_info_assigned$Prob_Chromosome,2)
full_info_assigned$Prob_Plasmid <- round(full_info_assigned$Prob_Plasmid,2)
full_info_assigned$coverage <- round(full_info_assigned$coverage,2)


assembly_nodes <- readDNAStringSet(filepath = path_nodes)
df_nodes <- data.frame(Contig_name = names(assembly_nodes), Sequence = paste(assembly_nodes))

df_nodes <- merge(df_nodes, full_info_assigned, by = 'Contig_name')

for(component in unique(df_nodes$Component))
{
  nodes_component <- subset(df_nodes, df_nodes$Component == component)
  component_complete_name <- paste(snakemake@params[["sample"]], 'component', sep = '_')
  filename <- paste('results/', component_complete_name, sep = '')
  filename <- paste(filename,component, sep = '_')
  filename <- paste(filename,'.fasta',sep = '')
  write.fasta(sequences = as.list(nodes_component$Sequence), names = nodes_component$Contig_name, file.out = filename)

}

suppressWarnings(write.table(x = full_info_assigned, file = snakemake@output[["results"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
suppressWarnings(write.table(x = results_subgraph, file = snakemake@output[["components"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))





