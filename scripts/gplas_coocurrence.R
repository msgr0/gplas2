#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(cooccur))
suppressMessages(library(ggrepel))

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
co_ocurrence <- t(co_ocurrence)

ocurrence_large_nodes <- co_ocurrence[rownames(co_ocurrence) %in% unique(solutions[,1]),]
ocurrence_large_nodes <- as.matrix(ocurrence_large_nodes)

suppressMessages(test <- cooccur(mat = ocurrence_large_nodes, type = "spp_site", thresh = FALSE, spp_names = TRUE))

df_combinations <- prob.table(test)

significant_associations <- subset(df_combinations, df_combinations$p_gt < 0.05/nrow(df_combinations))

significant_associations$sp1_name <- gsub(pattern = '\\+',replacement = '', significant_associations$sp1_name)
significant_associations$sp1_name <- gsub(pattern = '\\-',replacement = '', significant_associations$sp1_name)

significant_associations$sp2_name <- gsub(pattern = '\\+',replacement = '', significant_associations$sp2_name)
significant_associations$sp2_name <- gsub(pattern = '\\-',replacement = '', significant_associations$sp2_name)


significant_associations <- significant_associations[! significant_associations$sp1_name %in% repeats$number,]
significant_associations <- significant_associations[! significant_associations$sp2_name %in% repeats$number,]

dimensions_sign <- dim(significant_associations)
  
  if(is.null(dimensions_sign) == TRUE)
  {
    print('No significant associations')
    break
  }
# Looking back in the solutions if there are circular graphs

circular_sequences <- NULL

  for(solution in 1:nrow(solutions))
  {
    iteration <- solutions[solution,]
    iteration <- t(iteration)
    nodes <- as.character(iteration[,1])
    nodes <- nodes[nodes != '']
    if(length(nodes) > 1 & length(nodes) < 4 & nodes[1] == nodes[length(nodes)])
    {
      circular_sequences <- rbind(circular_sequences, c(nodes[1],nodes[length(nodes)]))
    }
  }

  if(is.null(circular_sequences) == TRUE)
  {
    df_graph <- data.frame(from = significant_associations$sp1_name, to = significant_associations$sp2_name)
  }
  if(is.null(circular_sequences) == FALSE)
  {
    circular_sequences <- data.frame(from = circular_sequences[,1], to = circular_sequences[,2])
    circular_sequences <- circular_sequences[! duplicated(circular_sequences),]
    df_graph <- data.frame(from = significant_associations$sp1_name, to = significant_associations$sp2_name)
    df_graph <- rbind(df_graph, circular_sequences)
  }

# Removing directionality from the graph 

df_graph$from <- gsub(pattern = '\\+',replacement = '', df_graph$from)
df_graph$from <- gsub(pattern = '\\-',replacement = '', df_graph$from)
df_graph$to <- gsub(pattern = '\\+',replacement = '', df_graph$to)
df_graph$to <- gsub(pattern = '\\-',replacement = '', df_graph$to)


hairball <- graph_from_data_frame(df_graph)

V(hairball)$name <- names(as.table(V(hairball)))

graph_viz <- ggraph(hairball, layout = 'nicely') + 
  geom_edge_link() +
  geom_node_point(size = 1.0) +
  geom_node_text(aes(label = name), size = 7.0) +
  geom_edge_loop()

suppressMessages(ggsave(filename = snakemake@output[["plot_graph"]], plot = graph_viz))

# Retrieving the clustering of the contigs into different components 

results_subgraph <- data.frame(number = labels(components(hairball)$membership),
                               Component = as.numeric(as.character(components(hairball)$membership)))

results_subgraph <- rbind(results_subgraph, results_subgraph[which(results_subgraph$number %in% circular_sequences[,1]),])

pl_nodes <- subset(clean_pred, clean_pred$Prediction == 'Plasmid' & clean_pred$Prob_Plasmid > as.numeric(as.character(threshold))) # Selecting only contigs predicted as plasmid-derived 

pl_notassigned <- subset(pl_nodes,! pl_nodes$number %in% results_subgraph$number)


pl_repeats <- subset(pl_notassigned, pl_notassigned$number %in% repeats$number)

pl_unassigned <- subset(pl_notassigned,! pl_notassigned$number %in% repeats$number)

pl_assigned <- subset(pl_nodes, pl_nodes$number %in% results_subgraph$number)
full_info_assigned <- merge(pl_assigned, results_subgraph, by = 'number')

if(exists('pl_unassigned') == FALSE)
{
  pl_unassigned$Component <- 'Unbinned'
  full_info_assigned <- rbind(full_info_assigned, pl_unassigned)
}

if(exists('pl_repeats') == FALSE)
{
  pl_repeats$Component <- 'Repeat-like'
  full_info_assigned <- rbind(full_info_assigned, pl_repeats)
}

full_info_assigned$Contig_length <- NULL
full_info_assigned$Prob_Chromosome <- round(full_info_assigned$Prob_Chromosome,2)
full_info_assigned$Prob_Plasmid <- round(full_info_assigned$Prob_Plasmid,2)
full_info_assigned$coverage <- round(full_info_assigned$coverage,2)

suppressWarnings(write.table(x = full_info_assigned, file = snakemake@output[["results"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
suppressWarnings(write.table(x = results_subgraph, file = snakemake@output[["components"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
