#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 

library(igraph)
library(ggraph)
library(mlplasmids)
library(tidyverse)
library(spatstat)
library(cooccur)
library(cowplot)
library(ggrepel)

# Inputs

path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["links"]]
path_prediction <- snakemake@input[["prediction"]]
path_graph_contigs <- snakemake@input[["graph_contigs"]]
path_graph_repeats <- snakemake@input[["graph_repeats"]]
path_init_nodes <- snakemake@input[["initialize_nodes"]]
path_cov_variation <- snakemake@input[["coverage"]]

graph_contigs <- read.table(file = path_graph_contigs)
graph_repeats <- read.table(file = path_graph_repeats)
initialize_nodes <- read.table(file = path_init_nodes, header = TRUE)
initialize_nodes <- initialize_nodes[,1]

max_variation <- read.table(file = path_cov_variation, header = TRUE)
max_variation <- as.numeric(max_variation[1,1])
max_variation_small <- max_variation*5.0

#Outputs

output_path <- snakemake@output[["solutions"]]
output_connections <- snakemake@output[["connections"]]

plasmid_graph <- function(path_nodes, path_links, output_path, classifier, verbose = TRUE, number_iterations = 1, number_nodes = 10, initial_seed, max_variation = max_variation, prob_small_repeats = prob_small_repeats)
{
  initial_links <- links #Links to start the function 
  
  for(iterations in 1:number_iterations) # Number of times we repeat this process
  {
    path <- NULL # We declare an empty vector
    path <- append(x = path, values = initial_seed, after = length(path)) # We add the initial seed to the path, first element in the vector
    
    if(verbose == TRUE)
    {
      print('INITIALISING THE PATH') # Adding extra information, useful to debug and observe the behaviour of the function 
    }
    
    seed <- initial_seed 
    links <- initial_links
    
    record_connections <- NULL
    
    record_connections <- data.frame(factor = factor_mad,
                                     number_iterations = number_iterations,
                                     iteration = iterations,
                                     elongation = 0, 
                                     first_node = initial_seed, 
                                     ingoing_node = seed, 
                                     outgoing_node = seed, 
                                     Probability_pl_chr = 1,
                                     Probability_cov = 1,
                                     Probability = 1,
                                     Probability_freq = 1,
                                     Veredict = 'selected')
    
    write.table(x = record_connections, 
                file = output_connections, 
                append = TRUE, 
                row.names = FALSE, 
                quote = FALSE, 
                col.names = FALSE)
    
    
    
    #################################### Coverage of the current path #################################
    
    # Extracting the info from our current path                     
    info_path <- graph_contigs[graph_contigs$number %in% path,]
    info_path <- info_path[! info_path$number %in% repeats_graph$number,] # Removing the contigs corresponding to transposases 
    length_path <- sum(info_path$length) # Length of the path
    info_path$contribution <- info_path$length/length_path # Shorter contigs should have a lower contribution to the coverage. k-mer coverage on these contigs fluctuates drastically
    
    path_mean <- weighted.mean(x = info_path$coverage, w = info_path$contribution) # Coverage of the current path 
    
    initial_path_mean <- path_mean
    
    ##################### Elongating the path ###########################################
    
    for(elongation in 1:number_nodes)
    {
      evaluation <- 'nope' # Default integer assigned to the evaluation, if different then exit the loop
      if(verbose == TRUE)
      {
        print(paste(elongation,'elongation:','possible connections'))
      }
      
      current_links <- subset(links, links$V1 == seed) # Consider the last element present in our path and observe all the possible links
      
      if(verbose == TRUE)
      {
        print(current_links) # Possible links
      }
      
      if(nrow(current_links) == 0)
      {
        if(verbose == TRUE)
        {
          print(path, quote = FALSE) 
        }
        
        output <- paste(path, collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed # There are no connections possible from this contig
        break # Exiting the loop 
      }
      
      list_connections <- unique(current_links$V3) # All the possible connections 
      
      # We do not allow that a node which is not a repeat appears more than 1 time in any solution but we exclude the initial seed from this consideration 
      
      if(length(path) != 1) # If the path has more than one element  
      {
        remove_nodes <- path[2:length(path)]
        
        remove_nodes <- gsub(x = remove_nodes, pattern = '\\+', replacement = '')
        remove_nodes <- gsub(x = remove_nodes, pattern = '\\-', replacement = '')
        
        positive_remove_nodes <- paste(remove_nodes, '+', sep = '')
        negative_remove_nodes <- paste(remove_nodes, '-', sep = '')
        
        ommit_nodes <- c(positive_remove_nodes, negative_remove_nodes)
        
        # We need to remove directionality from the path to avoid paths e.g. 54-,161-,54+
        
        ommit_nodes <- ommit_nodes[! ommit_nodes %in% repeats_graph$number] # Transposases can occurr more than one time 
        list_connections <- list_connections[! list_connections %in% ommit_nodes] # Avoiding loops within the solution
      }
      
      possible_connections <- as.character(list_connections) # Connections that can take place  
      total_connections <- length(possible_connections) # Number of connections
      base_probabilities <- rep(0, total_connections) # Generating a vector with the connection probabilities
      
      
      if(length(possible_connections) < 1)
      {
        output <- paste(path, collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed
        break 
      }
      
      ################################# Probabilities of the connections #####################################################
      
      # We generate a dataframe with the connection probabilities
      prob_df <- data.frame(number = possible_connections,
                            prob = base_probabilities)
      
      if(verbose == TRUE)
      {
        print(prob_df)
      }
      
      # Replacing the base probabilities with the probabilities generated by mlplasmids 
      
      prob_df$number <- gsub(pattern = '\\+',replacement = '',x = prob_df$number)
      prob_df$number <- gsub(pattern = '\\-',replacement = '',x = prob_df$number)
      prob_df$number <- as.character(prob_df$number)
      
      if(classifier == 'mlplasmids')
      {
        subset_probs <- mlplasmids_clean_pred[which(mlplasmids_clean_pred$number %in% prob_df$number),]
        final_probs <- subset_probs$Prob_Plasmid[match(prob_df$number,subset_probs$number)]
      }
      else
      {
        subset_probs <- plasflow_clean_prediction[which(plasflow_clean_prediction$number %in% prob_df$number),]
        final_probs <- subset_probs$Probability[match(prob_df$number,subset_probs$number)]
      }
      
      
      
      # Short contigs do not have a reliable probability, we assign them a predefined probability (passed via the argument 'prob_small_repeats')
      final_probs[is.na(final_probs)] <- prob_small_repeats # OVERLAP!  
      
      
      # Transposases are also corner-cases for the machine-learning algorithm, we follow the same principle as done with short contigs 
      index_repeats <- which(prob_df$number %in% repeats$number)
      final_probs[index_repeats] <- prob_small_repeats
      
      if(verbose == TRUE)
      {
        print(final_probs)
      }
      
      # Inner loop to evaluate all the sensible connections
      
      times_sampling <- 0
      
      
      # From the connections we sample a possible connection based on the vector with the probabilities
      
      record_connections <- NULL
      
      record_connections <- data.frame(factor = factor_mad,
                                       number_iterations, number_iterations,
                                       iteration = iterations,
                                       elongation = elongation, 
                                       first_node = initial_seed, 
                                       ingoing_node = seed, 
                                       outgoing_node = list_connections, 
                                       Probability_pl_chr = final_probs,
                                       Probability_cov = 0)
      
      cov_connections_info <- graph_contigs[graph_contigs$number %in% record_connections$outgoing_node,]
      
      up_cutoff <- cov_connections_info$coverage + max_variation
      down_cutoff <- cov_connections_info$coverage - max_variation
      
      up_threshold <- pnorm(up_cutoff, mean = path_mean, sd = 0.1 , lower.tail = TRUE)
      
      down_threshold <- pnorm(down_cutoff, mean = path_mean, sd = 0.1 , lower.tail= FALSE)
      
      window <- up_threshold - down_threshold
      
      cov_connections_info$Probability_cov <- 1-abs(window)
      
      record_connections$Probability_cov <- cov_connections_info$Probability_cov[match(record_connections$outgoing_node, cov_connections_info$number)]
      
      record_connections$Probability_cov[which(record_connections$outgoing_node %in% repeats_graph$number)] <- prob_small_repeats
      
      record_connections$Probability_cov[which(record_connections$outgoing_node %in% small_contigs$number)] <- prob_small_repeats
      
      
      
      record_connections$Probability <- record_connections$Probability_pl_chr * record_connections$Probability_cov
      
      record_connections$Probability_freq <- record_connections$Probability/(sum(record_connections$Probability))
      
      record_connections$Veredict <- 'non-selected'
      
      
      if(length(which(record_connections$Probability > 0.1)) == 0)
      {
        output <- paste(path, collapse = ',')
        write.table(x = output, 
                    file = output_path, 
                    append = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE,
                    col.names = FALSE)
        
        path <- initial_seed
        
        write.table(x = record_connections, 
                    file = output_connections, 
                    append = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE, 
                    col.names = FALSE)
        
        break 
      }
      
      # Filter step to avoid going into really bad connections 
      
      filter_connections <- subset(record_connections, record_connections$Probability > 0.1)
      
      random_connection <- sample(x = filter_connections$outgoing_node, size = 1, prob = filter_connections$Probability) # Choose one connection 
      
      record_connections$Veredict[which(record_connections$outgoing_node == random_connection)] <- 'selected'
      
      write.table(x = record_connections, 
                  file = output_connections, 
                  append = TRUE, 
                  row.names = FALSE, 
                  quote = FALSE, 
                  col.names = FALSE)
      
      
      if(random_connection == path[1])
      {
        path <- append(values = as.character(random_connection), after = length(path), x = path)
        output <- paste(path, collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed
        break 
      }
      
      
      path <- append(values = as.character(random_connection), after = length(path), x = path)
      seed <- as.character(random_connection)  
      
      if(length(path) == number_nodes) # We only exit the function if we have reached the maximum number of nodes allowed per path
      {
        output <- paste(path, collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed
        break 
      }
      
      
      info_path <- graph_contigs[graph_contigs$number %in% path,]
      
      info_path <- info_path[! info_path$number %in% repeats_graph$number,]
      
      length_path <- sum(info_path$length)
      info_path$contribution <- info_path$length/length_path
      
      path_mean <- weighted.mean(x = info_path$coverage, w = info_path$contribution)
      
      if(verbose == TRUE)
      {
        print(paste('in this',elongation,'elongation the path finishes like this:'))
        print(path)
        print(paste('and the new path has the updated following coverage:',path_mean))
        cat("\n")
      }
      
      
    }
    
  }
}

for(seed in initialize_nodes)
{
  set.seed(123)
  positive_seed <- paste(seed, '+', sep = '')
  negative_seed <- paste(seed, '-', sep = '')
  
  plasmid_graph(path_nodes = path_nodes, path_links =  path_links, output_path = output_path, initial_seed = positive_seed, number_iterations = 10, verbose = TRUE,  number_nodes = 20, prob_small_repeats = sqrt(0.6), max_variation = max_variation, max_variation_small = max_variation_small, classifier = 'mlplasmids')
  plasmid_graph(path_nodes = path_nodes, path_links =  path_links, output_path = output_path, initial_seed = negative_seed, number_iterations = 10, verbose = TRUE,  number_nodes = 20, prob_small_repeats = sqrt(0.6), max_variation = max_variation, max_variation_small = max_variation_small, classifier = 'mlplasmids')
}  



