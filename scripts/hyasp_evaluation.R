# hyasp evaluation 

#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(cooccur))
suppressMessages(library(ggrepel))

# Inputs

path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["clean_links"]]
path_graph_contigs <- snakemake@input[["graph_contigs"]]
path_graph_repeats <- snakemake@input[["graph_repeats"]]
path_cov_variation <- NA
path_alignments <- snakemake@input[["alignments"]]
path_circularity <- snakemake@input[["circularity"]]

path_contig_chains <- snakemake@input[["chains"]]
path_contig_bins <- snakemake@input[["bins"]]

classifier <- snakemake@params[["classifier"]]
number_iterations <- snakemake@params[["iterations"]]
species <- snakemake@params[["species"]]
sample <- snakemake@params[["name"]]


# Use to debug 
# 
 # path_nodes <- '/home/sergi/gplas/hyasp/E0139_hyasp_raw_nodes.fasta'
 # path_links <- '/home/sergi/gplas/hyasp/E0139_hyasp_clean_links.tab'
 # path_graph_contigs <- '/home/sergi/gplas/hyasp/E0139_hyasp_graph_contigs.tab'
 # path_graph_repeats <- '/home/sergi/gplas/hyasp/E0139_hyasp_repeats_graph.tab'
 # path_alignments <- '/home/sergi/gplas/hyasp/E0139_hyasp_alignment_test.txt'
 # path_contig_bins <- '/home/sergi/gplas/hyasp/E0139_hyasp_plasmid_bins_putative.csv'
 # path_contig_chains <- '/home/sergi/gplas/hyasp/E0139_hyasp_contig_chains.csv'
 # path_circularity <- '/home/sergi/gplas/hyasp/E0139_hyasp_circular_seq.txt'
 # classifier <- 'hyasp'
 # number_iterations <- NA
 # species <- 'Enterococcus faecium'
 # sample <- 'E0139'


# Recovering some info about the graph and contigs
max_variation <- NA

links <- read.table(file = path_links, header = TRUE)
graph_contigs <- read.table(file = path_graph_contigs, header = TRUE)

small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

repeats_graph <- read.table(file = path_graph_repeats)
repeats <- repeats_graph
repeats$number <- gsub(pattern = '\\+',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)
repeats$number <- gsub(pattern = '\\-',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)

# First reading the file with all the contig chains 'contig_chains.csv'
contig_chains <- read.table(file = path_contig_chains, sep = ';')


chains_info <- NULL
for(component in unique(contig_chains$V1))
{
  info_component <- subset(contig_chains, contig_chains$V1 == component)
  contigs_component <- as.character(info_component$V2)
  list_contigs <- unlist(strsplit(contigs_component, ","))
  for(contig in list_contigs)
  {
    contig_info <- data.frame(Contig_number = contig,
               Plasmid = component)
    chains_info <- rbind(chains_info, contig_info)
  }

}

putative_plasmid_bins <- read.table(file = path_contig_bins)

plasmid_bins <- NULL
bin_number <- 1
for(bin in unique(putative_plasmid_bins$V1))
{
  bin_info <- unlist(strsplit(bin, ","))
  for(plasmid in bin_info)
  {
    print(plasmid)
    associating_plasmids_bin <- data.frame(Plasmid = plasmid,
               Bin = bin_number )
    plasmid_bins <- rbind(plasmid_bins, associating_plasmids_bin)
  }
  bin_number <- bin_number + 1
}

complete_hyasp_info <- merge(chains_info, plasmid_bins, by = 'Plasmid')

# Reading if some of the components with a small number of contigs are actually predicted as circular by hyasp 

circularity <- read.table(path_circularity)
circularity <- subset(circularity, circularity$V8 == 'circular=1')
circularity$V1 <- gsub(pattern = '>', replacement = '', x = circularity$V1)

# Second reading the results of quast to observe which contigs belong to each plasmid unit 

alignment_quast <- read.table(file = path_alignments, sep = ' ', header = TRUE)

alignment_quast$Reference_number <- str_split_fixed(alignment_quast$Reference, "_", 6)[,1] # Extracting the reference of the number
alignment_quast$Reference_length <- as.numeric(str_split_fixed(alignment_quast$Reference, "_", 6)[,3]) # Extracting the length of the reference genome

alignment_quast$Contig_number <- str_split_fixed(alignment_quast$Contig, "_", 6)[,1] # Extracting the contig number
alignment_quast$Contig_number <- gsub(pattern = 'S', replacement = '', x = alignment_quast$Contig_number)

alignment_quast$Contig_length <- as.numeric(str_split_fixed(alignment_quast$Contig, "_", 6)[,4]) # Extracting the contig number

alignment_quast$Type <- ifelse(alignment_quast$Reference_length > 4e5, 'Chromosome','Plasmid')
alignment_quast$Type <- paste(alignment_quast$Type,alignment_quast$Reference_number, sep = '')

colnames(alignment_quast)[2] <- 'Contig_name'

# Remove duplicates present in the alignment (e.g. based on the fact that plasmids are circular)

non_duplicates_quast <<- alignment_quast[! duplicated(alignment_quast),]

non_transposases_quast <<- non_duplicates_quast[! non_duplicates_quast$Contig_number %in% repeats$number,]

gold_standard <- non_transposases_quast

# Including information about the number of contigs mapping to the reference genome and the total number of base-pairs 

gold_standard <- gold_standard %>%
  group_by(Type) %>%
  mutate(sum(Contig_length),
         count = n())

colnames(gold_standard)[c(8,9)] <- c('Reference_sum_bp','Reference_sum_contigs')
evaluation_ref_genomes <- gold_standard

# Dataframe with all the stats regarding the reference genomes (Number of contigs and number of base-pairs)
stats_gold_standard <- gold_standard %>% 
  group_by(Type) %>%
  summarise(sum(Contig_length),
            count = n())# Adding the contig direction into the results (we use this information to merge our gold_standard with the prediction) 

colnames(stats_gold_standard) <- c('Type','Reference_sum_bp','Reference_sum_contigs')


# Creating a dataframe to merge the quast mapping against the prediction from plasgraph. We need to include the directionality of the contigs 

contigs_reference_genome <- gold_standard %>%
  group_by(Type) %>%
  summarise(count = n(),
            sum(Contig_length))

colnames(contigs_reference_genome) <- c('Type','Reference_sum_contigs','Reference_sum_bp')


truth_set <- subset(gold_standard, select = c('Contig_number','Type'))
truth_set <- truth_set[order(truth_set$Type, decreasing = TRUE),]

truth_set <- truth_set[! duplicated(truth_set$Contig_number),]

# Merging prediction together with the truth set that we defined using quast

complete_hyasp_info$Contig_number <- gsub(pattern = '\\+', replacement = '', x = complete_hyasp_info$Contig_number)
complete_hyasp_info$Contig_number <- gsub(pattern = '-', replacement = '', x = complete_hyasp_info$Contig_number)


results <- merge(complete_hyasp_info, gold_standard, by = 'Contig_number')

results <- results[!duplicated(results),]

circular_results <- subset(results, results$Plasmid %in% circularity$V1)

benchmark <- results %>%
  group_by(Bin, Type) %>%
  summarise(count = n(), sum_comp = sum(Contig_length))

benchmark$Component <- benchmark$Bin


sort_benchmark <- benchmark[order(benchmark$Component),]
sort_benchmark <- sort_benchmark[order(as.numeric(sort_benchmark$count), decreasing = TRUE),]

total_component_info <- sort_benchmark

reference_component <- sort_benchmark[! duplicated(sort_benchmark$Component),]
reference_component <- reference_component[order(reference_component$Component, decreasing = FALSE),]
completeness_df <- NULL

for(component in unique(reference_component$Component))
{
  comp_info <- subset(reference_component, reference_component$Component == component)
  total_component <- subset(total_component_info, total_component_info$Component == component)
  
  chr_component <- subset(total_component, total_component$Type == 'Chromosome1')
  
  if(nrow(chr_component) > 0)
  {
    chr_contigs <- chr_component$count
    chr_bp <- chr_component$sum_comp
  }
  else
  {
    chr_contigs <- 0
    chr_bp <- 0
  }
  
  
  total_component_contigs <- sum(as.numeric(total_component$count))
  total_component_bp <- sum(as.numeric(total_component$sum_comp))
  
  if(component == 'No_component')
  {
    completeness_info <- data.frame(sample = sample,
                                    species = species,
                                    classifier = classifier, 
                                    factor = 1.0,
                                    number_iterations = number_iterations,
                                    variation = max_variation,
                                    component = component,
                                    total_contigs = comp_info$count,
                                    total_bp = comp_info$sum_comp,
                                    reference = 0,
                                    contigs_ref_component = 0,
                                    bp_ref_component = 0,
                                    contigs_reference = 0,
                                    bp_reference = 0,
                                    completeness_contigs = 0,
                                    completeness_bp = 0,
                                    chr_contigs = 0,
                                    chr_bp = 0)
  }
  else
  {
    
    count_truth <- subset(contigs_reference_genome$Reference_sum_contigs, contigs_reference_genome$Type == comp_info$Type)
    count_truth_bp <- subset(contigs_reference_genome$Reference_sum_bp, contigs_reference_genome$Type == comp_info$Type)
    
    completeness_contigs <- as.numeric(as.character(comp_info$count))/count_truth
    completeness_bp <- as.numeric(as.character(comp_info$sum_comp))/count_truth_bp
    
    
    completeness_info <- data.frame(sample = sample,
                                    species = species,
                                    classifier = classifier, 
                                    factor = 1.0,
                                    number_iterations = number_iterations,
                                    variation = max_variation,
                                    component = component,
                                    total_contigs = as.numeric(as.character(total_component_contigs)),
                                    total_bp = as.numeric(as.character(total_component_bp)),
                                    reference = comp_info$Type,
                                    contigs_ref_component = comp_info$count,
                                    bp_ref_component = comp_info$sum_comp,
                                    contigs_reference = count_truth,
                                    bp_reference = count_truth_bp,
                                    completeness_contigs = completeness_contigs,
                                    completeness_bp = completeness_bp,
                                    chr_contigs = chr_contigs,
                                    chr_bp = chr_bp)
    
  }
  
  completeness_df <- rbind(completeness_df, completeness_info)
  
}

multiple_contigs_comp <- subset(completeness_df, completeness_df$total_contigs > 1)
single_contigs_comp <- subset(completeness_df, completeness_df$total_contigs == 1)
single_contigs_eval_comp <- subset(single_contigs_comp, single_contigs_comp$completeness_contigs > 0.5)

circular_contigs_comp <- subset(single_contigs_comp, single_contigs_eval_comp$component %in% circular_results$Bin)

completeness_df <- rbind(multiple_contigs_comp, circular_contigs_comp)

write.table(x = completeness_df, 
            file = snakemake@output[["completeness"]],
            append = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            col.names = FALSE)


all_truth <- NULL

for(seed in truth_set$Contig_number)
{
  contig_seed <<- subset(truth_set, truth_set$Contig_number == seed)
  other_contigs <<- subset(truth_set, truth_set$Contig_number != seed)
  
  evaluation <- data.frame(Contig_number = seed,
                           Pair_contig = other_contigs$Contig_number,
                           Contig_number_type = contig_seed$Type,
                           Pair_contig_type = other_contigs$Type)
  
  
  evaluation$Truth <- ifelse(as.character(evaluation$Contig_number_type) == as.character(evaluation$Pair_contig_type), 'Same','Different')
  
  all_truth <- rbind(all_truth, evaluation)
}

all_truth$Contig_pair <- paste(all_truth$Contig_number, all_truth$Pair_contig, sep = '-')


## Extracting the number of connections from each plasmid and chromosome

true_connections <- subset(all_truth, all_truth$Contig_number_type == all_truth$Pair_contig_type)

count_true_connections <- true_connections %>%
  group_by(Contig_number_type) %>%
  count()

duplicated_hyasp_results <- complete_hyasp_info[duplicated(complete_hyasp_info$Contig_number),]
complete_hyasp_info <- subset(complete_hyasp_info,! complete_hyasp_info$Contig_number %in% duplicated_hyasp_results$Contig_number)

all_evaluation <- NULL

for(seed in complete_hyasp_info$Contig_number)
{
  contig_seed <- subset(complete_hyasp_info, complete_hyasp_info$Contig_number == seed)
  other_contigs <- subset(complete_hyasp_info, complete_hyasp_info$Contig_number != seed)
  
  evaluation <- data.frame(Contig_number = seed,
                           Pair_contig = other_contigs$Contig_number,
                           Contig_number_component = contig_seed$Bin,
                           Pair_contig_component = other_contigs$Bin)
  
  
  evaluation$Prediction <- ifelse(as.character(evaluation$Contig_number_component) == as.character(evaluation$Pair_contig_component), 'Same','Different')
  
  all_evaluation <- rbind(all_evaluation, evaluation)
}

all_evaluation$Contig_pair <- paste(all_evaluation$Contig_number, all_evaluation$Pair_contig, sep = '-' )


# Creating a confusion matrix 

eval_truth <<- merge(all_evaluation, all_truth, by = 'Contig_pair')


eval_truth$Evaluation <- ifelse(eval_truth$Truth == eval_truth$Prediction, 'Positive', 'Negative')

positive_eval_truth <- subset(eval_truth, eval_truth$Evaluation == 'Positive')
positive_eval_truth$Evaluation <- ifelse(positive_eval_truth$Truth == 'Same', 'True_positive','True_negative')

negative_eval_truth <- subset(eval_truth, eval_truth$Evaluation == 'Negative')
negative_eval_truth$Evaluation <- ifelse(negative_eval_truth$Truth == 'Different', 'False_positive','False_negative')

eval_truth <- rbind(positive_eval_truth, negative_eval_truth)

purity_components <- eval_truth %>% 
  group_by(eval_truth$Contig_number_component, eval_truth$Evaluation) %>%
  count(Pair_contig_component)

purity_components <- subset(purity_components, purity_components$`eval_truth$Evaluation` %in% c('True_positive','False_positive'))

precision_components <- NULL

for(component in unique(purity_components$Pair_contig_component))
{
  component_df <- subset(purity_components, purity_components$Pair_contig_component == component)
  test <<- subset(component_df, component_df$`eval_truth$Evaluation` == 'True_positive')
  
  tp_component <- subset(component_df$n, component_df$`eval_truth$Evaluation` == 'True_positive')
  fp_component <- subset(component_df$n, component_df$`eval_truth$Evaluation` == 'False_positive')
  
  if(identical(fp_component, integer(0)))
  {
    fp_component <- 0
  }
  
  if(identical(tp_component, integer(0)))
  {
    tp_component <- 0
  }
  
  precision <- tp_component/(tp_component+fp_component)
  
  record_results <- data.frame( sample = sample,
                                species = species,
                                classifier = classifier, 
                                factor = 1.0,
                                number_iterations = number_iterations,
                                variation = max_variation,
                                component = component, 
                                tp_component = tp_component,
                                fp_component = fp_component,
                                precision = precision)
  
  
  precision_components <- rbind(precision_components, record_results)
  
}

write.table(x = precision_components, 
            file = snakemake@output[["precision"]],
            append = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            col.names = FALSE)