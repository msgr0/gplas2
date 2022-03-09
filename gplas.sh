#!/bin/bash

##to debug
#set -e
#set -v
#set -x

########################################
## A bash script to run the gplas pipeline
## This script has been influenced by and adapted from the code in 'bactofida', by aschuerch:
## https://gitlab.com/aschuerch/bactofidia

version="1.0.0"

## Help message
usage() {
cat <<HELP_USAGE
Welcome to the user guide of gplas (version ${version}).

BASIC USAGE
  ./$(basename ${0}) [-i <file>] [-c <string>] [-s <string>] [...]

Example:
  ./$(basename ${0}) -i mygraph.gfa -c mlplasmids -s 'Enterococcus faecium'

USER OPTIONS
Input:
  -i  [MANDATORY] Path to the graph file in GFA (.gfa) format, used to extract nodes and links. 
  -c  [MANDATORY] Classifier used to predict the contigs extracted from the input graph.
      String value: 'plasflow', 'mlplasmids', 'extract' or 'predict'.
      Values 'extract' and 'predict' are used in the case of external binary classifiers, according to the steps:
        1) Classifier option must be set to 'extract'. This will generate a fasta file containing the extracted nodes
        sequences, which will be saved to gplas_input/\${name}_raw_nodes.fasta.
        2) Extracted nodes will be used as input for the binary classification tool selected by the user.
        After binary classifcation, predictions will need to be formated appropriately and saved to
        independent_prediction/\${name}_plasmid_prediction.tab.
        3) Classifier option must be set to 'predict'
  -s  [Optional; MANDATORY for 'mlplasmids' classifier] Bacterial species from the graph file.
      If you have specified mlplasmids as classifier you need to provide one of the following three bacterial species:
      'Enterococcus faecium', 'Enterococcus faecalis', 'Klebsiella pneumoniae', 'Acinetobacter baumannii' or
      'Escherichia coli'.

Output:
  -n  [Optional] Output name used in the files generated by gplas. Default: 'unnamed'.
  -k  [Optional] Keeps intermediary files (i.e. plasmid-walks).

Parameters:
  -t  [Optional] Threshold to predict plasmid-derived sequences. Integer value ranging from 0 to 1. 
      Only valid for classifier options 'plasflow' (default: 0.5) and 'mlplasmids' (default: 0.7).
  -b  [Optional] Coverage variance allowed for bold walks to recover unbinned plasmid-predicted nodes.
      Numeric value: X times coverage variance of the chromsome. Default: 5
  -x  [Optional] Number of times gplas finds plasmid walks per each plasmid starting node.
      Integer value ranging from 1 to infinite. Default: 20
  -f  [Optional] Gplas filtering threshold score to reject possible outcoming edges.
      Integer value ranging from 0 to 1. Default: 0.1
  -q  [Optional] Modularity threshold to split components present in the plasmidome network.
      Integer value ranging from 0 to 1. Default: 0.2

Other:
  -h  Print this help message.
  -v  Print gplas version.

HELP_USAGE
}

## Argument parsing
while getopts ":i:n:s:c:t:x:f:e:q:b:khv" opt; do
  case $opt in
    h)
      cat figures/logo.txt
      echo -e ""
      usage
      exit 0
      ;;
    i)
      input=$OPTARG
      ;;
    c)
      classifier=$OPTARG
      ;;
    n)
      name=$OPTARG
      ;;
    k)
      keep='true'   
      ;;
    s)
      species=$OPTARG
      ;;
    t)
      threshold_prediction=$OPTARG
      ;;
    x)
      number_iterations=$OPTARG
      ;;
    f)
      filt_gplas=$OPTARG
      ;;
    e)
      edge_gplas=$OPTARG
      ;;
    q)
      modularity_threshold=$OPTARG
      ;;
    b)
      bold_sd_coverage=$OPTARG
      ;;
    v)
      echo -e "gplas version ${version}"
      exit 0
      ;;
    \?)
      usage
      echo -e "Invalid option: -$OPTARG\n  gplas requires a single input corresponding to a GFA file.\n" >&2
      exit 1
      ;;
    :)
      usage
      echo -e "Error: Option -$OPTARG requires an argument.\n" >&2
      exit 1
      ;;
    esac
done

## Test conditions
if [ -z "$input" ]; then
  usage
  echo -e "\n Error: Oops, it seems that you are missing the input graph.\n"
  exit 1
fi

# better classifier test conditions
if [ -z "$classifier" ]; then
  # No classifier was chosen; throw error and exit
  usage
  echo -e "\n Error: Oops, it seems that you are missing the classifier that you want to use.\n"
  exit 1
elif [ "$classifier" == "mlplasmids" ]; then
  # set snakefile for analysis
  snakeFile="snakefiles/mlplasmidssnake.smk"
  # assert that species should be from list of valid species
  list_species=(Enterococcus faecium Enterococcus faecalis Klebsiella pneumoniae Acinetobacter baumannii Escherichia coli);
  if [ -z species ]; then
    usage
    echo -e "\n Error: You have specified mlplasmids as classifier but you have not indicated one of the following three bacterial species:
    \t'Enterococcus faecium','Enterococcus faecalis', 'Klebsiella pneumoniae', 'Acinetobacter baumannii' or 'Escherichia coli'\n"
    exit 1
  elif [[ " "${list_species[@]}" " == *" "$species" "* ]]; then
    # do nothing
    :
  else
    usage
    echo -e "Ups! Something went wrong\n"
    echo -e "The provided species:" "$species" "is not included in mlplasmids. Valid species names are:\n"
    echo "${list_species[@]/%/,}"
    exit 1
  fi
elif [ "$classifier" == "plasflow" ]; then
  # set snakefile for analysis
  snakeFile="snakefiles/plasflowsnake.smk"
elif [ "$classifier" == "predict" ] || [ "$classifier" == "extract" ]; then
  # set snakefile for analysis
  snakeFile="snakefiles/otherclassifiers.smk"
else
  # classifier provided is from an invalid set of classifiers
  usage
  echo -e "Error: Classifier provided is invalid. Please use one of the correct set of classifier values:"
  echo -e "'plasflow', 'mlplasmids', 'extract' or 'predict'"
  exit 1
fi

## Set default parameters if they are not provided
if [ -z "$species" ] && [ "$classifier" != "mlplasmids" ]; then
  # if classifier NOT mlplasmids, then species can be set to 'unknown'
  # if classifier is mlplasmids, species must be provided from valid set
  species="unknown"
fi

if [ -z "$name" ]; then
  name="unnamed"
fi

if [ -z "$threshold_prediction" ]; then
  if [ "$classifier" == "plasflow" ]; then
    # default threshold for plasflow
    threshold_prediction=0.7
  else
    # default threshold for other classifiers
    threshold_prediction=0.5
  fi #else nothing?
fi

if [ -z "$filt_gplas" ]; then
  filt_gplas=0.1
fi

if [ -z "$edge_gplas" ]; then
  edge_gplas=0.1
fi

if [ -z "$number_iterations" ]; then
  number_iterations=20
fi

if [ -z "$modularity_threshold" ]; then
  modularity_threshold=0.2
fi

if [ -z "$bold_sd_coverage" ]; then
  bold_sd_coverage=5
fi

if [ -z "$reference" ]; then
  reference="No reference provided"
else
  mkdir -p reference_genome
  rm reference_genome/*
  cp $reference reference_genome/
  mv reference_genome/*.fasta reference_genome/"$name"_ref_genome.fasta
fi

########################################
## Start gplas
echo -e "\n"
cat figures/logo.txt
echo -e "\n"
echo "##################################################################"

## Print chosen parameters
cat <<START
Your results will be named ${name}
Input graph: ${input}
Bacterial species: ${species}
Classifier: ${classifier}
Threshold for predicting plasmid-derived contigs: ${threshold_prediction}
Number of plasmid walks created per node: ${number_iterations}
Threshold of gplas scores: ${filt_gplas}
Minimum frequency to consider an edge: ${edge_gplas}
Modularity threshold used to partition the network: ${modularity_threshold}
Coverage SD for bold mode: ${bold_sd_coverage}
Reference genome to evaluate the results of gplas: ${reference}
START

## Set up snakemake config templates
cat <<EOF >templates/"${name}"_assembly.yaml
samples:
  "${name}": "${input}"
reference:
  "${reference}": "${reference}"
species: "${species}"
threshold_prediction: "${threshold_prediction}"
number_iterations: "${number_iterations}"
classifier: "${classifier}"
filt_gplas: "${filt_gplas}"
edge_gplas: "${edge_gplas}"
name: "${name}"
modularity_threshold: "${modularity_threshold}"
bold_sd_coverage: "${bold_sd_coverage}"
EOF

sleep 1s

## Test for conda installation
if command -v conda > /dev/null; then
  echo -e "\nConda is present\n"
else
  echo -e "\nError: Conda is needed to run gplas.\n Please install conda before running gplas."
  exit 1
fi

# Activate gplas conda env
eval "$(conda shell.bash hook)"
{
  # TRY activate existing gplas env
  conda activate gplas
} || { # CATCH create gplas env with snakemake and activate it
  echo -e "Creating a conda environment to install and run snakemake"
  conda create --name gplas --file envs/spec-snakemake.txt && \
  conda activate gplas
}

## Run analysis
if [ "$classifier" == "extract" ]; then
  # If classifier is extract, then unlock folder, perform extraction mode and quit gplas
  echo -e "We need to extract the contigs first from the assembly graph, use later those contigs for your binary prediction."
  snakemake --unlock --use-conda --configfile templates/"$name"_assembly.yaml -d $PWD -s $snakeFile gplas_input/"$name"_raw_nodes.fasta
  snakemake --use-conda --configfile templates/"$name"_assembly.yaml -d $PWD -s $snakeFile gplas_input/"$name"_raw_nodes.fasta
  echo -e "Next step is to predict the extracted contigs with your desired binary classifier."
  exit 0

else
  # If analysis mode classifier is provided, then do analysis
  if [ "$classifier" == "predict" ]; then
    # If classifier is external ('predict'), then verify that the provided input prediction file is valid
    echo -e "Resuming gplas using the prediction given by the user."
    echo -e "Checking if prediction file is correctly formatted."
    check_file_output=$(Rscript scripts/check_independent_prediction_format.R "$name")

    if [ -z "$check_file_output" ]; then
      # If input file is not valid, quit
      echo -e "Please modify format on input files and re-run gplas"
      exit 1 
    fi
  fi

  # Unlock folder and perform analysis
  snakemake --unlock --configfile templates/"$name"_assembly.yaml --use-conda -d $PWD -s $snakeFile results/normal_mode/"$name"_results.tab 
  snakemake --use-conda --configfile templates/"$name"_assembly.yaml -d $PWD -s $snakeFile results/normal_mode/"$name"_results.tab

  if [ -f results/normal_mode/"$name"_bin_Unbinned.fasta ]; then
    # Unlock & Rerun in bold mode if there are unbinned plasmids
    echo -e "Some plasmid contigs were left unbinned, running gplas in bold mode"
    snakemake --unlock --use-conda --configfile templates/"$name"_assembly.yaml -d $PWD -s $snakeFile results/"$name"_results.tab
    snakemake --use-conda --configfile templates/"$name"_assembly.yaml -d $PWD -s $snakeFile results/"$name"_results.tab
  else
    #move the results
    echo -e "No contigs were left unbinned, is not necessary to run gplas in bold mode"
    mv results/normal_mode/"$name"* results/
  fi
fi

# if the keep flag is not specified, remove all intermediate files
if [ -z "$keep" ] && [ "$classifier" != "extract" ]; then
  echo -e "Intermediate files will be deleted. If you want to keep these files, use the -k flag"
  bash scripts/remove_intermediate_files.sh -n "$name"
fi

## Check that output has been correctly created
if [ -f "results/${name}_results.tab" ]; then
  # gplas was able to detect plasmids in the assembly graph
  cat figures/logo.txt
  cat <<YES_PREDICTION
Congratulations! Prediction succesfully done.
Your results are in results/ and walks/

We hope it helps your research, thanks for using gplas!
If you have used plasflow as a classifier please cite:
  Pawel S Krawczyk et al. PlasFlow: predicting plasmid sequences in metagenomic data using genome signatures, Nucleic Acids Research, doi: 10.1093/nar/gkx1321
If you have used mlplasmids as a classifier please cite:
  Arredondo-Alonso et al. mlplasmids: a user-friendly tool to predict plasmid- and chromosome-derived sequences for single species, Microbial Genomics, doi: 10.1099/mgen.0.000224

Thank you for using gplas version ${version}! - https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
YES_PREDICTION
else
  # gplas was not able to detect plasmids in the assembly graph
  cat <<NO_PREDICTION
Looks like no plasmids could be detected in your assembly graph
Please check the file:   coverage/*_clean_prediction.tab    If all contigs were predicted as chromosome, gplas probably failed at the step to create random walks starting from plasmid seeds
If that's the case, probably your isolate does not carry any plasmid(s)
If you don't see any files present at:   gplas_input/  or  coverage/   most likely the installation of gplas failed at some point

If you have used plasflow as a classifier please cite:
  Pawel S Krawczyk et al. PlasFlow: predicting plasmid sequences in metagenomic data using genome signatures, Nucleic Acids Research, doi: 10.1093/nar/gkx1321
If you have used mlplasmids as a classifier please cite:
  Arredondo-Alonso et al. mlplasmids: a user-friendly tool to predict plasmid- and chromosome-derived sequences for single species, Microbial Genomics, doi: 10.1099/mgen.0.000224

Thank you for using gplas version ${version}! - https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
NO_PREDICTION
fi

exit 0