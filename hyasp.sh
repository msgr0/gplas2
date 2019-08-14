#!/bin/bash

##to debug
#set -e
#set -v
#set -x

########################################
## A bash script to run the gplas pipeline
## This script has been converted and transformed from the script present in the gitlab repo 'bactofidia' by aschuerch

while getopts ":i:n:s:t:x:r:h" opt; do
 case $opt in
   h)
   echo "Welcome to the user guide to run gplas.sh:"
   echo -e "Basic usage: ./hyasp.sh -i mygraph.gfa"
   echo -e "Input:\n \t -i \t graph file in *.gfa format used to extract nodes and links"
   echo -e "Projectname/Output:\n \t -n \t project name given to gplas. Default: 'unnamed'"
   echo -e "Settings: \n \t -s \t bacterial species from the graph file. If bacterial species corresponds to:
                'Enterococcus faecium','Klebsiella pneumoniae' or 'Escherichia coli' then prediction will be perfomed using mlplasmids. Default: 'unknown'"
   echo -e "\t -t \t threshold to predict plasmid-derived sequences. Default: 0.5"
   echo -e "\t -x \t Number of times gplas finds plasmid paths per each plasmid seed. Default: 10"
   echo -e "Benchmarking: \n \t -r Complete reference genome corresponding to the graph given"
   exit
   ;;
   i)
     input=$OPTARG
     ;;
   n)
     name=$OPTARG
     ;;
   s)
      species=$OPTARG
       ;;
   h)
     help=$OPTARG
     ;;
   t)
     threshold_prediction=$OPTARG
     ;;
   x)
     number_iterations=$OPTARG
     ;;
   r)
     reference=$OPTARG
     ;;
   \?)
     echo "Invalid option: -$OPTARG" >&2
     echo "Ups shomething went wrong!"
     echo "hyasp requires a single input corresponding to a gfa file, use the following syntax:"
     echo "./hyasp.sh -i mygraph.gfa"
     exit 1
     ;;
   :)
     echo "Option -$OPTARG requires an argument." >&2
     echo "Ups shomething went wrong!"
     echo "hyasp requires a single input corresponding to a gfa file, use the following syntax:"
     echo "./hyasp.sh -i mygraph.gfa"
     exit 1
     ;;
 esac
done
echo -e "\n"
echo -e "!!!!!!Welcome to the script to run and evaluate hyasp!!!!!!\n"
echo -e "####### Preparation of the files for snakemake #########################\n"

if [ -z "$input" ];
then
    echo -e "Ups, it seems that you are missing the input graph. Use gplas in the following way:\n"
    echo "./hyasp.sh -i mygraph.gfa"
    exit
else
  echo -e "This is your INPUT graph:" $input "\n"
fi


if [ -z "$species" ];
then
    echo -e "You have not indicated the SPECIES that you are trying to predict. Your genome will predicted using plasflow\n"
    species="metagenome"
else
  echo -e "This is the SPECIES that you are trying to predict:" $species "\n"
  species=$(echo "'"$species"'")
fi

classifier="hyasp"

if [ -z "$name" ];
then
    echo -e "You have not passed the NAME of the project. Your results will be named as 'unnamed_project''\n"
    name="unnamed_project"
else
  echo -e "Your results will be named using" $name "\n"
fi

if [ -z "$threshold_prediction" ];
then
    echo -e "You have not passed a minimum threshold to carry out the plasmid prediction, using 0.5 as default\n"
    threshold_prediction=NA
else
  echo -e "You have indicated a threshold prediction of:" $threshold_prediction "\n"
fi

if [ -z "$number_iterations" ];
then
    echo -e "You have not passed the number of times to look for paths based on each plasmid seed, using 10 as default\n"
    number_iterations=NA
else
  echo -e "You have indicated a number of iterations of:" $number_iterations "\n"
fi

if [ -z "$reference" ];
then
    reference="No reference provided"
else
  echo -e "You have given a reference genome to evaluate the results of gplas:" $reference
  mkdir -p reference_genome
  rm reference_genome/*
  cp $reference reference_genome/
  mv reference_genome/*.fasta reference_genome/"$name"_ref_genome.fasta
fi


echo "##################################################################"

sleep 10s

( echo "cat <<EOF >final.yaml";
  cat template.yml;
  echo "EOF";
) >temp.yaml
. temp.yaml
cat final.yaml


echo -e "Let's check if snakemake is present in a previous conda environment"
source activate snakeplas

snakemake --use-conda  -s hyasp.smk hyasp/"$name"_completeness.tab
