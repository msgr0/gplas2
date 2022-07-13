#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 13:16:56 2022

@author: jpaganini
"""

"""
A python wrapped to run the gplas pipeline.
## This script has been influenced by and adapted from the code in 'bactofida', by aschuerch:
## https://gitlab.com/aschuerch/bactofidia
"""

from cProfile import run
import shutil
import glob
import os
import sys
import argparse
import time
import subprocess
from pathlib import Path
from .version import version as VERSION
# version
#VERSION="1.0.0"

# Directories
pkgdir = os.path.dirname(__file__)
snkdir = f"{pkgdir}/snakefiles"
scriptdir = f"{pkgdir}/scripts"
#******************************#
#*                            *#
#* Command line parsing       *#
#*                            *#
#******************************#

class PriorityPrinting(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if option_string == "-h" or option_string == "--help":
            parser.print_help()
        elif option_string == "-v" or option_string == "--version":
            print(f"gplas version {VERSION}")    
        parser.exit()

#create a function to pass float ranges
parser = argparse.ArgumentParser(description='gplas (A tool for binning plasmid-predicted contigs into individual predictions.',formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=False)
parser.register('action','printing',PriorityPrinting)
parser.add_argument('-i','--input',type=str, required=True, help="Path to the graph file in GFA (.gfa) format, used to extract nodes and links")
parser.add_argument('-c','--classifier',type=str, required=True, choices=['mlplasmids','extract','predict'], help="Classifier used to predict the contigs extracted from the input graph.")
parser.add_argument('-s','--species',type=str, required='mlplasmids' in sys.argv, default='Unknown', help="Bacterial species from the graph file.\
      If you have specified mlplasmids as classifier you need to provide one of the following three bacterial species:\
      'Enterococcus faecium', 'Enterococcus faecalis', 'Klebsiella pneumoniae', 'Acinetobacter baumannii' or\
      'Escherichia coli'.")
parser.add_argument('-n','--name',type=str, default='unnamed', help="Output name used in the gplas files")
parser.add_argument('-k','--keep', action='store_true', help="Keep intermediary files")
parser.add_argument('-t','--threshold_prediction',type=float,help="Prediction threshold for plasmid-derived sequences")
parser.add_argument('-b','--bold_walks',type=int, default=5, help="Coverage variance allowed for bold walks to recover unbinned plasmid-predicted nodes")
parser.add_argument('-x','--number_iterations',type=int, default=20,help="Number of walk iterations per starting node")
parser.add_argument('-f','--filt_gplas',type=float, default=0.1, help="filtering threshold to reject outgoing edges")
parser.add_argument('-e','--edge_threshold', type=float, default=0.1,help="Edge threshold")
parser.add_argument('-q','--modularity_threshold',type=float, default=0.2, help="Modularity threshold to split components in the plasmidome network")
parser.add_argument('-l','--length_filter',type=int, default=1000, help="Filtering threshold for sequence length")
parser.add_argument('-h','--help',action='printing',nargs=0,help="Prints this message")
parser.add_argument('-v','--version',action='printing',nargs=0,help="Prints gplas version")
parser.add_argument('-P','--prediction',default='independent_prediction/',help="Location of independent prediction input")
args = parser.parse_args()





#Get help function
def help_message():
    message="""
Welcome to the user guide of gplas (version 1.0.0).

BASIC USAGE
  gplas.py [-i <file>] [-c <string>] [-s <string>] [...]

Example:
  gplas.py -i mygraph.gfa -c mlplasmids -s 'Enterococcus faecium'

USER OPTIONS
Input:
  -i  [MANDATORY] Path to the graph file in GFA (.gfa) format, used to extract nodes and links. 
  -c  [MANDATORY] Classifier used to predict the contigs extracted from the input graph.
      String value: 'mlplasmids', 'extract' or 'predict'.
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
  -t  [Optional] Threshold to predict plasmid-derived sequences. Float value ranging from 0 to 1. 
      Only valid for classifier options. 'mlplasmids' (default: 0.7).
  -b  [Optional] Coverage variance allowed for bold walks to recover unbinned plasmid-predicted nodes.
      Integer value: X times coverage variance of the chromsome. Default: 5
  -x  [Optional] Number of times gplas finds plasmid walks per each plasmid starting node.
      Integer value ranging from 1 to infinite. Default: 20
  -f  [Optional] Gplas filtering threshold score to reject possible outcoming edges.
      Float value ranging from 0 to 1. Default: 0.1
  -q  [Optional] Modularity threshold to split components present in the plasmidome network.
      Float value ranging from 0 to 1. Default: 0.2
  -l  [Optional] Filter threshold for minimum length of sequences to be considered.
      Integer value ranging from 0 to infinite. Default: 1000

Other:
  -h  Print this help message.
  -v  Print gplas version.
          """
    return print(message)

## Print version
          
#Success Message

def success_message():
    print ('\n')
    print(read_logo)
    print ('\n')
    print(f"""

Congratulations! Prediction succesfully done.

Your results are in results/

We hope it helps your research, thanks for using gplas!
If you have used smlplasmids as a classifier please cite:
  Arredondo-Alonso et al. mlplasmids: a user-friendly tool to predict plasmid- and chromosome-derived sequences for single species, Microbial Genomics, doi: 10.1099/mgen.0.000224

Thank you for using gplas version {VERSION}! - https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
""")
    
    sys.exit(0)
    
#Prediction not successful
def error_message():    
    print('\n')
    print("""
Looks like no plasmids could be detected in your assembly graph

Please check the file:   coverage/*_clean_prediction.tab. 
If all contigs were predicted as chromosome, gplas probably failed at the step to create random walks starting from plasmid seeds. If that's the case, probably your isolate does not carry any plasmid(s)
If you don't see any files present at:   gplas_input/  or  coverage/  most likely the installation of gplas failed at some point

If you have used mlplasmids as a classifier please cite:
  Arredondo-Alonso et al. mlplasmids: a user-friendly tool to predict plasmid- and chromosome-derived sequences for single species, Microbial Genomics, doi: 10.1099/mgen.0.000224

Thank you for using gplas version 1.0.0! - https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
          
        """)
    sys.exit(0)

#******************************#
#*                            *#
#*  Initial exploration of    *#
#*         arguments          *#
#*                            *#
#******************************#

# Runs a (snakemake) command and checks if it returned without error
def runSnake(command):
    run = subprocess.run(command, shell=True, text=True, executable='/bin/bash')
    if run.returncode == 1:
        print("It looks like Snakemake encountered an error. Please refer to the Snakemake output to see what went wrong")
        sys.exit(1)
    else:
        return 0

if args.help==True:
    help_message()
    sys.exit(0)
    
if args.version==True:
    print_version()
    sys.exit(0)
 
if args.classifier=='mlplasmids':
   snakeFile=f'{snkdir}/mlplasmidssnake.smk'
   list_species=['Escherichia coli', 'Klebsiella pneumoniae', 'Acinetobacter baumannii','Enterococcus faecium','Enterococcus faecalis']
   if args.species not in list_species:
       print("Error: You have specified mlplasmids as classifier but you have not indicated a valid species. Please select one of the following:'Enterococcus faecium', 'Enterococcus faecalis', 'Klebsiella pneumoniae', 'Acinetobacter baumannii' or 'Escherichia coli'.")
       sys.exit(1)
  


elif args.classifier=='plasflow':
    print("PlasFlow has been deprecated for use with gplas since the tool is no longer actively maintained. If you insist on using PlasFlow, please use the 'extract' and 'predict' options of gplas")
    sys.exit(1)
    #snakeFile=f"{snkdir}/plasflowsnake.smk"
    #if args.threshold_prediction is None:
    #    args.threshold_prediction=0.7
    
elif args.classifier=='extract' or args.classifier=='predict':
    snakeFile=f"{snkdir}/otherclassifiers.smk"
    
    
    
if args.threshold_prediction is None:
    args.threshold_prediction=0.5

#******************************#
#*                            *#
#*        Start gplas         *#
#*                            *#
#******************************#

#Print messages
logo_file=open(f'{pkgdir}/figures/logo.txt', 'r')
read_logo=logo_file.read()
print ('\n')
print(read_logo)
print("##################################################################")

#Print chosen parameters
print("Your results will be named", args.name)
print("Input graph: ", args.input)
print("Bacterial species: ", args.species)
print("Classifier: ", args.classifier)
print("Threshold for predicting plasmid-derived contigs: ", args.threshold_prediction)
print("Number of plasmid walks created per node: ", args.number_iterations)
print("Threshold of gplas scores: ", args.filt_gplas)
print("Minimum frequency to consider an edge: ", args.edge_threshold)
print("Modularity threshold used to partition the network: ", args.modularity_threshold)
print("Coverage SD for bold mode: ", args.bold_walks)
print("Minimum sequence length: ", args.length_filter)
print("##################################################################")    

#Set up snakemake config files

template_file=f'templates/{args.name}_assembly.yaml'
if not os.path.exists("templates"):
    os.mkdir('templates')


with open(template_file, 'w+') as template:
    template.write('samples:'+'\n')
    template.write(f'  "{args.name}": "{args.input}"\n')
    template.write(f'species: "\'{args.species}\'"\n')
    template.write(f'threshold_prediction: "{str(args.threshold_prediction)}"\n')
    template.write(f'number_iterations: "{str(args.number_iterations)}"\n')
    template.write(f'classifier: "{args.classifier}"\n')
    template.write(f'filt_gplas: "{str(args.filt_gplas)}"\n')
    template.write(f'edge_gplas: "{str(args.edge_threshold)}"\n')
    template.write(f'name: "{args.name}"\n')
    template.write(f'modularity_threshold: "{str(args.modularity_threshold)}"\n')
    template.write(f'bold_sd_coverage: "{str(args.bold_walks)}"\n')
    template.write(f'min_node_length: "{str(args.length_filter)}"\n')
    template.write(f'predict_dir: "{str(args.prediction)}"\n')

time.sleep(1)
   
#3. Run analysis

    

if args.classifier=='predict':
        ##3.2 If classifier is external ('predict'), then verify that the provided input prediction file is the correct format
        print("Resuming gplas using the prediction given by the user.\n")
        print("Checking if prediction file is correctly formatted.\n")
        print(f"Checking for {args.prediction}")
        check_file_output_command=f'Rscript {scriptdir}/check_independent_prediction_format.R {args.name} {args.prediction}'
        check_file_run=subprocess.run(check_file_output_command, shell=True, text=True, executable='/bin/bash',capture_output=True)
        print(check_file_run.stdout)
        print(check_file_run.stderr)
        if check_file_run.returncode == 0:
            print(check_file_run.stdout)
            
        else:
            print(check_file_run.stderr)
            print("Please modify format on input files and re-run gplas")
            sys.exit(1)

if args.classifier=='extract':
    ##3.1  If classifier is extract, then unlock folder, perform extraction mode and quit gplas
    print("We need to extract the contigs first from the assembly graph, use later those contigs for your binary prediction.\n")
    command_snakemake_unlock=f'snakemake --unlock --use-conda --configfile {template_file} -d $PWD -s {snakeFile} gplas_input/{args.name}_raw_nodes.fasta'
    command_snakemake_run=f'snakemake --use-conda --configfile {template_file} -d $PWD -s {snakeFile} gplas_input/{args.name}_raw_nodes.fasta'
    runSnake(command_snakemake_unlock)
    runSnake(command_snakemake_run)
    

else:
    ##3.3 Run snakemake workflows        

    command_snakemake_unlock=f'snakemake --unlock --use-conda --configfile  {template_file} -d $PWD -s {snakeFile} results/normal_mode/{args.name}_results.tab'
    command_snakemake_run=f'snakemake --use-conda --configfile {template_file} -d $PWD -s {snakeFile} results/normal_mode/{args.name}_results.tab'
    runSnake(command_snakemake_unlock)
    runSnake(command_snakemake_run)

    ##3.4 Check if there are Unbinned contigs
    unbinned_path=f'results/normal_mode/{args.name}_bin_Unbinned.fasta'
    if os.path.exists(unbinned_path):
        ##3.5 If there are contigs left unbinned, unlock and run gplas in bold-mode
        print('\n')
        print('Some contigs were left Unbinned, running gplas in bold mode')
        print('\n')
        command_snakemake_unlock=f' snakemake --unlock --use-conda --configfile {template_file} -d $PWD -s {snakeFile} results/{args.name}_results.tab'
        command_snakemake_run=f'snakemake --use-conda --configfile {template_file} -d $PWD -s {snakeFile} results/{args.name}_results.tab'
        runSnake(command_snakemake_unlock)
        runSnake(command_snakemake_run)
    else:
        ##3.6 If there was not unbinned contigs, just move results files to the final location.
        for file in glob.glob(f"results/normal_mode/{args.name}*"):
            shutil.copy(file, "results/")
           
##3.7 If the -k flag was not selected, delete intermediary files
if args.keep==False and args.classifier!='extract':
  print("Intermediate files will be deleted. If you want to keep these files, use the -k flag")
  remove_command=f'bash {scriptdir}/remove_intermediate_files.sh -n '+args.name
  subprocess.run(remove_command, shell=True, text=True, executable='/bin/bash')
      
##3.8 Check that output has been correctly created
final_results_path='results/'+args.name+'_results.tab'
if os.path.exists(final_results_path):
    success_message()
    sys.exit(0)
    
else:
    error_message()
    sys.exit(1)


def start():
    print("Starting gplas")
