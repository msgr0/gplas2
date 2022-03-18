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

import os
import sys
import argparse


parser = argparse.ArgumentParser(description='gplas (A tool for binning plasmid-predicted contigs into individual predictions.')
parser.add_argument('-i','--input',type=str, required=True)
parser.add_argument('-c','--classifier',type=str, required=True)
parser.add_argument('-s','--species',type=str)
parser.add_argument('-n','--name',type=int, default='unnamed')
parser.add_argument('-k','--keep', action='store_true')
parser.add_argument('-t','--threshold_probability',type=float, default=0.5)
parser.add_argument('-b','--bold_walks',type=int, default=5)
parser.add_argument('-x','--walks',type=int, default=20)
parser.add_argument('-f','--filtering_threshold',type=float, default=0.1)
parser.add_argument('-q','--modularity_threshold',type=float, default=0.2)
parser.add_argument('-l','--length_filter',type=int, default=1000)
parser.add_argument('-h','--help',action='store_true')
parser.add_argument('-v','--version',action='store_true')

args = parser.parse_args()



version="1.0.0"



#help function

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
  -t  [Optional] Threshold to predict plasmid-derived sequences. Float value ranging from 0 to 1. 
      Only valid for classifier options 'plasflow' (default: 0.5) and 'mlplasmids' (default: 0.7).
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
          
def print_version():
    print("gplas version", version)
    
