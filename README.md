gplas: binning plasmid-predicted contigs
================

<div align="center"><img src="figures/logo.png" alt="gplas" width="600"/></div>

gplas is a tool to bin plasmid-predicted contigs based on sequence
composition, coverage and assembly graph information. Gplas is a new
tool that extends the possibility of accurately binning predicted
plasmid contigs into several discrete plasmid components.

# Table of Contents
- [gplas: binning plasmid-predicted contigs](#gplas-binning-plasmid-predicted-contigs)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  - [Installation using conda (to be implemented)](#installation-using-conda-to-be-implemented)
  - [Installation using pip and conda](#installation-using-pip-and-conda)
- [Usage](#usage)
    - [Step 1 - Binary classification of nodes](#binary-classification-of-nodes-using-an-external-tool)
        - [Using plasmidEC](#using-plasmidec)
        - [Using a different tool](#using-a-different-tool)
    - [ Step 2 - Predict plasmids](#predict-plasmids)
- [Output files](#main-output-files)
- [Complete usage](#complete-usage)
    - [Intermediary results files](#intermediary-results-files)
- [Issues and Bugs](#issues-and-bugs)

# Installation

## Installation using conda (to be implemented)

## Installation using pip and conda

While the conda recipe is under construction, the prefered way of installing gplas is using pip and a conda environment. Please follow the instructions below:

Clone the repository and enter the directory
``` bash
git clone https://gitlab.com/mmb-umcu/gplas.git
cd gplas
```

Create a new conda environment and activate it
``` bash
conda env create --name gplas --file envs/gplas.yaml
conda activate gplas
```

Install gplas using pip
``` bash
pip install -e .
```
When this has finished, test the installation using 
``` bash
gplas --help
```
This should should show the help page of gplas.

# Input

Gplas needs two inputs:

1) An assembly graph in **.gfa** format. Such an assembly graph can be obtained
with [SPAdes genome assembler](https://github.com/ablab/spades) or with [Unicycler](https://github.com/rrwick/Unicycler). 

2) A **.tab** file containing a binary classification (plasmid/chromsome) of each node in the assembly graph. See below for instructions on how to generate such a file.

# Usage

### Step 1 - Binary classification of nodes <a name="binary-classification-of-nodes-using-an-external-tool"></a>

To predict individual plasmids, gplas requires that nodes in assembly graph are classified as either plasmid or chromosome. This step has to be completed by using an external classification tool.

We strongly recommend using [plasmidEC](https://github.com/lisavader/plasmidEC) for this step. 

##### <ins>Using plasmidEC</ins> <a name="using-plasmidec"></a>

PlasmidEC outperforms most available tools, and it offers two extra-advantages:
1) It uses an assembly graph in **.gfa** format as input.
2) It outputs a **.tab** classification file that is automatically compatible with gplas. 

Currently, plasmidEC can be used for binary classification of 8 species: *E. coli, K. pneumoniae, Acinetobacter baummannii, P. aeruginosa, S. enterica, S. aureus, E. faecalis, E. faecium*

Follow the instructions on the [plasmidEC](https://github.com/lisavader/plasmidEC) repository to 
classify your nodes and move to step 2.

##### <ins>Using a different tool</ins> <a name="using-a-different-tool"></a>

Other binary classification tools exist, and we've recently listed and reviewed several of these  [here](https://www.mdpi.com/2076-2607/9/8/1613). Although they are all compatible with gplas, extra steps are required to complete the plasmid predictions. 

1) Use gplas to convert the nodes from the assembly graph to FASTA format. For this, the **-c** flag should be set to **extract**.

``` bash
gplas -i test/test_ecoli.gfa -c extract -n 'my_isolate'
```

The output FASTA file, containing the nodes sequences, will be located in: __gplas_input/__*my_isolate*_raw_nodes.fasta. 

2) Use this FASTA file as an input for the binary classification tool of your choice. 

3) Format the output as indicated below:

The output from the selected binary classification tools has to be formatted as a tab separated file containing the following columns and headers (case sensitive):

| Prob\_Chromosome | Prob\_Plasmid |  Prediction  | Contig\_name                             | Contig\_length|
|-----------------:|--------------:|:-------------|:-----------------------------------------|--------------:|
|       0.40       |      0.60     |    Plasmid   |  S1\_LN:i:4240\_dp:f:1.936810327946946   |      4240     |
|       0.65       |      0.35     |  Chromosome  | S18\_LN:i:147394\_dp:f:1.05847808445255  |     147394    |
|       0.12       |      0.88     |  Chromosome  |  S25\_LN:i:7135\_dp:f:2.03512069877433   |      7135     |

Once formatted, save this file with the name: **my_isolate_plasmid_prediction.tab**. 
The prefix of the file-name (in this example: **my_isolate**) must match with the argument passed to **-n** in Step 1. 

Once you've formatted the output file as above, move to Step 2.

### Step 2 - Predict plasmids <a name="predict-plasmids"></a>
Gplas will now predict individual plasmids in your sample. For this, you will run gplas setting **-c** flag to **predict**. Also, with the **-P** flag, you will indicate the path to the directory holding the binary classification file (obtained in Step 1). 

``` bash
gplas -i test/test_ecoli.gfa -c predict -n 'my_isolate' -P ${binary_classifcation_directory}
```

# Output files

Gplas will create a folder called ‘results’ with the following files:

``` bash
ls results/my_isolate*
```

    ## results/my_isolate_bin_1.fasta
    ## results/my_isolate_bins.tab
    ## results/my_isolate_plasmidome_network.png
    ## results/my_isolate_results.tab

##### results/\*.fasta

Fasta files with the nodes belonging to each predicted component.

``` bash
grep '>' results/my_isolate*.fasta
```

    ## >S18_LN:i:54155_dp:f:1.0514645940835776
    ## >S31_LN:i:21202_dp:f:1.194722937126809
    ## >S33_LN:i:18202_dp:f:1.1628830074648842
    ## >S46_LN:i:8487_dp:f:1.2210058174026983
    ## >S47_LN:i:8177_dp:f:0.9996798934685464
    ## >S50_LN:i:4993_dp:f:1.1698997426343487
    ## >S52_LN:i:4014_dp:f:0.9783821389091624
    ## >S54_LN:i:3077_dp:f:1.1553028848000615
    ## >S57_LN:i:2626_dp:f:0.9929149754371588
    ## >S60_LN:i:1589_dp:f:1.0577429501871556

##### results/\*plasmidome\_network.png

Png file of the plasmidome network generated by gplas after creating an
undirected graph from edges between plasmid unitigs co-existing in the
walks created by gplas.

![](figures/my_isolate_plasmidome_network.png)<!-- -->

##### results/\*results.tab

Tab delimited file containing the prediction given by plasmidEC (or
other binary classification tool) together with the bin prediction by gplas. The file contains
the following information: contig number, probability of being
chromosome-derived, probability of being plasmid-derived, class
prediction, contig name, k-mer coverage, length, bin assigned.

| number | Contig\_name                             | Prob\_Chromosome | Prob\_Plasmid | Prediction | length | coverage | Bin |
|-------:|:-----------------------------------------|-----------------:|--------------:|:-----------|-------:|---------:|----:|
|     18 | S18\_LN:i:54155\_dp:f:1.0514645940835776 |             0.01 |          0.99 | Plasmid    |  54155 |     1.05 |   1 |
|     31 | S31\_LN:i:21202\_dp:f:1.194722937126809  |             0.15 |          0.85 | Plasmid    |  21202 |     1.19 |   1 |
|     33 | S33\_LN:i:18202\_dp:f:1.1628830074648842 |             0.40 |          0.60 | Plasmid    |  18202 |     1.16 |   1 |
|     46 | S46\_LN:i:8487\_dp:f:1.2210058174026983  |             0.03 |          0.97 | Plasmid    |   8487 |     1.22 |   1 |
|     47 | S47\_LN:i:8177\_dp:f:0.9996798934685464  |             0.04 |          0.96 | Plasmid    |   8177 |     1.00 |   1 |
|     50 | S50\_LN:i:4993\_dp:f:1.1698997426343487  |             0.02 |          0.98 | Plasmid    |   4993 |     1.17 |   1 |
|     52 | S52\_LN:i:4014\_dp:f:0.9783821389091624  |             0.03 |          0.97 | Plasmid    |   4014 |     0.98 |   1 |
|     54 | S54\_LN:i:3077\_dp:f:1.1553028848000615  |             0.08 |          0.92 | Plasmid    |   3077 |     1.16 |   1 |
|     57 | S57\_LN:i:2626\_dp:f:0.9929149754371588  |             0.03 |          0.97 | Plasmid    |   2626 |     0.99 |   1 |
|     60 | S60\_LN:i:1589\_dp:f:1.0577429501871556  |             0.00 |          1.00 | Plasmid    |   1589 |     1.06 |   1 |


# Complete usage

``` bash
gplas -h
```
``` bash
usage: gplas -i INPUT -c {mlplasmids,extract,predict} [-s SPECIES] [-n NAME]
             [-k] [-t THRESHOLD_PREDICTION] [-b BOLD_WALKS]
             [-x NUMBER_ITERATIONS] [-f FILT_GPLAS] [-e EDGE_THRESHOLD]
             [-q MODULARITY_THRESHOLD] [-l LENGTH_FILTER] [-h] [-v]
             [-P PREDICTION]

gplas (A tool for binning plasmid-predicted contigs into individual
predictions.

optional arguments:
  -i INPUT, --input INPUT
                        Path to the graph file in GFA (.gfa) format, used to
                        extract nodes and links (default: None)
  -c {extract,predict}, --classifier {extract,predict}
                        Classifier used to predict the contigs extracted from
                        the input graph. (default: None)
  -n NAME, --name NAME  Output name used in the gplas files (default: unnamed)
  -k, --keep            Keep intermediary files (default: False)
  -t THRESHOLD_PREDICTION, --threshold_prediction THRESHOLD_PREDICTION
                        Prediction threshold for plasmid-derived sequences
                        (default: None)
  -b BOLD_WALKS, --bold_walks BOLD_WALKS
                        Coverage variance allowed for bold walks to recover
                        unbinned plasmid-predicted nodes (default: 5)
  -x NUMBER_ITERATIONS, --number_iterations NUMBER_ITERATIONS
                        Number of walk iterations per starting node (default:
                        20)
  -f FILT_GPLAS, --filt_gplas FILT_GPLAS
                        filtering threshold to reject outgoing edges (default:
                        0.1)
  -e EDGE_THRESHOLD, --edge_threshold EDGE_THRESHOLD
                        Edge threshold (default: 0.1)
  -q MODULARITY_THRESHOLD, --modularity_threshold MODULARITY_THRESHOLD
                        Modularity threshold to split components in the
                        plasmidome network (default: 0.2)
  -l LENGTH_FILTER, --length_filter LENGTH_FILTER
                        Filtering threshold for sequence length (default:
                        1000)
  -h, --help            Prints this message (default: None)
  -v, --version         Prints gplas version (default: None)
  -P PREDICTION, --prediction PREDICTION
                        Location of independent prediction input (default:
                        independent_prediction/)
```

### Intermediary results files

If the **-k** flag is selected, gplas will also **keep** all intermediary files needed to construct the plasmid predictions. For example:

##### walks/normal_mode/\*solutions.csv

gplas generates plasmid-like walks per each plasmid starting node. These
paths are used later to generate the edges from the plasmidome network
but they can also be useful to observe all the different walks starting
from a single node (plasmid unitig). These walk can be directly given to
Bandage to visualize and manually inspect a walk.

In this case, we find different possible plasmid walks starting from the
node 18+. These paths may contain inversions and rearrangements since
repeats units such as transposases which can be present several times in
the same plasmid sequence. In these cases, gplas can traverse the
sequence in different ways generating different plasmid-like paths.

``` bash
head -n 10 walks/normal_mode/my_isolate_solutions.csv
```

    ## 18+,76-,102+,33+,76-,102+,92+,47+,115-,64+,31-,79+,60-,70-,50+,64-,116+,61-,88-,89+,69-,96-,119+,64-,116+,61-,88-,90+,69-,100+,119+,64-,116+,63+,115-,64+,119-,100-,69+
    ## 18+,76-,52+,94+,57-,77+,18+
    ## 18+,76-,52+,94+,57-,77+,87-,65+,54-,94+
    ## 18+,76-,52+,94+,57-,77+,87-,65+,54-,94+
    ## 18+,76-,52+,94+,57-,77+,87-,65+,54-,94+
    ## 18+,76-,102+,33+,76-,52+,94+,57-,77+,18+
    ## 18+,76-,52+,94+,57-,77+,18+
    ## 18+,76-,102+,92+,47+,115-,64+,31-,79+,60-,70-,50+,64-,113+
    ## 18+,76-,52+,94+,57-,77+,18+
    ## 18+,76-,102+,92+,47+,115-,64+,31-,79+,46-,79+,60-,70-,50+,64-,114+

For example, we can inspect in Bandage the path:
18+,76-,52+,94+,57-,77+,18+

This path forms a circular sequence since there is overlap between the
initial and end node of the path.

![](figures/bandage_path.jpg)<!-- -->

# Issues and Bugs

You can report any issues or bugs that you find while installing/running
gplas using the [issue tracker](https://gitlab.com/sirarredondo/gplas/issues).
