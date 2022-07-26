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
    - [Input files](#input-files)
    - [Binary classification of nodes - Preprocessing](#binary-classification-of-nodes-using-an-external-tool)
        - [Using plasmidEC](#using-plasmidec)
        - [Using a different binary classifier](#using-a-different-tool)
    - [Predict plasmids](#predict-plasmids)
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

# Usage

### Input files

Gplas needs two inputs:

1) An assembly graph in **.gfa** format. Such an assembly graph can be obtained
with [SPAdes genome assembler](https://github.com/ablab/spades) or with [Unicycler](https://github.com/rrwick/Unicycler). 

2) A tab-separated file containing a binary classification (plasmid/chromsome) of each node in the assembly graph. See [below](https://gitlab.com/mmb-umcu/gplas#binary-classification-of-nodes-using-an-external-tool) for instructions on how to obtain this file. 

### Binary classification of nodes - Preprocessing <a name="binary-classification-of-nodes-using-an-external-tool"></a>

-To predict individual plasmids, some preprocessing is needed. Gplas requires that nodes in assembly graph are classified as either plasmid or chromosome, and these classifications should be summarised in a **tab-separated** file with an specific format.

-This classification step has to be completed by using an **external classification tool**. We strongly recommend using [plasmidEC](https://github.com/lisavader/plasmidEC) for this step. However, all binary classification tools are compatible with gplas.

##### <ins>Using plasmidEC</ins> <a name="using-plasmidec"></a>

PlasmidEC outperforms most available binary classification tools, and it offers two extra-advantages:
1) It uses assembly graphs in **.gfa** format as input (most tools can't). 
2) It outputs a **classification file** that is automatically compatible with gplas (Other tools will require extra processing of the output). 

Currently, plasmidEC can be used for binary classification of 8 species: *E. coli, K. pneumoniae, Acinetobacter baummannii, P. aeruginosa, S. enterica, S. aureus, E. faecalis, E. faecium*. Although, it is possible to train models for new species.

Follow the instructions on the [plasmidEC](https://github.com/lisavader/plasmidEC) repository to 
classify the nodes in your .gfa file. After obtaining your **classification file**, move to [Predict plasmids](https://gitlab.com/mmb-umcu/gplas/#predict-plasmids).

##### <ins>Using a different binary classifier</ins> <a name="using-a-different-tool"></a>

Other binary classification tools exist, and we've recently listed and reviewed several of these  [here](https://www.mdpi.com/2076-2607/9/8/1613). Although they are all compatible with gplas, extra preprocessing steps are required:

1) Use gplas to convert the nodes from the assembly graph to FASTA format (most binary classifiers only accept FASTA files as input). To do this, the **-c** flag should be set to **extract**.

``` bash
gplas -i test/test_ecoli.gfa -c extract -n 'my_isolate'
```

The output FASTA file will be located in: __gplas_input/__*my_isolate*_raw_nodes.fasta. 

2) Use this FASTA file as an input for the binary classification tool of your choice. 

3) Format the output file: 

The output from the binary classification tool has to be formatted as a tab separated file containing specific columns and headers (case sensitive). See a preloaded example below:

``` bash
head -n 4 gplas/independent_prediction/test_ecoli_plasmid_prediction.tab
```

| Prob\_Chromosome | Prob\_Plasmid |  Prediction  | Contig\_name                             | Contig\_length|
|-----------------:|--------------:|:-------------|:-----------------------------------------|--------------:|
|       1       |      0     |  Chromosome  |  S1\_LN:i:374865\_dp:f:1.0749885035087077   |      374865     |
|       1       |      0     |  Chromosome  | S10\_LN:i:198295\_dp:f:0.8919341045340952  |     198295    |
|       0       |      1     |    Plasmid   |  S20\_LN:i:91233\_dp:f:0.5815421095375989   |      91233     |


Once you've formatted the output file as above, move to [Predict plasmids](https://gitlab.com/mmb-umcu/gplas/#predict-plasmids).

### Predict plasmids <a name="predict-plasmids"></a>
After pre-processing, we are now ready to predict individual plasmids. 

Run gplas and set the **-c** flag to **predict**. Provide the paths to your assembly graph, using the **-i** flag, and to your binary classification file, with the **-P** flag. Set the name of your output with the **-n** flag. See example below: 

``` bash
gplas -c predict -i test/test_ecoli.gfa -P gplas/independent_prediction/test_ecoli_plasmid_prediction.tab -n 'my_isolate'
```
*Note: If you didn't use plasmidEC for preprocessing, make sure that the **-n** argument (in this example: **my_isolate**) matches for both the 'extract' and 'predict' commands.*

# Output files

Gplas will create a folder called ‘results’ with the following files:

``` bash
ls results/my_isolate*
```

    ## results/my_isolate_bin_1.fasta
    ## results/my_isolate_bin_2.fasta
    ## results/my_isolate_bins.tab
    ## results/my_isolate_plasmidome_network.png
    ## results/my_isolate_results.tab

##### results/\*.fasta

Fasta files with the nodes belonging to each predicted component.

``` bash
grep '>' results/my_isolate*.fasta
```

``` bash
>S32_LN:i:42460_dp:f:0.6016122804021161
>S47_LN:i:17888_dp:f:0.5893320957724726
>S50_LN:i:11225_dp:f:0.6758514700227541
>S56_LN:i:6837_dp:f:0.5759570101860518
>S59_LN:i:5519_dp:f:0.5544497698217399
>S67_LN:i:2826_dp:f:0.6746421335091037
>S20_LN:i:91233_dp:f:0.5815421095375989
```

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
| ------ | ---------------------------------------- | ---------------- | ------------- | ---------- | ------ | -------- | --- |
| 20     | S20\_LN:i:91233\_dp:f:0.5815421095375989 | 0                | 1             | Plasmid    | 91233  | 0.58     | 2   |
| 32     | S32\_LN:i:42460\_dp:f:0.6016122804021161 | 0                | 1             | Plasmid    | 42460  | 0.6      | 1   |
| 47     | S47\_LN:i:17888\_dp:f:0.5893320957724726 | 0                | 1             | Plasmid    | 17888  | 0.59     | 1   |
| 50     | S50\_LN:i:11225\_dp:f:0.6758514700227541 | 0                | 1             | Plasmid    | 11225  | 0.68     | 1   |
| 56     | S56\_LN:i:6837\_dp:f:0.5759570101860518  | 0.33             | 0.67          | Plasmid    | 6837   | 0.58     | 1   |
| 59     | S59\_LN:i:5519\_dp:f:0.5544497698217399  | 0                | 1             | Plasmid    | 5519   | 0.55     | 1   |
| 67     | S67\_LN:i:2826\_dp:f:0.6746421335091037  | 0.33             | 0.67          | Plasmid    | 2826   | 0.67     | 1   |

# Complete usage

``` bash
gplas --help
```
``` bash
usage: gplas -i INPUT -c {extract,predict} [-n NAME] [-P PREDICTION] [-k]
             [-t THRESHOLD_PREDICTION] [-b BOLD_WALKS] [-x NUMBER_ITERATIONS]
             [-f FILT_GPLAS] [-e EDGE_THRESHOLD] [-q MODULARITY_THRESHOLD]
             [-l LENGTH_FILTER] [-h] [-v]

gplas: A tool for binning plasmid-predicted contigs into individual
predictions.

optional arguments:
  -i INPUT, --input INPUT
                        Path to the graph file in GFA (.gfa) format, used to
                        extract nodes and links (default: None)
  -c {extract,predict}, --classifier {extract,predict}
                        Select to extract nodes from the assembly graph or to
                        predict individual plasmids. (default: None)
  -n NAME, --name NAME  Output name used in the gplas files (default: unnamed)
  -P PREDICTION, --prediction PREDICTION
                        Path to the binary classification file (default: None)
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
node 67-. These paths may contain inversions and rearrangements since
repeats units such as transposases which can be present several times in
the same plasmid sequence. In these cases, gplas can traverse the
sequence in different ways generating different plasmid-like paths.

``` bash
tail -n 10 walks/normal_mode/my_isolate_solutions.csv
```

``` bash
67-,70-,50-,143-
67-,70-,50-,143-
67-,70-,50-,143-
67-,70-,47+,117-,84-,59+,70-,50-,143-
67-,70-,50-,143-
67-,70-,50-,143-
67-,70-,47+,117-,84-,59+,70-,50-,143-
67-,70-,47+,117-,84-,59+,70-,50-,143-
67-,70-,50-,143-
67-,70-,50-,143-
```

For example, we can inspect in Bandage the path:
67-,70-,47+,117-,84-,59+,70-,50-,143-

![](figures/bandage_path.jpg)<!-- -->

# Issues and Bugs

You can report any issues or bugs that you find while installing/running
gplas using the [issue tracker](https://gitlab.com/sirarredondo/gplas/issues).
