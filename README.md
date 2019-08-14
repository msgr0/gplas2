gplas: binning plasmid-predicted contigs
================

gplas is a tool to bin plasmid-predicted contigs using co-occurence
networks. Gplas is a new implementation of mlplasmids which adds the
possibility of accurately binning predicted plasmid contigs into several
independent components. We rely on the structure of the assembly graph
and k-mer information to create a new plasmidome network and predict
plasmid contigs belonging to different genomic units.

![](figures/logo.png)<!-- -->

## Installation

``` bash

git clone https://gitlab.com/sirarredondo/gplas.git
cd gplas
./gplas.sh -i test/faecium_graph.gfa
```

First-time installation can take some time (you don’t have to do
anything at all, so go for a coffee and come back in 10m). gplas needs
to check whether:

1.  **Conda** is installed.

2.  **Snakemake** is present in the conda environment.

3.  Install all the **R libraries** required to run gplas.

      **igraph** version 1.2.4.1

      **ggraph** version 1.0.2

      **Biostrings** version 2.50.2

      **seqinr** version 3.4-5

      **tidyverse** version 1.2.1

      **spatstat** version 1.59-0

      **cooccur** version 1.3

      **ggrepel** version 0.8.0

4.  **mlplasmids** (version 1.0.0) is present in the conda environment

5.  **Plasflow** (version 1.1) is present in the conda environment.

## Usage

### Quick usage

**Running gplas providing uniquely the assembly graph**

Gplas only requires a single argument ‘-i’ corresponding to an assembly
graph in gfa format.

``` bash
./gplas.sh -i test/faecium_graph.gfa
```

### Complete usage

Gplas can take the following arguments:

Mandatory arguments:

  - **-i**: Path to the graph file in \*.gfa format used to extract
    nodes and links. Gfa file format

Optional arguments:

  - **-n**: Project name given to gplas. Default: ‘unnamed’
  - **-s**: Bacterial species from the graph file. If bacterial species
    corresponds to: ‘Enterococcus faecium’,‘Klebsiella pneumoniae’ or
    ‘Escherichia coli’ then prediction will be perfomed using
    mlplasmids. Default: ‘unknown’
  - **-t**: Threshold to predict plasmid-derived sequences. Integer
    value ranging from 0 to 1. Default: 0.5
  - **-x**: Number of times gplas finds plasmid paths per each plasmid
    starting node. Integer value ranging from 1 to infinite. Default: 10
  - **-m**: Mode to run gplas: ‘normal’ or ‘bold’. Bold mode increases
    the acceptance of connections to enlogate the path. String value.
    Default: ‘normal’

For benchmarking purposes you can pass a complete genome to gplas and
will generate a precision and completeness. Using this you can assess
the performance of gplas on a small set of genomes in which perhaps you
have generated long-reads.

  - **r**: Path to the complete reference genome corresponding to the
    graph given. Fasta file format

For example, if we want to use mlplasmids with the *Enterococcus
faecium* model on the same test set we only need to specify the
following:

``` bash
./gplas.sh -i test/faecium_graph.gfa -n 'using_mlplasmids' -s 'Enterococcus_faecium'
```

## Main output files

Gplas will create a folder called ‘results’ with the following files:

``` bash
ls results/using_mlplasmids*
```

    ## results/using_mlplasmids_component_1.fasta
    ## results/using_mlplasmids_components.tab
    ## results/using_mlplasmids_plasmidome_network.png
    ## results/using_mlplasmids_results.tab

### results.tab

Tab delimited file containing the prediction given by mlplasmids or
plasflow together with the bin prediction assigned by gplas. The file
contains the following information: contig number, probability of being
chromosome-derived, probability of being plasmid-derived, class
prediction, contig name, k-mer coverage, length, component
assigned.

| number | Prob\_Chromosome | Prob\_Plasmid | Prediction | Contig\_name                             | coverage | length | Component |
| -----: | ---------------: | ------------: | :--------- | :--------------------------------------- | -------: | -----: | --------: |
|     18 |             0.12 |          0.88 | Plasmid    | S18\_LN:i:54155\_dp:f:1.0514645940835776 |     1.05 |  54155 |         1 |
|     33 |             0.35 |          0.65 | Plasmid    | S33\_LN:i:18202\_dp:f:1.1628830074648842 |     1.16 |  18202 |         1 |
|     42 |             0.16 |          0.84 | Plasmid    | S42\_LN:i:10840\_dp:f:1.123936712804688  |     1.12 |  10840 |         1 |
|     47 |             0.48 |          0.52 | Plasmid    | S47\_LN:i:8177\_dp:f:0.9996798934685464  |     1.00 |   8177 |         1 |
|     49 |             0.45 |          0.55 | Plasmid    | S49\_LN:i:5022\_dp:f:1.1574796092139463  |     1.16 |   5022 |         1 |
|     50 |             0.08 |          0.92 | Plasmid    | S50\_LN:i:4993\_dp:f:1.1698997426343487  |     1.17 |   4993 |         1 |
|     52 |             0.00 |          1.00 | Plasmid    | S52\_LN:i:4014\_dp:f:0.9783821389091624  |     0.98 |   4014 |         1 |
|     54 |             0.10 |          0.90 | Plasmid    | S54\_LN:i:3077\_dp:f:1.1553028848000615  |     1.16 |   3077 |         1 |
|     55 |             0.04 |          0.96 | Plasmid    | S55\_LN:i:2927\_dp:f:1.1906170373500302  |     1.19 |   2927 |         1 |
|     56 |             0.29 |          0.71 | Plasmid    | S56\_LN:i:2716\_dp:f:1.1248842281377909  |     1.12 |   2716 |         1 |
|     57 |             0.29 |          0.71 | Plasmid    | S57\_LN:i:2626\_dp:f:0.9929149754371588  |     0.99 |   2626 |         1 |
|     60 |             0.37 |          0.63 | Plasmid    | S60\_LN:i:1589\_dp:f:1.0577429501871556  |     1.06 |   1589 |         1 |
|     66 |             0.01 |          0.99 | Plasmid    | S66\_LN:i:1102\_dp:f:0.8307959555606772  |     0.83 |   1102 |         1 |

### components.tab

Tab delimited file containing the bin prediction reported by gplas with
the following information: contig number, component assignation

| number | Component |
| -----: | --------: |
|     18 |         1 |
|     33 |         1 |
|     52 |         1 |
|     57 |         1 |
|     42 |         1 |
|     56 |         1 |
|     50 |         1 |
|     49 |         1 |
|     47 |         1 |
|     60 |         1 |
|     54 |         1 |
|     55 |         1 |
|     66 |         1 |

### plasmidome\_network.png

Png file of the plasmidome network generated by gplas after creating an
undirected graph using the significant co-occurrence links corresponding
to plasmid starting nodes.

![](results/using_mlplasmids_plasmidome_network.png)<!-- -->

### components.fasta

Fasta files with the nodes belonging to each predicted component.

``` bash
grep '>' results/using_mlplasmids*.fasta
```

    ## >S18_LN:i:54155_dp:f:1.0514645940835776
    ## >S33_LN:i:18202_dp:f:1.1628830074648842
    ## >S42_LN:i:10840_dp:f:1.123936712804688
    ## >S47_LN:i:8177_dp:f:0.9996798934685464
    ## >S49_LN:i:5022_dp:f:1.1574796092139463
    ## >S50_LN:i:4993_dp:f:1.1698997426343487
    ## >S52_LN:i:4014_dp:f:0.9783821389091624
    ## >S54_LN:i:3077_dp:f:1.1553028848000615
    ## >S55_LN:i:2927_dp:f:1.1906170373500302
    ## >S56_LN:i:2716_dp:f:1.1248842281377909
    ## >S57_LN:i:2626_dp:f:0.9929149754371588
    ## >S60_LN:i:1589_dp:f:1.0577429501871556
    ## >S66_LN:i:1102_dp:f:0.8307959555606772
