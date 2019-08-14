gplas documentation
================

# Gplas

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
```

### Installation of all the dependencies using a test graph

``` bash
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
prediction, contig name, k-mer coverage, length, component assigned.

``` r
results <- read.table(file = 'results/using_mlplasmids_results.tab')
results
```

    ##        V1              V2           V3         V4
    ## 1  number Prob_Chromosome Prob_Plasmid Prediction
    ## 2      18            0.12         0.88    Plasmid
    ## 3      33            0.35         0.65    Plasmid
    ## 4      42            0.16         0.84    Plasmid
    ## 5      47            0.48         0.52    Plasmid
    ## 6      49            0.45         0.55    Plasmid
    ## 7      50            0.08         0.92    Plasmid
    ## 8      52               0            1    Plasmid
    ## 9      54             0.1          0.9    Plasmid
    ## 10     55            0.04         0.96    Plasmid
    ## 11     56            0.29         0.71    Plasmid
    ## 12     57            0.29         0.71    Plasmid
    ## 13     60            0.37         0.63    Plasmid
    ## 14     66            0.01         0.99    Plasmid
    ##                                        V5       V6     V7        V8
    ## 1                             Contig_name coverage length Component
    ## 2  S18_LN:i:54155_dp:f:1.0514645940835776     1.05  54155         1
    ## 3  S33_LN:i:18202_dp:f:1.1628830074648842     1.16  18202         1
    ## 4   S42_LN:i:10840_dp:f:1.123936712804688     1.12  10840         1
    ## 5   S47_LN:i:8177_dp:f:0.9996798934685464        1   8177         1
    ## 6   S49_LN:i:5022_dp:f:1.1574796092139463     1.16   5022         1
    ## 7   S50_LN:i:4993_dp:f:1.1698997426343487     1.17   4993         1
    ## 8   S52_LN:i:4014_dp:f:0.9783821389091624     0.98   4014         1
    ## 9   S54_LN:i:3077_dp:f:1.1553028848000615     1.16   3077         1
    ## 10  S55_LN:i:2927_dp:f:1.1906170373500302     1.19   2927         1
    ## 11  S56_LN:i:2716_dp:f:1.1248842281377909     1.12   2716         1
    ## 12  S57_LN:i:2626_dp:f:0.9929149754371588     0.99   2626         1
    ## 13  S60_LN:i:1589_dp:f:1.0577429501871556     1.06   1589         1
    ## 14  S66_LN:i:1102_dp:f:0.8307959555606772     0.83   1102         1

### components.tab

Tab delimited file containing the bin prediction reported by gplas with
the following information: contig number, component
assignation

``` r
components <- read.table(file = 'results/using_mlplasmids_components.tab')
components
```

    ##        V1        V2
    ## 1  number Component
    ## 2      18         1
    ## 3      33         1
    ## 4      52         1
    ## 5      57         1
    ## 6      42         1
    ## 7      56         1
    ## 8      50         1
    ## 9      49         1
    ## 10     47         1
    ## 11     60         1
    ## 12     54         1
    ## 13     55         1
    ## 14     66         1

### plasmidome\_network.png

Png file of the plasmidome network generated by gplas after creating an
undirected graph using the significant co-occurrence links corresponding
to plasmid starting nodes.

![](results/using_mlplasmids_plasmidome_network.png)<!-- -->

### components.fasta

Fasta files with the nodes belonging to each predicted
    component.

``` bash
grep '>' results/unnamed_project*.fasta
```

    ## results/unnamed_project_component_1.fasta:>S18_LN:i:54155_dp:f:1.0514645940835776
    ## results/unnamed_project_component_1.fasta:>S47_LN:i:8177_dp:f:0.9996798934685464
    ## results/unnamed_project_component_1.fasta:>S52_LN:i:4014_dp:f:0.9783821389091624
    ## results/unnamed_project_component_1.fasta:>S54_LN:i:3077_dp:f:1.1553028848000615
    ## results/unnamed_project_component_1.fasta:>S57_LN:i:2626_dp:f:0.9929149754371588
    ## results/unnamed_project_component_2.fasta:>S31_LN:i:21202_dp:f:1.194722937126809
    ## results/unnamed_project_component_2.fasta:>S46_LN:i:8487_dp:f:1.2210058174026983
    ## results/unnamed_project_component_2.fasta:>S50_LN:i:4993_dp:f:1.1698997426343487
    ## results/unnamed_project_component_2.fasta:>S60_LN:i:1589_dp:f:1.0577429501871556
