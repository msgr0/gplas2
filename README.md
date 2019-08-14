-   [Installation](#installation)
    -   [Installation of all the dependencies using a test
        graph](#installation-of-all-the-dependencies-using-a-test-graph)
-   [Usage](#usage)
    -   [Quick usage](#quick-usage)
-   [Help page](#help-page)
-   [Using mlplasmids](#using-mlplasmids)

\#gplas

gplas is a tool to bin plasmid-predicted contigs using co-occurence
networks. Gplas is a new implementation of mlplasmids which adds the
possibility of accurately binning predicted plasmid contigs into several
independent components. We rely on the structure of the assembly graph
and k-mer information to create a new plasmidome network and predict
plasmid contigs belonging to different genomic units.

Installation
------------

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

1.  Conda is installed.
2.  Snakemake is present in the conda environment.
3.  Install all the R libraries required to run gplas.

      igraph version 1.2.4.1

      ggraph version 1.0.2

      Biostrings version 2.50.2

      seqinr version 3.4-5

      tidyverse version 1.2.1

      spatstat version 1.59-0

      cooccur version 1.3

      ggrepel version 0.8.0

1.  mlplasmids (version 1.0.0) is present in the conda environment  
2.  Plasflow (version 1.1) is present in the conda environment.

Usage
-----

### Quick usage

*Running gplas providing uniquely the assembly graph*

``` bash
./gplas.sh -i test/faecium_graph.gfa
```

Gplas only requires a single argument ‘-i’ corresponding to an assembly
graph in gfa format.

Running gplas without any other arguments will asume that:

-   Classifier to predict plasmid-derived contigs: plasflow
-   Probability threshold to assign plasmid-derived contigs: 0.5
-   Number of iterations to look for plasmid paths per each plasmid
    seed: 10
-   Output files will be called as ‘unnamed\_project’

gplas will automatically create different folders containing all the
information generated. The final result is present in
’results/\*\_results.tab’

``` bash
ls results/unnamed_project*
```

    ## results/unnamed_project_component_1.fasta
    ## results/unnamed_project_component_2.fasta
    ## results/unnamed_project_components.tab
    ## results/unnamed_project_plasmidome_network.png
    ## results/unnamed_project_results.tab

``` r
results <- read.table(file = 'results/unnamed_project_results.tab')
results
```

    ##        V1              V2           V3         V4
    ## 1  number Prob_Chromosome Prob_Plasmid Prediction
    ## 2      18            0.01         0.99    Plasmid
    ## 3      31            0.15         0.85    Plasmid
    ## 4      46            0.03         0.97    Plasmid
    ## 5      47            0.04         0.96    Plasmid
    ## 6      50            0.02         0.98    Plasmid
    ## 7      52            0.03         0.97    Plasmid
    ## 8      54            0.08         0.92    Plasmid
    ## 9      57            0.03         0.97    Plasmid
    ## 10     60               0            1    Plasmid
    ##                                        V5       V6     V7        V8
    ## 1                             Contig_name coverage length Component
    ## 2  S18_LN:i:54155_dp:f:1.0514645940835776     1.05  54155         1
    ## 3   S31_LN:i:21202_dp:f:1.194722937126809     1.19  21202         2
    ## 4   S46_LN:i:8487_dp:f:1.2210058174026983     1.22   8487         2
    ## 5   S47_LN:i:8177_dp:f:0.9996798934685464        1   8177         1
    ## 6   S50_LN:i:4993_dp:f:1.1698997426343487     1.17   4993         2
    ## 7   S52_LN:i:4014_dp:f:0.9783821389091624     0.98   4014         1
    ## 8   S54_LN:i:3077_dp:f:1.1553028848000615     1.16   3077         1
    ## 9   S57_LN:i:2626_dp:f:0.9929149754371588     0.99   2626         1
    ## 10  S60_LN:i:1589_dp:f:1.0577429501871556     1.06   1589         2

Help page
---------

If you want to explore more arguments, go into the help page of gplas.

``` bash
./gplas.sh -h
```

    ## Welcome to the user guide of gplas:
    ## 
    ## Basic usage: ./gplas.sh -i mygraph.gfa
    ## 
    ## Input:
    ##       -i      Mandatory: Path to the graph file in *.gfa format used to extract nodes and links. Gfa file format
    ## Projectname/Output:
    ##       -n      Optional: Project name given to gplas. Default: 'unnamed'
    ## Settings: 
    ##       -s      Optional: Bacterial species from the graph file. If bacterial species corresponds to:
    ##                 'Enterococcus faecium','Klebsiella pneumoniae' or 'Escherichia coli' then prediction will be perfomed using mlplasmids.
    ##                  Default: 'unknown'
    ##   -t      Optional: Threshold to predict plasmid-derived sequences. Integer value ranging from 0 to 1. Default: 0.5
    ##   -x      Optional: Number of times gplas finds plasmid paths per each plasmid starting node. Integer value ranging from 1 to infinite.
    ##                  Default: 10
    ##   -m      Optional: Mode to run gplas: 'normal' or 'bold'. Bold mode increases the acceptance of connections to enlogate the path.
    ##                  String value. Default: 'normal'
    ## Benchmarking purposes: 
    ##       -r      Optional: Path to the complete reference genome corresponding to the graph given. Fasta file format

Using mlplasmids
----------------

We can specify one of the species present in mlplasmids if we want to
predict the contigs using this algorithm. In this case, we can also be
more stringent setting the probability threshold to call a contig as
‘plasmid-derived’ to 0.7.

``` bash
./gplas.sh -i test/faecium_graph.gfa -n usingmlplasmids -s 'Enterococcus faecium' -t 0.7 -x 10
```
