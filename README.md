gplas documentation
================

# gplas

gplas is a tool to bin plasmid-predicted contigs. MORE INFO coming soon

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

1.  Conda is installed.
2.  Snakemake is present in the conda environment.
3.  Install many R libraries required to run gplas.
4.  Plasflow is present in the conda environment.
5.  mlplasmids is present in the conda environment

Running gplas without any other arguments will asume that:

  - Classifier to predict plasmid-derived contigs: plasflow
  - Probability threshold to assign plasmid-derived contigs: 0.5
  - Number of iterations to look for plasmid paths per each plasmid
    seed: 20
  - Output files will be called as ‘unnamed’

gplas will automatically create different folders containing all the
information generated. The final result is present in
’results/\*\_results.tab’

``` r
results <- read.table(file = 'results/unnamed_project_results.tab')
results
```

    ##        V1              V2           V3         V4
    ## 1  number Prob_Chromosome Prob_Plasmid Prediction
    ## 2      18            0.12         0.88    Plasmid
    ## 3      33            0.35         0.65    Plasmid
    ## 4      47            0.48         0.52    Plasmid
    ## 5      49            0.45         0.55    Plasmid
    ## 6      50            0.08         0.92    Plasmid
    ## 7      52               0            1    Plasmid
    ## 8      54             0.1          0.9    Plasmid
    ## 9      55            0.04         0.96    Plasmid
    ## 10     56            0.29         0.71    Plasmid
    ## 11     57            0.29         0.71    Plasmid
    ## 12     60            0.37         0.63    Plasmid
    ##                                        V5       V6     V7        V8
    ## 1                             Contig_name coverage length Component
    ## 2  S18_LN:i:54155_dp:f:1.0514645940835776     1.05  54155         1
    ## 3  S33_LN:i:18202_dp:f:1.1628830074648842     1.16  18202         1
    ## 4   S47_LN:i:8177_dp:f:0.9996798934685464        1   8177         1
    ## 5   S49_LN:i:5022_dp:f:1.1574796092139463     1.16   5022         2
    ## 6   S50_LN:i:4993_dp:f:1.1698997426343487     1.17   4993         1
    ## 7   S52_LN:i:4014_dp:f:0.9783821389091624     0.98   4014         1
    ## 8   S54_LN:i:3077_dp:f:1.1553028848000615     1.16   3077         1
    ## 9   S55_LN:i:2927_dp:f:1.1906170373500302     1.19   2927         2
    ## 10  S56_LN:i:2716_dp:f:1.1248842281377909     1.12   2716         1
    ## 11  S57_LN:i:2626_dp:f:0.9929149754371588     0.99   2626         1
    ## 12  S60_LN:i:1589_dp:f:1.0577429501871556     1.06   1589         1

## Help page

If you want to explore more arguments, go into the help page of gplas.

``` bash
./gplas.sh -h
```

    ## Welcome to the user guide of gplas:
    ## Basic usage: ./gplas.sh -i mygraph.gfa
    ## Input:
    ##       -i      graph file in *.gfa format used to extract nodes and links
    ## Projectname/Output:
    ##       -n      project name given to gplas. Default: 'unnamed'
    ## Settings: 
    ##       -s      bacterial species from the graph file. If bacterial species corresponds to:
    ##                 'Enterococcus faecium','Klebsiella pneumoniae' or 'Escherichia coli' then prediction will be perfomed using mlplasmids. Default: 'unknown'
    ##   -t      threshold to predict plasmid-derived sequences. Default: 0.5
    ##   -x      Number of times gplas finds plasmid paths per each plasmid seed. Default: 10

## Using mlplasmids

We can specify one of the species present in mlplasmids if we want to
predict the contigs using this algorithm. In this case, we can also be
more stringent setting the probability threshold to call a contig as
‘plasmid-derived’ to
0.7.

``` bash
./gplas.sh -i test/faecium_graph.gfa -n usingmlplasmids -s 'Enterococcus faecium' -t 0.7 -x 10
```

    ## 
    ## 
    ## !!!!!!Welcome to GPLAS!!!!!!
    ## 
    ## ####### Preparation of the files for snakemake #########################
    ## 
    ## This is your INPUT graph: test/faecium_graph.gfa 
    ## 
    ## This is the SPECIES that you are trying to predict: Enterococcus faecium 
    ## 
    ## 'Enterococcus faecium' is included as one of the species in mlplasmids
    ## 
    ## using mlplasmids as classifier
    ## 
    ## Plasmid prediction will be performed using mlplasmids
    ## 
    ## Your results will be named using usingmlplasmids 
    ## 
    ## You have indicated a threshold prediction of: 0.7 
    ## 
    ## You have indicated a number of iterations of: 10 
    ## 
    ## ##################################################################
    ## samples:
    ##   "usingmlplasmids": "test/faecium_graph.gfa"
    ## species: "'Enterococcus faecium'"
    ## threshold_prediction: "0.7"
    ## number_iterations: "10"
    ## classifier: "mlplasmids"
    ## Conda is present, so there is no need to install it. Well done!
    ## 
    ## Let's check if snakemake is present in a previous conda environment, otherwise will proceed to the installation
    ## Building DAG of jobs...
    ## Creating conda environment envs/r_packages.yaml...
    ## Downloading remote packages.
    ## Environment for envs/r_packages.yaml created (location: .snakemake/conda/be5729a6)
    ## Using shell: /bin/bash
    ## Provided cores: 1
    ## Rules claiming more threads will be scaled down.
    ## Job counts:
    ##  count   jobs
    ##  1   awk_links
    ##  1   awk_nodes
    ##  1   gplas_coocurr
    ##  1   gplas_coverage
    ##  1   gplas_paths
    ##  1   mlplasmids
    ##  6
    ## 
    ## [Wed Jul 10 13:04:46 2019]
    ## Job 5: Extracting the links from the graph test/faecium_graph.gfa
    ## 
    ## Activating conda environment: /home/sergi/gplas/.snakemake/conda/be5729a6
    ## [Wed Jul 10 13:04:47 2019]
    ## Finished job 5.
    ## 1 of 6 steps (17%) done
    ## 
    ## [Wed Jul 10 13:04:47 2019]
    ## Job 4: Extracting the nodes from the graph test/faecium_graph.gfa
    ## 
    ## Activating conda environment: /home/sergi/gplas/.snakemake/conda/be5729a6
    ## [Wed Jul 10 13:04:49 2019]
    ## Finished job 4.
    ## 2 of 6 steps (33%) done
    ## 
    ## [Wed Jul 10 13:04:49 2019]
    ## rule mlplasmids:
    ##     input: gplas_input/usingmlplasmids_raw_nodes.fasta
    ##     output: mlplasmids_prediction/usingmlplasmids_plasmid_prediction.tab
    ##     log: logs/usingmlplasmids_error_log_mlplasmids.txt, logs/usingmlplasmids_normal_log_mlplasmids.txt
    ##     jobid: 1
    ##     wildcards: sample=usingmlplasmids
    ## 
    ## Activating conda environment: /home/sergi/gplas/.snakemake/conda/be5729a6
    ## [Wed Jul 10 13:05:05 2019]
    ## Finished job 1.
    ## 3 of 6 steps (50%) done
    ## 
    ## [Wed Jul 10 13:05:05 2019]
    ## rule gplas_coverage:
    ##     input: mlplasmids_prediction/usingmlplasmids_plasmid_prediction.tab, gplas_input/usingmlplasmids_raw_links.txt, gplas_input/usingmlplasmids_raw_nodes.fasta
    ##     output: coverage/usingmlplasmids_graph_contigs.tab, coverage/usingmlplasmids_estimation.txt, coverage/usingmlplasmids_repeats_graph.tab, coverage/usingmlplasmids_clean_links.tab, coverage/usingmlplasmids_clean_prediction.tab, coverage/usingmlplasmids_initialize_nodes.tab
    ##     jobid: 3
    ##     wildcards: sample=usingmlplasmids
    ## 
    ## Activating conda environment: /home/sergi/gplas/.snakemake/conda/be5729a6
    ## WARNING: ignoring environment value of R_HOME
    ## [Wed Jul 10 13:05:14 2019]
    ## Finished job 3.
    ## 4 of 6 steps (67%) done
    ## 
    ## [Wed Jul 10 13:05:14 2019]
    ## rule gplas_paths:
    ##     input: mlplasmids_prediction/usingmlplasmids_plasmid_prediction.tab, coverage/usingmlplasmids_graph_contigs.tab, coverage/usingmlplasmids_estimation.txt, coverage/usingmlplasmids_repeats_graph.tab, coverage/usingmlplasmids_clean_links.tab, coverage/usingmlplasmids_clean_prediction.tab, coverage/usingmlplasmids_initialize_nodes.tab, gplas_input/usingmlplasmids_raw_nodes.fasta
    ##     output: paths/usingmlplasmids_connections.csv, paths/usingmlplasmids_solutions.csv
    ##     jobid: 2
    ##     wildcards: sample=usingmlplasmids
    ## 
    ## Activating conda environment: /home/sergi/gplas/.snakemake/conda/be5729a6
    ## WARNING: ignoring environment value of R_HOME
    ## [Wed Jul 10 13:05:29 2019]
    ## Finished job 2.
    ## 5 of 6 steps (83%) done
    ## 
    ## [Wed Jul 10 13:05:29 2019]
    ## rule gplas_coocurr:
    ##     input: mlplasmids_prediction/usingmlplasmids_plasmid_prediction.tab, coverage/usingmlplasmids_graph_contigs.tab, coverage/usingmlplasmids_estimation.txt, coverage/usingmlplasmids_repeats_graph.tab, coverage/usingmlplasmids_clean_links.tab, coverage/usingmlplasmids_clean_prediction.tab, coverage/usingmlplasmids_initialize_nodes.tab, paths/usingmlplasmids_solutions.csv, gplas_input/usingmlplasmids_raw_nodes.fasta
    ##     output: network/usingmlplasmids_plot_coocurrence_network.png, results/usingmlplasmids_results.tab, network/usingmlplasmids_components.tab
    ##     jobid: 0
    ##     wildcards: sample=usingmlplasmids
    ## 
    ## Activating conda environment: /home/sergi/gplas/.snakemake/conda/be5729a6
    ## WARNING: ignoring environment value of R_HOME
    ## 
      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |                                                                      |   1%
      |                                                                            
      |=                                                                     |   1%
      |                                                                            
      |=                                                                     |   2%
      |                                                                            
      |==                                                                    |   2%
      |                                                                            
      |==                                                                    |   3%
      |                                                                            
      |==                                                                    |   4%
      |                                                                            
      |===                                                                   |   4%
      |                                                                            
      |===                                                                   |   5%
      |                                                                            
      |====                                                                  |   5%
      |                                                                            
      |====                                                                  |   6%
      |                                                                            
      |=====                                                                 |   6%
      |                                                                            
      |=====                                                                 |   7%
      |                                                                            
      |=====                                                                 |   8%
      |                                                                            
      |======                                                                |   8%
      |                                                                            
      |======                                                                |   9%
      |                                                                            
      |=======                                                               |   9%
      |                                                                            
      |=======                                                               |  10%
      |                                                                            
      |=======                                                               |  11%
      |                                                                            
      |========                                                              |  11%
      |                                                                            
      |========                                                              |  12%
      |                                                                            
      |=========                                                             |  12%
      |                                                                            
      |=========                                                             |  13%
      |                                                                            
      |==========                                                            |  14%
      |                                                                            
      |==========                                                            |  15%
      |                                                                            
      |===========                                                           |  15%
      |                                                                            
      |===========                                                           |  16%
      |                                                                            
      |============                                                          |  17%
      |                                                                            
      |============                                                          |  18%
      |                                                                            
      |=============                                                         |  18%
      |                                                                            
      |=============                                                         |  19%
      |                                                                            
      |==============                                                        |  19%
      |                                                                            
      |==============                                                        |  20%
      |                                                                            
      |===============                                                       |  21%
      |                                                                            
      |===============                                                       |  22%
      |                                                                            
      |================                                                      |  22%
      |                                                                            
      |================                                                      |  23%
      |                                                                            
      |=================                                                     |  24%
      |                                                                            
      |=================                                                     |  25%
      |                                                                            
      |==================                                                    |  25%
      |                                                                            
      |==================                                                    |  26%
      |                                                                            
      |===================                                                   |  27%
      |                                                                            
      |====================                                                  |  28%
      |                                                                            
      |====================                                                  |  29%
      |                                                                            
      |=====================                                                 |  30%
      |                                                                            
      |======================                                                |  31%
      |                                                                            
      |======================                                                |  32%
      |                                                                            
      |=======================                                               |  32%
      |                                                                            
      |=======================                                               |  33%
      |                                                                            
      |========================                                              |  34%
      |                                                                            
      |========================                                              |  35%
      |                                                                            
      |=========================                                             |  35%
      |                                                                            
      |=========================                                             |  36%
      |                                                                            
      |==========================                                            |  37%
      |                                                                            
      |===========================                                           |  38%
      |                                                                            
      |===========================                                           |  39%
      |                                                                            
      |============================                                          |  40%
      |                                                                            
      |=============================                                         |  41%
      |                                                                            
      |=============================                                         |  42%
      |                                                                            
      |==============================                                        |  43%
      |                                                                            
      |===============================                                       |  44%
      |                                                                            
      |================================                                      |  45%
      |                                                                            
      |================================                                      |  46%
      |                                                                            
      |=================================                                     |  47%
      |                                                                            
      |==================================                                    |  48%
      |                                                                            
      |==================================                                    |  49%
      |                                                                            
      |===================================                                   |  50%
      |                                                                            
      |====================================                                  |  51%
      |                                                                            
      |====================================                                  |  52%
      |                                                                            
      |=====================================                                 |  53%
      |                                                                            
      |======================================                                |  54%
      |                                                                            
      |======================================                                |  55%
      |                                                                            
      |=======================================                               |  56%
      |                                                                            
      |========================================                              |  57%
      |                                                                            
      |=========================================                             |  58%
      |                                                                            
      |=========================================                             |  59%
      |                                                                            
      |==========================================                            |  60%
      |                                                                            
      |===========================================                           |  61%
      |                                                                            
      |===========================================                           |  62%
      |                                                                            
      |============================================                          |  63%
      |                                                                            
      |=============================================                         |  64%
      |                                                                            
      |=============================================                         |  65%
      |                                                                            
      |==============================================                        |  65%
      |                                                                            
      |==============================================                        |  66%
      |                                                                            
      |===============================================                       |  67%
      |                                                                            
      |===============================================                       |  68%
      |                                                                            
      |================================================                      |  68%
      |                                                                            
      |================================================                      |  69%
      |                                                                            
      |=================================================                     |  70%
      |                                                                            
      |==================================================                    |  71%
      |                                                                            
      |==================================================                    |  72%
      |                                                                            
      |===================================================                   |  73%
      |                                                                            
      |====================================================                  |  74%
      |                                                                            
      |====================================================                  |  75%
      |                                                                            
      |=====================================================                 |  75%
      |                                                                            
      |=====================================================                 |  76%
      |                                                                            
      |======================================================                |  77%
      |                                                                            
      |======================================================                |  78%
      |                                                                            
      |=======================================================               |  78%
      |                                                                            
      |=======================================================               |  79%
      |                                                                            
      |========================================================              |  80%
      |                                                                            
      |========================================================              |  81%
      |                                                                            
      |=========================================================             |  81%
      |                                                                            
      |=========================================================             |  82%
      |                                                                            
      |==========================================================            |  82%
      |                                                                            
      |==========================================================            |  83%
      |                                                                            
      |===========================================================           |  84%
      |                                                                            
      |===========================================================           |  85%
      |                                                                            
      |============================================================          |  85%
      |                                                                            
      |============================================================          |  86%
      |                                                                            
      |=============================================================         |  87%
      |                                                                            
      |=============================================================         |  88%
      |                                                                            
      |==============================================================        |  88%
      |                                                                            
      |==============================================================        |  89%
      |                                                                            
      |===============================================================       |  89%
      |                                                                            
      |===============================================================       |  90%
      |                                                                            
      |===============================================================       |  91%
      |                                                                            
      |================================================================      |  91%
      |                                                                            
      |================================================================      |  92%
      |                                                                            
      |=================================================================     |  92%
      |                                                                            
      |=================================================================     |  93%
      |                                                                            
      |=================================================================     |  94%
      |                                                                            
      |==================================================================    |  94%
      |                                                                            
      |==================================================================    |  95%
      |                                                                            
      |===================================================================   |  95%
      |                                                                            
      |===================================================================   |  96%
      |                                                                            
      |====================================================================  |  96%
      |                                                                            
      |====================================================================  |  97%
      |                                                                            
      |====================================================================  |  98%
      |                                                                            
      |===================================================================== |  98%
      |                                                                            
      |===================================================================== |  99%
      |                                                                            
      |======================================================================|  99%
      |                                                                            
      |======================================================================| 100%
    ## [Wed Jul 10 13:05:35 2019]
    ## Finished job 0.
    ## 6 of 6 steps (100%) done
    ## Complete log: /home/sergi/gplas/.snakemake/log/2019-07-10T130219.341147.snakemake.log
    ## 
    ## 
    ## Congratulations! We have finished!!!!
    ## 
    ## Your settings were the following:
    ## 
    ## samples:
    ##   "usingmlplasmids": "test/faecium_graph.gfa"
    ## species: "'Enterococcus faecium'"
    ## threshold_prediction: "0.7"
    ## number_iterations: "10"
    ## classifier: "mlplasmids"
    ## 
    ## 
    ## Thank you for using gplas, hasta la vista amig@ :)
