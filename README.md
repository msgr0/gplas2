gplas
=====

gplas is a tool to bin plasmid-predicted contigs. MORE INFO coming soon

Installation
------------

    git clone https://gitlab.com/sirarredondo/gplas.git
    cd gplas

### Installation of all the dependencies using a test graph

    ./gplas.sh -i test/faecium_graph.gfa

Help page
---------

    ./gplas.sh -h

    ## Welcome to the user guide of gplas:
    ## Basic usage: ./gplas.sh -i mygraph.gfa
    ## Input:
    ##       -i      graph file in *.gfa format used to extract nodes and links
    ## Output:
    ##       -o      folder name in which the results will be stored. Default: 'gplas_results/'
    ## Settings:
    ##       -n      project name given to gplas. Default: 'unnamed'
    ##   -s      bacterial species corresponding to the graph file. If bacterial species corresponds to:
    ##                 'Enterococcus faecium','Klebsiella pneumoniae' or 'Escherichia coli' then prediction will be perfomed using mlplasmids. Default: 'unknown'
    ##   -t      threshold to predict plasmid-derived sequences. Default: 0.5
    ##   -x      Number of times gplas finds plasmid paths per each plasmid seed. Default: 10

Advanced usage
--------------

    ./gplas.sh -i mygraph.gfa -o foldername -n myfirstprediction -s 'Helicobacter pylori' -t 0.5 -x 10
