wildcard_constraints:
  sample="[^/]+"

rule awk_links:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "gplas_input/{sample}_raw_links.txt"
    conda:
        "../envs/r_packages.yaml"
    log:
        "logs/{sample}_log_links.txt"
    message:
        "Extracting the links from the graph {input}"
    shell:
        """
        awk -F "\\t" '{{if($1 == "L") print $N}}' {input}  1>> {output} 2>> {log}
        """

rule awk_nodes:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "gplas_input/{sample}_raw_nodes.fasta"
    params:
        min_node_length=config["min_node_length"]
    conda:
        "../envs/r_packages.yaml"
    log:
        "logs/{sample}_log_nodes.txt"
    message:
        "Extracting the nodes from the graph {input}"
    shell:
        """
        awk '{{if($1 == "S") print ">"$1$2"_"$4"_"$5"\\n"$3}}' {input}  1>> {wildcards.sample}_raw_nodes_unfiltered.fasta 2>> {log}
        awk -v min={params.min_node_length} 'BEGIN {{RS = ">" ; ORS = ""}} length($2) >= min {{print ">"$0}}' {output.nodes_unfiltered} > {output}
        rm {wildcards.sample}_raw_nodes_unfiltered.fasta
        """

rule gplas_coverage:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    links="gplas_input/{sample}_raw_links.txt",
        prediction="independent_prediction/{sample}_plasmid_prediction.tab"
    output:
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
    	graph_repeats="coverage/{sample}_repeats_graph.tab",
    	clean_links="coverage/{sample}_clean_links.tab",
    	clean_prediction="coverage/{sample}_clean_prediction.tab",
    	initialize_nodes="coverage/{sample}_initialize_nodes.tab"
    params:
        classifier = config["classifier"],
        threshold = config["threshold_prediction"]
    conda:
        "../envs/r_packages.yaml"
    message:
        "Extracting the sd k-mer coverage from the chromosome-predicted contigs"
    script:
        "../scripts/gplas_coverage.R"

rule gplas_paths:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    clean_links="coverage/{sample}_clean_links.tab",
        prediction="independent_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab"
    output:
        solutions="walks/normal_mode/{sample}_solutions.csv",
        connections="walks/normal_mode/{sample}_connections.tab"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        filt_gplas = config["filt_gplas"]
    conda:
        "../envs/r_packages.yaml"
    threads: 1
    message:
        "Searching for plasmid-like walks using a greedy approach"
    script:
        "../scripts/gplas_paths.R"

rule gplas_paths_bold:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
        clean_links="coverage/{sample}_clean_links.tab",
        prediction="independent_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab"
    output:
        solutions="walks/bold_mode/{sample}_solutions_bold.csv",
        connections="walks/bold_mode/{sample}_connections_bold.tab"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        filt_gplas = config["filt_gplas"],
        bold_sd_coverage = config["bold_sd_coverage"]
    conda:
        "../envs/r_packages.yaml"
    threads: 1
    message:
        "Searching for plasmid-like walks using a greedy approach"
    script:
        "../scripts/gplas_paths_bold.R"

rule gplas_coocurr:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    clean_links="coverage/{sample}_clean_links.tab",
        prediction="independent_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab",
        solutions="walks/normal_mode/{sample}_solutions.csv"
    output:
        plot_graph="results/normal_mode/{sample}_plasmidome_network.png",
        components="results/normal_mode/{sample}_bins.tab",
        results="results/normal_mode/{sample}_results.tab"
    params:
        threshold = config["threshold_prediction"],
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        edge_gplas = config["edge_gplas"],
        sample = config["name"],
        modularity_threshold = config["modularity_threshold"]
    conda:
        "../envs/r_packages.yaml"
    message:
        "Generating weights for the set of new edges connecting plasmid unitigs"
    script:
        "../scripts/gplas_coocurrence.R"

rule extract_unbinned_solutions:
    input:
        results="results/normal_mode/{sample}_results.tab",
        bold_walks="walks/bold_mode/{sample}_solutions_bold.csv"
    output:
        unbinned_walks="walks/unbinned_nodes/{sample}_solutions_unbinned.csv"
    message:
        "Extracting unbinned nodes from the initial run"
    shell:
        """
        for node in $(grep Unbinned {input.results} | cut -f 1 -d ' '); do \
        grep -w "^${{node}}" {input.bold_walks} >> {output.unbinned_walks} || continue; \
        done
        """

rule combine_solutions:
    input:
        unbinned_walks="walks/unbinned_nodes/{sample}_solutions_unbinned.csv",
       	normal_walks="walks/normal_mode/{sample}_solutions.csv"
    output:
        combined_walks="walks/{sample}_solutions.csv"
    message:
        "Combinning walks from normal and bold modes"
    shell:
        """
        cat {input.unbinned_walks} {input.normal_walks} > {output.combined_walks}
        """ 

rule gplas_coocurr_final:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
        clean_links="coverage/{sample}_clean_links.tab",
        prediction="independent_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab",
        solutions="walks/{sample}_solutions.csv"
    output:
        plot_graph="results/{sample}_plasmidome_network.png",
        components="results/{sample}_bins.tab",
        results="results/{sample}_results.tab"
    params:
        threshold = config["threshold_prediction"],
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        edge_gplas = config["edge_gplas"],
        sample = config["name"],
        modularity_threshold = config["modularity_threshold"]
    conda:
        "../envs/r_packages.yaml"
    message:
        "Generating weights for the set of new edges connecting plasmid unitigs"
    script:
        "../scripts/gplas_coocurrence_final.R"
