configfile: "final.yaml"

rule awk_links:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "gplas_input/{sample}_raw_links.txt"
    conda:
        "envs/r_packages.yaml"
    message:
        "Extracting the links from the graph {input}"
    log:
        "logs/{sample}_log_links.txt"
    shell:
        """awk -F "\\t" '{{if($1 == "L") print $N}}' {input}  1>> {output} 2>> {log}"""

rule awk_nodes:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "gplas_input/{sample}_raw_nodes.fasta"
    conda:
        "envs/r_packages.yaml"
    message:
        "Extracting the nodes from the graph {input}"
    log:
        "logs/{sample}_log_nodes.txt"
    message:
        "Extracting the nodes from the graph {input}"
    shell:
        """awk '{{if($1 == "S") print "\>"$1$2"_"$4"_"$5"\\n"$3}}' {input}  1>> {output} 2>> {log}"""

rule mlplasmids:
    input:
        "gplas_input/{sample}_raw_nodes.fasta"
    output:
        "mlplasmids_prediction/{sample}_plasmid_prediction.tab"
    params:
        species = config["species"],
        threshold = config["threshold_prediction"]
    conda:
        "envs/r_packages.yaml"
    log:
        normalmessage="logs/{sample}_normal_log_mlplasmids.txt",
        errormessage="logs/{sample}_error_log_mlplasmids.txt"
    message:
        "Running mlplasmids to obtain the plasmid prediction using the nodes extracted from the graph"
    shell:
        "Rscript scripts/run_mlplasmids.R {input} {output} {params.threshold} {params.species} 1>> {log.normalmessage} 2>> {log.errormessage}"

rule gplas_coverage:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    links="gplas_input/{sample}_raw_links.txt",
        prediction="mlplasmids_prediction/{sample}_plasmid_prediction.tab"
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
        "envs/r_packages.yaml"
    message:
        "Extracting the sd k-mer coverage from the chromosome-predicted contigs"
    script:
        "scripts/gplas_coverage.R"

rule gplas_paths:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    clean_links="coverage/{sample}_clean_links.tab",
        prediction="mlplasmids_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab"
    output:
        solutions="paths/{sample}_solutions.csv",
        connections="paths/{sample}_connections.tab"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        mode = config["mode"]
    conda:
        "envs/r_packages.yaml"
    threads: 1
    message:
        "Searching for paths with a congruent probability of being plasmid and having a similar k-mer coverage"

    script:
        "scripts/gplas_paths.R"

rule gplas_coocurr:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    clean_links="coverage/{sample}_clean_links.tab",
        prediction="mlplasmids_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab",
        solutions="paths/{sample}_solutions.csv"
    output:
        plot_graph="results/{sample}_plasmidome_network.png",
        components="results/{sample}_components.tab",
        results="results/{sample}_results.tab"
    params:
        threshold = config["threshold_prediction"],
        classifier = config["classifier"],
        sample = config["name"]
    conda:
        "envs/r_packages.yaml"
    message:
        "Creating a co-occurrence network and selecting significant associations between nodes."
    script:
        "scripts/gplas_coocurrence.R"

rule quast_alignment:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
        reference="reference_genome/{sample}_ref_genome.fasta"
    output:
        align=directory("evaluation/{sample}_alignments")
    conda:
        "envs/quast.yaml"
    shell:
        "quast.py -R {input.reference} -a all -m 1000 -o {output.align} {input.nodes}"

rule awk_parsing_alignment:
    input:
        alignment=directory("evaluation/{sample}_alignments")
    output:
        "evaluation/{sample}_alignment_test.txt"
    shell:
        """awk '{{print $5,$6}}' {input.alignment}/contigs_reports/*_raw_nodes.tsv > {output}"""

rule gplas_evaluation:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	    clean_links="coverage/{sample}_clean_links.tab",
        prediction="mlplasmids_prediction/{sample}_plasmid_prediction.tab",
        coverage="coverage/{sample}_estimation.txt",
        graph_contigs="coverage/{sample}_graph_contigs.tab",
        graph_repeats="coverage/{sample}_repeats_graph.tab",
        clean_prediction="coverage/{sample}_clean_prediction.tab",
        initialize_nodes="coverage/{sample}_initialize_nodes.tab",
        alignments="evaluation/{sample}_alignment_test.txt",
        solutions="paths/{sample}_solutions.csv",
        components="results/{sample}_components.tab"
    output:
        completeness="evaluation/{sample}_completeness.tab",
        precision="evaluation/{sample}_precision.tab"
    conda:
        "envs/r_packages.yaml"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        species = config["species"],
        name = config["name"]
    script:
        "scripts/gplas_evaluation.R"
