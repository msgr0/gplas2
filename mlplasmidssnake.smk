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
    shell:
        """awk '{{if($1 == "S") print "\>"$1$2"_"$4"_"$5"\\n"$3}}' {input}  1>> {output} 2>> {log}"""

rule mlplasmids:
    input:
        "gplas_input/{sample}_raw_nodes.fasta"
    output:
        "mlplasmids_prediction/{sample}_plasmid_prediction.tab"
    params:
        species = config["species"],
        prob = config["threshold_prediction"]
    conda:
        "envs/r_packages.yaml"
    log:
        normalmessage="logs/{sample}_normal_log_mlplasmids.txt",
        errormessage="logs/{sample}_error_log_mlplasmids.txt"
    shell:
        "Rscript /home/sergi/mlplasmids/scripts/run_mlplasmids.R {input} {output} {params.prob} {params.species} 1>> {log.normalmessage} 2>> {log.errormessage}"

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
        classifier = config["classifier"]
    conda:
        "envs/r_packages.yaml"
    script:
        "/home/sergi/plasgraph/Scripts/gplas_coverage.R"

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
        connections="paths/{sample}_connections.csv"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"]
    conda:
        "envs/r_packages.yaml"
    threads: 1
    script:
        "/home/sergi/plasgraph/Scripts/gplas_paths.R"

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
        plot_graph="network/{sample}_plot_coocurrence_network.png",
        components="network/{sample}_components.csv"
    conda:
        "envs/r_packages.yaml"
    script:
        "/home/sergi/plasgraph/Scripts/gplas_coocurrence.R"


rule quast_alignment:
    input:
        nodes="gplas_input/{sample}_raw_nodes.fasta",
	reference="data/{sample}_reference.fasta"
    output:
        align=directory("evaluation/{sample}_alignments"),
    shell:
        "quast.py -R {input.reference} -a all -m 1000 -o {output.align} {input.nodes}"

rule awk_parsing_alignment:
    input:
        "evaluation/{sample}_alignments/contigs_reports/all_alignments_{sample}_raw_nodes.tsv"
    output:
        "evaluation/{sample}_alignment_test.txt"
    shell:
        """awk '{{print $5,$6}}' {input}'"""

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
        alignments="data/{sample}_parsing_alignments.tsv",
        solutions="paths/{sample}_solutions.csv",
        components="network/{sample}_components.csv"
    output:
        completeness="evaluation/{sample}_completeness.tab",
        precision="evaluation/{sample}_precision.tab"
    params:
        iterations = config["number_iterations"]
    script:
        "/home/sergi/plasgraph/Scripts/gplas_evaluation.R"
