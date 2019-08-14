configfile: "final.yaml"

rule awk_links:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "hyasp/{sample}_raw_links.txt"
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
        "hyasp/{sample}_raw_nodes.fasta"
    conda:
        "envs/r_packages.yaml"
    message:
        "Extracting the nodes from the graph {input}"
    log:
        "logs/{sample}_log_nodes.txt"
    shell:
        """awk '{{if($1 == "S") print "\>"$1$2"_"$4"_"$5"\\n"$3}}' {input}  1>> {output} 2>> {log}"""

rule hyasp_map:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "hyasp/{sample}.csv"
    params:
        species = config["species"],
        prob = config["threshold_prediction"]
    message:
        "Running hyasp map"
    shell:
        "hyasp map /home/sergi/Programs/hyasp/databases/ncbi_database_genes.fasta {output} -g {input}"

rule hyasp_filter:
    input:
        "hyasp/{sample}.csv"
    output:
        "hyasp/{sample}_filtered.csv"
    params:
        species = config["species"],
        prob = config["threshold_prediction"]
    message:
        "Running hyasp filter"
    shell:
        "hyasp filter /home/sergi/Programs/hyasp/databases/ncbi_database_genes.fasta {input} {output}"

rule hyasp_find:
    input:
        graph=lambda wildcards: config["samples"][wildcards.sample],
        mapfiltered="hyasp/{sample}_filtered.csv"
    output:
        results=directory("hyasp/{sample}_results")
    params:
        species = config["species"],
        prob = config["threshold_prediction"]
    message:
        "Running hyasp find"
    shell:
        "hyasp find {input.graph} /home/sergi/Programs/hyasp/databases/ncbi_database_genes.fasta {input.mapfiltered} {output.results} -b 1.0"



rule quast_alignment:
    input:
        nodes="hyasp/{sample}_raw_nodes.fasta",
        reference="reference_genome/{sample}_ref_genome.fasta"
    output:
        align=directory("hyasp/{sample}_alignments")
    conda:
        "envs/quast.yaml"
    shell:
        "quast.py -R {input.reference} -a all -m 1000 -o {output.align} {input.nodes}"

rule awk_parsing_alignment:
    input:
        alignment=directory("hyasp/{sample}_alignments")
    output:
        "hyasp/{sample}_alignment_test.txt"
    shell:
        """awk '{{print $5,$6}}' {input.alignment}/contigs_reports/*_raw_nodes.tsv > {output}"""

rule cp_files:
    input:
        results=directory("hyasp/{sample}_results")
    output:
        chains="hyasp/{sample}_contig_chains.csv",
        bins="hyasp/{sample}_plasmid_bins_putative.csv",
        seq="hyasp/{sample}_putative_plasmids.fasta"
    shell:
        "cp {input.results}/contig_chains.csv {output.chains} && cp {input.results}/plasmid_bins_putative.csv {output.bins} && cp {input.results}/putative_plasmids.fasta {output.seq}"

rule circularity:
    input:
        seq="hyasp/{sample}_putative_plasmids.fasta"
    output:
        circular="hyasp/{sample}_circular_seq.txt"
    shell:
        "grep '>' {input.seq} > {output.circular}"


rule hyasp_coverage:
    input:
        nodes="hyasp/{sample}_raw_nodes.fasta",
	    links="hyasp/{sample}_raw_links.txt",
    output:
        graph_contigs="hyasp/{sample}_graph_contigs.tab",
    	graph_repeats="hyasp/{sample}_repeats_graph.tab",
    	clean_links="hyasp/{sample}_clean_links.tab"
    params:
        classifier = config["classifier"],
        threshold = config["threshold_prediction"]
    conda:
        "envs/r_packages.yaml"
    message:
        "Cleaning some of the inputs required for the evaluation of the tool"
    script:
        "scripts/hyasp_coverage.R"

rule gplas_evaluation:
    input:
        nodes="hyasp/{sample}_raw_nodes.fasta",
	    clean_links="hyasp/{sample}_clean_links.tab",
        chains="hyasp/{sample}_contig_chains.csv",
        bins="hyasp/{sample}_plasmid_bins_putative.csv",
        graph_contigs="hyasp/{sample}_graph_contigs.tab",
        graph_repeats="hyasp/{sample}_repeats_graph.tab",
        alignments="hyasp/{sample}_alignment_test.txt",
        circularity="hyasp/{sample}_circular_seq.txt"
    output:
        completeness="hyasp/{sample}_completeness.tab",
        precision="hyasp/{sample}_precision.tab"
    conda:
        "envs/r_packages.yaml"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        species = config["species"],
        name = config["name"]
    script:
        "scripts/hyasp_evaluation.R"
