configfile: "final.yaml"

rule awk_nodes:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mobrecon/{sample}_raw_nodes.fasta"
    conda:
        "envs/r_packages.yaml"
    message:
        "Extracting the nodes from the graph {input}"
    log:
        "mobrecon/{sample}_log_nodes.txt"
    shell:
        """awk '{{if($1 == "S") print "\>"$1$2"_"$4"_"$5"\\n"$3}}' {input}  1>> {output} 2>> {log}"""

rule awk_links:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mobrecon/{sample}_raw_links.txt"
    conda:
        "envs/r_packages.yaml"
    message:
        "Extracting the links from the graph {input}"
    log:
        "mobrecon/{sample}_log_links.txt"
    shell:
        """awk -F "\\t" '{{if($1 == "L") print $N}}' {input}  1>> {output} 2>> {log}"""

rule mobrecon:
    input:
        "mobrecon/{sample}_raw_nodes.fasta"
    output:
        directory("mobrecon/{sample}_results")
    message:
        "Running mob recon"
    conda:
        "envs/mobrecon.yaml"
    shell:
        "mob_recon --infile {input} --outdir {output} --run_typer"

rule quast_alignment:
    input:
        nodes="mobrecon/{sample}_raw_nodes.fasta",
        reference="reference_genome/{sample}_ref_genome.fasta"
    output:
        align=directory("mobrecon/{sample}_alignments")
    conda:
        "envs/quast.yaml"
    shell:
        "quast.py -R {input.reference} -a all -m 1000 -o {output.align} {input.nodes}"


rule awk_parsing_alignment:
    input:
        alignment=directory("mobrecon/{sample}_alignments")
    output:
        "mobrecon/{sample}_alignment_test.txt"
    shell:
        """awk '{{print $5,$6}}' {input.alignment}/contigs_reports/*_raw_nodes.tsv > {output}"""


rule cp_files:
    input:
        results=directory("mobrecon/{sample}_results")
    output:
        contigreport="mobrecon/{sample}_contig_report.txt"
    shell:
        "cp {input.results}/contig_report.txt {output.contigreport}"


rule hyasp_coverage:
    input:
        nodes="mobrecon/{sample}_raw_nodes.fasta",
	    links="mobrecon/{sample}_raw_links.txt",
    output:
        graph_contigs="mobrecon/{sample}_graph_contigs.tab",
    	graph_repeats="mobrecon/{sample}_repeats_graph.tab",
    	clean_links="mobrecon/{sample}_clean_links.tab"
    params:
        classifier = config["classifier"],
        threshold = config["threshold_prediction"]
    conda:
        "envs/r_packages.yaml"
    message:
        "Cleaning some of the inputs required for the evaluation of the tool"
    script:
        "scripts/mob_coverage.R"


rule mob_evaluation:
    input:
        nodes="mobrecon/{sample}_raw_nodes.fasta",
	    clean_links="mobrecon/{sample}_clean_links.tab",
        bins="mobrecon/{sample}_contig_report.txt",
        graph_contigs="mobrecon/{sample}_graph_contigs.tab",
        graph_repeats="mobrecon/{sample}_repeats_graph.tab",
        alignments="mobrecon/{sample}_alignment_test.txt"
    output:
        completeness="mobrecon/{sample}_completeness.tab",
        precision="mobrecon/{sample}_precision.tab"
    conda:
        "envs/r_packages.yaml"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        species = config["species"],
        name = config["name"]
    script:
        "scripts/mob_evaluation.R"
