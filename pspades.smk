configfile: "final.yaml"

rule quast_alignment:
    input:
        nodes=lambda wildcards: config["samples"][wildcards.sample],
        reference="reference_genome/{sample}_ref_genome.fasta"
    output:
        align=directory("pspades/{sample}_alignments")
    conda:
        "envs/quast.yaml"
    shell:
        "quast.py -R {input.reference} -a all -m 1000 -o {output.align} {input.nodes}"


rule awk_parsing_alignment:
    input:
        alignment=directory("pspades/{sample}_alignments")
    output:
        "pspades/{sample}_alignment_test.txt"
    shell:
        """awk '{{print $5,$6}}' {input.alignment}/contigs_reports/*_contigs.tsv > {output}"""

rule pspades_evaluation:
    input:
        nodes=lambda wildcards: config["samples"][wildcards.sample],
        alignments="pspades/{sample}_alignment_test.txt"
    output:
        completeness="pspades/{sample}_completeness.tab",
        precision="pspades/{sample}_precision.tab"
    conda:
        "envs/r_packages.yaml"
    params:
        iterations = config["number_iterations"],
        classifier = config["classifier"],
        species = config["species"],
        name = config["name"]
    script:
        "scripts/pspades_evaluation.R"
