"""``snakemake`` file that does library design."""


configfile: "config.yaml"


rule all:
    input:
        config["sequential_to_reference"],
        config["alignment_counts"],
        "results/usher/translated_muts.tsv",
        "results/usher/translated_recent_muts.tsv",


rule sequential_to_reference:
    input:
        extended_spike=config["extended_spike"],
        reference_spike=config["reference_spike"],
    output:
        config["sequential_to_reference"],
    log:
        notebook="results/notebooks/sequential_to_reference.ipynb",
    notebook:
        "notebooks/sequential_to_reference.py.ipynb"


rule alignment_counts:
    params:
        table_url=config["alignment_count_url"],
    output:
        alignment_counts=config["alignment_counts"],
    script:
        "scripts/alignment_counts.py"


rule get_usher_mat:
    params:
        mat=config["usher_mat"],
        fasta=config["usher_fasta"],
        gtf=config["usher_gtf"],
    output:
        mat="results/usher/mat.pb.gz",
        fasta="results/usher/ref.fa",
        gtf="results/usher/ref.gtf",
    shell:
        """
        wget -O - {params.fasta} | gunzip -c > {output.fasta}
        wget -O - {params.gtf} | gunzip -c > {output.gtf}
        wget -O {output.mat} {params.mat}
        """


rule translate_mat:
    input:
        mat=rules.get_usher_mat.output.mat,
        fasta=rules.get_usher_mat.output.fasta,
        gtf=rules.get_usher_mat.output.gtf,
    output:
        tsv="results/usher/translated_muts.tsv",
    shell:
        """
        matUtils summary \
            -i {input.mat} \
            -g {input.gtf} \
            -f {input.fasta} \
            -t {output.tsv}
        """


rule translate_recent_mat:
    input:
        mat=rules.get_usher_mat.output.mat,
        fasta=rules.get_usher_mat.output.fasta,
        gtf=rules.get_usher_mat.output.gtf,
    output:
        recent_mat="results/usher/recent_mat.pb.gz",
        tsv="results/usher/translated_recent_muts.tsv",
    params:
        clades=",".join(config["usher_recent_clades"])
    shell:
        """
        matUtils extract -i {input.mat} -o {output.recent_mat} -c "{params.clades}"
        matUtils summary \
            -i {output.recent_mat} \
            -g {input.gtf} \
            -f {input.fasta} \
            -t {output.tsv}
        """
