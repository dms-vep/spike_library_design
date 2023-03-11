"""``snakemake`` file that does library design."""


configfile: "config.yaml"


rule all:
    input:
        config["sequential_to_reference"],
        alignment_counts=config["alignment_counts"],


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
