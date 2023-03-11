"""``snakemake`` file that does library design."""


configfile: "config.yaml"


rule all:
    input:
        config["sequential_to_reference"],


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
