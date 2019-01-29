from snakemake.utils import R

# This sort of thing can be split off into a configuration file for maximum
# configurableness.  But this is fine for now.
DATA_ROOT = "/media/lorax/users/marubel/CM_KEGG"

# A small example for now.
rule all:
    input: "KO_kegg_df_test.tsv"

# The R function from snakemake.utils allows R code inline as described here:
# https://snakemake.readthedocs.io/en/v3.11.0/snakefiles/utils.html#scripting-with-r
# Compare that with the "snakemake" S4 object available if run via a "script"
# directive as described here:
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts
#
# Notes with using R():
#  * Any text output gets caught as a RRuntimeWarning and wrapped up in a bulky
#    Python log message.  Turning off optional messages in the R code helps
#    clean things up.
#  * You still need to quote things like "{input.sample}" since R will see
#    those things as literal text

rule kegg_to_ko:
    output:
        tsv = "KO_kegg_df_{sample}.tsv"
    input:
        #sample = DATA_ROOT + "/{sample}.m8",
        sample = "{sample}.m8",
        ko_mapping_file_fp = DATA_ROOT + "/ko_genes.list"
    run: R("""
            source("kegg-r-ator.R")
            KEGG_to_ko(
              "{input.sample}",
              "{input.ko_mapping_file_fp}",
              "{output.tsv}")
            """)

rule test_m8:
    output: tsv = "test.m8"
    input: tsv = DATA_ROOT + "/D0041_1.m8"
    shell: "head {input.tsv} > {output.tsv}"
