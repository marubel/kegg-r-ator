from pathlib import Path
from snakemake.utils import R

# This sort of thing can be split off into a configuration file for maximum
# configurableness.  But this is fine for now.
DATA_ROOT = Path("/media/lorax/users/marubel/CM_KEGG")
METADATA_ROOT = DATA_ROOT
if not DATA_ROOT.exists():
    DATA_ROOT = Path("/home/rubel/01_CameroonShotgun/sunbeam_output/mapping/sbx_gene_family/species_prokaryotes")
    METADATA_ROOT = Path("/home/rubel/05_KEGG")


M8 = DATA_ROOT.glob('*.m8') # full path objects for all .m8 files
SAMPLES = [p.stem for p in M8] # sample names, figured out from .m8 files

# A small example for now.
rule example:
    input: "KO_kegg_df_test.tsv"

rule all_kegg_to_ko:
    input: expand("KO_kegg_df_{sample}.tsv", sample = SAMPLES)

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

# Use the (much smaller) deduplicated m8 files as input to the KEGG_to_ko R
# function.
rule kegg_to_ko:
    output:
        tsv = "KO_kegg_df_{sample}.tsv"
    input:
        #sample = str(DATA_ROOT / "{sample}.m8"),
        sample = "{sample}.dedup.m8",
        ko_mapping_file_fp = str(METADATA_ROOT / "ko_genes.list")
    run: R("""
            source("kegg-r-ator.R")
            KEGG_to_ko(
              "{input.sample}",
              "{input.ko_mapping_file_fp}",
              "{output.tsv}")
            """)

# Take only the first line in the input file for each qseqid (first field).
# Assumes original sort order, with the best BLAST result first for each
# qseqid, is preserved.
rule deduplicate_m8:
    output: tsv = "{sample}.dedup.m8"
    input: tsv = str(DATA_ROOT / "{sample}.m8")
    # https://stackoverflow.com/a/1916188
    # This is faster than the answer using sort, but even more importantly,
    # uses way less RAM.
    shell: "awk '!_[$1]++' < {input} > {output}"

rule test_m8:
    output: tsv = "test.m8"
    input: tsv = str(DATA_ROOT / "D0041_1.m8")
    shell: "head {input.tsv} > {output.tsv}"
