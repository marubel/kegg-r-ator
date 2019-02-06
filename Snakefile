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

rule all_kegg_to_ko:
    input: 
        expand("ko/KO_kegg_df_{sample}.tsv", sample = SAMPLES)

rule all_kegg_pathways:
    input:
        expand("pathways/ko_kegg_pathway_{sample}.tsv", sample = SAMPLES)
    output:
        rbound = "pathways/all_pathway_summary.tsv",
        matrix = "pathways/all_weighted_pathway_matrix.tsv"
    run:
        R("""
	source("kegg-r-ator.R")
	infiles = strsplit("{input}", " ", fixed = TRUE)[[1]]
	agg_to_pathways(infiles, "{output.rbound}","{output.matrix}")
	""")

rule all_kegg_modules:
    input:
        expand("modules/ko_module_df_{sample}.tsv", sample = SAMPLES)
    output:
        rbound = "modules/all_module_summary.tsv",
        matrix = "modules/all_weighted_module_matrix.tsv"
    run:
        R("""
        source("kegg-r-ator.R")
        infiles = strsplit("{input}", " ", fixed = TRUE)[[1]]
        agg_to_modules(infiles, "{output.rbound}","{output.matrix}")
        """)

rule all_kegg_enzymes:
    input:
        expand("enzymes/ko_enzymes_df_{sample}.tsv", sample = SAMPLES)
    output:
        rbound = "enzymes/all_enzyme_summary.tsv",
        matrix = "enzymes/all_weighted_enzyme_matrix.tsv"
    run:
        R("""
        source("kegg-r-ator.R")
        infiles = strsplit("{input}", " ", fixed = TRUE)[[1]]
        agg_to_enzymes(infiles, "{output.rbound}", "{output.matrix}")
        """)

rule all_kegg:
    input:
        rules.all_kegg_pathways.output, rules.all_kegg_modules.output, rules.all_kegg_enzymes.output

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
        tsv = "ko/KO_kegg_df_{sample}.tsv"
    input:
        sample = "deduplicated/{sample}.m8",
        ko_mapping_file_fp = str(METADATA_ROOT / "ko_genes.list")
    run: R("""
            source("kegg-r-ator.R")
            KEGG_to_ko(
              "{input.sample}",
              "{input.ko_mapping_file_fp}",
              "{output.tsv}")
            """)
 
# Use the KEGG_to_ko ouput_{sample}.tsv files as input and the ko_to_pathways R
# function
rule ko_to_pathways:
    output: 
        tsv = "pathways/ko_kegg_pathway_{sample}.tsv"
    input: 
        KO_kegg_df_fp = "ko/KO_kegg_df_{sample}.tsv",
        map_title_fp = str(METADATA_ROOT / "map_title.tab"),
        ko_enzyme_fp = str(METADATA_ROOT / "ko_enzyme.list"),
        module_parsed_fp = str(METADATA_ROOT / "module_parsed.tsv"),
        ko_pathway_fp = str(METADATA_ROOT / "ko_pathway.list")
    run: R("""
	    source("kegg-r-ator.R")
	    ko_to_pathways(
	      KO_kegg_df_fp = "{input.KO_kegg_df_fp}",
	      map_title_fp = "{input.map_title_fp}",
	      module_parsed_fp = "{input.module_parsed_fp}",
	      ko_enzyme_fp ="{input.ko_enzyme_fp}",
	      ko_pathway_fp = "{input.ko_pathway_fp}",
              out_path = "{output.tsv}")
	   """)

# Use the KEGG_to_ko ouput_{sample}.tsv files as input and the ko_to_modules R
# function
rule ko_to_modules:
    output:
        tsv = "modules/ko_module_df_{sample}.tsv"
    input:
        KO_kegg_df_fp = "ko/KO_kegg_df_{sample}.tsv",
        ko_module_fp = str(METADATA_ROOT / "ko_module.list"),
        map_title_fp = str(METADATA_ROOT / "map_title.tab"),
        module_parsed_fp = str(METADATA_ROOT / "module_parsed.tsv"),
        ko_enzyme_fp = str(METADATA_ROOT / "ko_enzyme.list")
    run: R("""
	    source("kegg-r-ator.R")
	    ko_to_modules(
    	      KO_kegg_df_fp = "{input.KO_kegg_df_fp}",
	      map_title_fp = "{input.map_title_fp}",
	      module_parsed_fp = "{input.module_parsed_fp}",
              ko_enzyme_fp = "{input.ko_enzyme_fp}",
              ko_module_fp = "{input.ko_module_fp}",
              out_path = "{output.tsv}")
            """)

# Use the KEGG_to_ko ouput_{sample}.tsv files as input and the ko_to_enzymes R
# function

rule ko_to_enzymes:
    output:
        tsv = "enzymes/ko_enzymes_df_{sample}.tsv"
    input:
        KO_kegg_df_fp = "ko/KO_kegg_df_{sample}.tsv",
        map_title_fp = str(METADATA_ROOT / "map_title.tab"),
        module_parsed_fp = str(METADATA_ROOT / "module_parsed.tsv"),
        ko_enzyme_fp = str(METADATA_ROOT / "ko_enzyme.list")
    run: R("""
            source("kegg-r-ator.R")
            ko_to_enzymes(
              KO_kegg_df_fp = "{input.KO_kegg_df_fp}",
              map_title_fp = "{input.map_title_fp}",
              module_parsed_fp = "{input.module_parsed_fp}",
              ko_enzyme_fp = "{input.ko_enzyme_fp}", 
              out_path = "{output.tsv}")
            """)

# Take only the first line in the input file for each qseqid (first field).
# Assumes original sort order, with the best BLAST result first for each
# qseqid, is preserved.
rule deduplicate_m8:
    output: tsv = "deduplicated/{sample}.m8"
    input: tsv = str(DATA_ROOT / "{sample}.m8")
    # https://stackoverflow.com/a/1916188
    # This is faster than the answer using sort, but even more importantly,
    # uses way less RAM.
    shell: "awk '!_[$1]++' < {input} > {output}"


