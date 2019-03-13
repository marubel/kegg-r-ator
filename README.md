# KEGG to R annotatOR (KEGG-R-ator) 
Automate stuff in [R] with [KEGG].

[R]: https://www.r-project.org
[KEGG]: https://www.kegg.jp
---
title: "Keggrator_README.rmd"
author: "MAR"
date: "3/13/2019"
output: 
  html_document:
    keep_md: yes
    toc: yes  
---

```{r setup, echo=FALSE, message=FALSE}
knitr::opts_chunk$set(
  fig.retina=2,
  fig.width=6,
  fig.height=4
)
```
## Introduction

KEGG-R-ator is a [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension written in [snakemake] (http://snakemake.readthedocs.io) that creates and aligns a DIAMOND-built database of metagenomic reads to KEGG, builds output files for KEGG enzymes, modules, and pathways, and visualizes the results. KEGG-R-ator uses [conda](http://condo.io) to manage dependencies. 

Specifically, KEGG-R-ator automates the following tasks: 
* [sbx_gene_clusters](https://github.com/sunbeam-labs/sbx_gene_clusters) is used to make a   blastdb from merged, paired-end reads using Diamond and a fasta file of reference   
  sequences. 
* R build of KEGG gene numbers to **KEGG orthologs**
* R build of KEGG orthologs to **KEGG pathways** in genes
* R build of **KEGG enzymes** from EC numbers in enzyme.list
* R build to **KEGG modules** from ko modules.list

KEGG-R-ator output can be used to generate heatmaps of differentially enriched pathways/enzymes/modules between grouups: [LEfSe](https://github.com/ressy/LEfSe)

------
### Contributors 
Meagan Rubel ([@marubel](https://github.com/marubel)
Jesse Connell ([@ressy](https://github.com/ressy))
Louis Taylor ([@louiejtaylor](https://github.com/louiejtaylor))
