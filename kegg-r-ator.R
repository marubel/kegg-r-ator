######################
###Define libraries###
######################

# conda only provides some of this.  we could also wrap these in a try/catch
# and handle package installs here.
library(magrittr,   warn.conflicts = FALSE, quietly = TRUE)
library(dplyr,      warn.conflicts = FALSE, quietly = TRUE)
# data.table prints a status line unless we turn it off
suppressPackageStartupMessages(library(data.table))

################################
###Define constants/filepaths###
################################

#ko_mapping_file_fp <- "filepath to ko_mapping_file_fp"
#data_path <- "filepath to *.m8 files"
#map_title_fp <- "Pathway to map_title.tab" 
#ko_pathway_fp <- "Pathway to ko_pathway.list"
#ko_enzyme_fp <- "Pathway to ko_enzyme.list"
#ko_module_fp <- "Pathway to ko_module.list"
#module_parsed_fp <- "Pathway to module_parsed.tsv" #to parse module use parse_module.bh script located in /media/lorax/users/marubel/CM_KEGG/ 
#out_path <- "Pathway to output directory"

ko_mapping_file_fp <- "/media/lorax/users/marubel/CM_KEGG/ko_genes.list"
data_path <- "/media/lorax/users/marubel/CM_KEGG" #there are currently two *m8 files in this directory
map_title_fp <- "/media/lorax/users/marubel/CM_KEGG/map_title.tab" 
ko_pathway_fp <- "/media/lorax/users/marubel/CM_KEGG/ko_pathway.list"
ko_enzyme_fp <- "/media/lorax/users/marubel/CM_KEGG/ko_enzyme.list"
ko_module_fp <- "/media/lorax/users/marubel/CM_KEGG/ko_module.list"
module_parsed_fp <- "/media/lorax/users/marubel/CM_KEGG/module_parsed.tsv" #to parse module use parse_module.bh script located in /media/lorax/users/marubel/CM_KEGG/ 
#out_path <- "/media/lorax/users/marubel/CM_KEGG"
out_path <- "/media/lorax/users/jesse/CM_KEGG"

######################
###Format constants###
######################
ko_format_constants <- function(map_title_fp, module_parsed_fp, ko_enzyme_fp, out_path){
pathway_names <- map_title_fp %>%
  read_tsv(col_names = c("PathwayID", "PathwayName")) %>%
  mutate(PathwayID = paste0("path:map", PathwayID))

module_names <- module_parsed_fp %>%
  read_tsv(col_names = c("ModuleID", "ModuleName")) %>%
  mutate(ModuleID = paste0("md:", ModuleID))
  
kegg_enzymes <- ko_enzyme_fp %>%
  read_tsv(col_names = c("ko_gene_id", "EcNumber"))
  
  list(pathway_names = pathway_names, module_names = module_names, kegg_enzymes = kegg_enzymes)
}
  
########################
###Blast6out function###
########################
read_blast6out <- function(filepath) {
  
  sample_id <- sub("*.m8", "", basename(filepath))
  
  blastx <- fread(filepath, sep ="\t", header = FALSE, showProgress = FALSE)

  blastx <- blastx %>% 
    set_colnames(c("qseqid", "sseqid", "pident","qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send","e_value", "bit_score"))
  
  data.frame(sample = sample_id, blastx)
}

#########################
###KEGG_to_ko function###
#########################
KEGG_to_ko <- function(filepath, ko_mapping_file_fp, outfile, plots=FALSE) {
  ko_genes <- fread(file = ko_mapping_file_fp, sep="\t", header=FALSE, showProgress=FALSE) %>%
    set_colnames(c("ko_gene_id","sseqid"))
  
  kegg_df <- read_blast6out(filepath)
  
  if (plots) {
    kegg_df %>%
     ggplot(aes(x = pident)) + geom_histogram()
    kegg_df %>%
      ggplot(aes(x = length)) + geom_histogram()
  }
  
  kegg_df <- kegg_df %>%
    group_by(qseqid) %>%
    slice(1) %>%
    ungroup() %>%
    filter(pident >= 65, length >= 35) %>%
    filter(!grepl("hsa", sseqid))
  
  KO_kegg_df <- ko_genes %>%
    right_join(kegg_df, by = c("sseqid"), all.x= TRUE, all.y = TRUE)
  
  KO_kegg_df <- na.omit(KO_kegg_df)
  
    KO_kegg_df <- KO_kegg_df %>%
    group_by(sseqid, qseqid) %>%
    mutate(weight_ko = 1/n()) %>%
    ungroup() 
  
    KO_kegg_df <- KO_kegg_df %>%
    group_by(sample, ko_gene_id) %>%
    summarize(num_ko = sum(weight_ko)) %>%
    ungroup()
  
  write.table(KO_kegg_df, file=outfile, sep='\t', quote=F, row.names = F)
  
}

##############################
####ko_to_pathways function###
##############################
ko_to_pathways <- function(module_parsed_fp, ko_enzyme_fp, KO_kegg_df_fp, map_title_fp, ko_pathway_fp, out_path) { 
  ko_constants <- ko_format_constants(map_title_fp, module_parsed_fp, ko_enzyme_fp) 

  KO_kegg_df <- read_tsv(KO_kegg_df_fp)
  
  kegg_pathways <- ko_pathway_fp %>%
    read_tsv(col_names = c("ko_gene_id", "PathwayID")) %>%
    filter(str_detect(PathwayID, "path:map")) %>%
    left_join(ko_constants$pathway_names, by="PathwayID") %>%
    group_by(ko_gene_id) %>%
    mutate(Weight = 1 / n()) %>%
    ungroup()

  ko_kegg_pathway <- KO_kegg_df%>%
     left_join(kegg_pathways, by="ko_gene_id") %>%
     group_by(ko_gene_id) %>%
     mutate(weighted_kegg_pathways = Weight*num_ko) %>%
     ungroup()
	 
  write.table(ko_kegg_pathway, file=out_path, sep='\t', quote=F, row.names = F)

}

############################
###ko_to_modules function###
############################
ko_to_modules <- function(module_names, KO_kegg_df_fp, ko_module_fp, out_path) {
  ko_constants <- ko_format_constants(map_title_fp, module_parsed_fp, ko_enzyme_fp) 
  
  KO_kegg_df <- read_tsv(KO_kegg_df_fp)

  kegg_modules <- ko_module_fp %>% 
    read_tsv(col_names = c("ko_gene_id", "ModuleID")) %>%
    left_join(module_names, by="ModuleID") %>%
    group_by(ko_gene_id) %>%
    mutate(Weight = 1 / n()) %>%
    ungroup()
  
  ko_kegg_module <-KO_kegg_df %>%
    left_join(kegg_modules, by="ko_gene_id", "ModuleID") %>%
    group_by(ko_gene_id) %>%
    mutate(weighted_kegg_modules = Weight*num_ko) %>%
    ungroup()
   
  write.table(ko_kegg_module, file=out_path, sep='\t', quote=F, row.names = F)

}

############################
###ko_to_enzymes function###
############################ 
ko_to_enzymes <- function(kegg_enzymes, KO_kegg_df_fp, ko_enzyme_fp, out_path){
  ko_constants <- ko_format_constants(map_title_fp, module_parsed_fp, ko_enzyme_fp) 

  KO_kegg_df <- read_tsv(KO_kegg_df_fp)

  kegg_enzymes <- kegg_enzymes %>%
    group_by(ko_gene_id) %>%
    mutate(Weight = 1 / n()) %>%
    ungroup()
  
  ko_kegg_enzyme <- KO_kegg_df %>%
    left_join(kegg_enzymes, by="ko_gene_id", "EcNumber") %>%
    group_by(ko_gene_id) %>%
    mutate(weighted_kegg_enzymes = Weight*num_ko) %>%
    ungroup() 
	
  write.table(ko_kegg_pathway, file=out_path, sep='\t', quote=F, row.names = F)
}
