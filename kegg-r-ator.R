######################
###Define libraries###
######################

# conda only provides some of this.  we could also wrap these in a try/catch
# and handle package installs here.
library(magrittr,   warn.conflicts = FALSE, quietly = TRUE)
library(plyr,       warn.conflicts = FALSE, quietly = TRUE)
library(dplyr,      warn.conflicts = FALSE, quietly = TRUE)
library(readr,      warn.conflicts = FALSE, quietly = TRUE)
library(stringr,    warn.conflicts = FALSE, quietly = TRUE)
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
ko_to_modules <- function(map_title_fp, KO_kegg_df_fp, module_parsed_fp, ko_enzyme_fp, ko_module_fp, out_path) {
  ko_constants <- ko_format_constants(map_title_fp, module_parsed_fp, ko_enzyme_fp) 
  
  KO_kegg_df <- read_tsv(KO_kegg_df_fp)

  kegg_modules <- ko_module_fp %>% 
    read_tsv(col_names = c("ko_gene_id", "ModuleID")) %>%
    left_join(ko_constants$module_names, by="ModuleID") %>%
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
ko_to_enzymes <- function(KO_kegg_df_fp, map_title_fp, module_parsed_fp, ko_enzyme_fp, out_path){
  ko_constants <- ko_format_constants(map_title_fp, module_parsed_fp, ko_enzyme_fp) 

  KO_kegg_df <- read_tsv(KO_kegg_df_fp)

  kegg_enzymes <- ko_constants$kegg_enzymes %>%
    group_by(ko_gene_id) %>%
    mutate(Weight = 1 / n()) %>%
    ungroup()
  
  ko_kegg_enzyme <- KO_kegg_df %>%
    left_join(kegg_enzymes, by="ko_gene_id", "EcNumber") %>%
    group_by(ko_gene_id) %>%
    mutate(weighted_kegg_enzymes = Weight*num_ko) %>%
    ungroup() 
	
  write.table(ko_kegg_enzyme, file=out_path, sep='\t', quote=F, row.names = F)
}

#######################################
###Aggregate ko_to_pathways function###
#######################################
# helper func to rename and select only the weighted abundance 
# and ID cols
rename_select_summarise_pathways <- function(df) {
    sampid <- df$sample[1]
    df <- dplyr::select(df, c("PathwayID", "weighted_kegg_pathways")) %>%
      group_by(PathwayID) %>% 
      summarise(total_weighted = sum(weighted_kegg_pathways)) %>% 
      dplyr::select(c("PathwayID", "total_weighted"))
    dplyr::rename(df, !!sampid:=total_weighted)
}

agg_to_pathways=function(paths, out_path, matrix_path) {
    pathwaylist= lapply(paths, function (x) fread(file = x, sep="\t", header = TRUE, showProgress=FALSE))
    #outdf = data.frame(ko_gene_id = c(NA))
    #rename(df,!!samp_id:=weighted_kegg_pathways)
    #lapply(pathwaylist, function(x))
    all_dfs <- lapply(pathwaylist, rename_select_summarise_pathways)
    p1 <- join_all(all_dfs, by="PathwayID", type = "full")
    #p1 <- Reduce(function(x,y) full_join(x,y,by="ko_gene_id"), all_dfs)
    write.table(p1, file=matrix_path, sep='\t', quote = F, row.names = F)

    p2 <- Reduce(function(x,y) merge(x,y, all=TRUE), pathwaylist)
    write.table(p2, file=out_path, sep='\t', quote=F, row.names = F)
}

######################################
###Aggregate ko_to_modules function###
######################################
rename_select_summarise_modules <- function(df) {
    sampid <- df$sample[1]
    df <- dplyr::select(df, c("ModuleID", "weighted_kegg_modules")) %>%
      group_by(ModuleID) %>%
      summarise(total_weighted = sum(weighted_kegg_modules)) %>%
      dplyr::select(c("ModuleID", "total_weighted"))
    dplyr::rename(df, !!sampid:=total_weighted)
}

agg_to_modules=function(paths, out_path, matrix_path) {
    modulelist= lapply(paths, function (x) fread(file = x, sep="\t", header = TRUE, showProgress=FALSE))

    all_dfs <- lapply(modulelist, rename_select_summarise_modules)
    p1 <- join_all(all_dfs, by="ModuleID", type = "full")
    #p1 <- Reduce(function(x,y) full_join(x,y,by="ko_gene_id"), all_dfs)
    write.table(p1, file=matrix_path, sep='\t', quote = F, row.names = F)

    p2 <- Reduce(function(x,y) merge(x,y, all=TRUE), modulelist)
    write.table(p2, file=out_path, sep='\t', quote=F, row.names = F)
}

######################################
###Aggregate ko_to_enzymes function###
######################################
rename_select_summarise_enzymes <- function(df) {
    sampid <- df$sample[1]
    df <- dplyr::select(df, c("EcNumber", "weighted_kegg_enzymes")) %>%
      group_by(EcNumber) %>%
      summarise(total_weighted = sum(weighted_kegg_enzymes)) %>%
      dplyr::select(c("EcNumber", "total_weighted"))
    dplyr::rename(df, !!sampid:=total_weighted)
}

agg_to_enzymes=function(paths, out_path, matrix_path) {
    enzymelist= lapply(paths, function (x) fread(file = x, sep="\t", header = TRUE, showProgress=FALSE))

    all_dfs <- lapply(enzymelist, rename_select_summarise_enzymes)
    p1 <- join_all(all_dfs, by="EcNumber", type = "full")
    #p1 <- Reduce(function(x,y) full_join(x,y,by="ko_gene_id"), all_dfs)
    write.table(p1, file=matrix_path, sep='\t', quote = F, row.names = F)

    p2 <- Reduce(function(x,y) merge(x,y, all=TRUE), enzymelist)
    write.table(p2, file=out_path, sep='\t', quote=F, row.names = F)
}

###############################################################
###Append metadata variables of interest to aggregated files###
###############################################################
#metadata_13019 <- "/home/rubel/01_CameroonShotgun/Shotgun_FULL_CM_metadata_2019-01-30.txt" %>%
#  read_tsv(col_names = TRUE) %>%
#  rename(Sample = sampleID) 

#metadata_13019 <- as.data.frame(metadata_13019)
 
#metadata_13019_2 <- metadata_13019[,-1] 
#rownames(metadata_13019_2) <- metadata_13019[,1] 
#metadata_13019 <- metadata_13019_2 
#metadata_13019 <- as.data.frame(t(metadata_13019))
#rm(metadata_13019_2)

#agg_to_[something]$[metadata_variable] <- metadata_13019[rownames(agg_to_[something]), "[metadata_variable]" 
