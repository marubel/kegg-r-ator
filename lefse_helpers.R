#
# These are functions that help manage input/output for LEfSe.
#

pkgs <- c("ggplot2", "pheatmap", "yaml")
notloaded <- ! sapply(pkgs, require, character.only = TRUE, quietly = TRUE)
if (sum(notloaded)) {
  install.packages(pkgs[notloaded])
}

# LEfSe Helpers -----------------------------------------------------------


# Load a LEfSe .res file into a data frame.
#
# From http://huttenhower.sph.harvard.edu/galaxy/:
#
# "The output consists of a tabular file listing all the features, the logarithm
# value of the highest mean among all the classes, and if the feature is
# discriminative, the class with the highest mean and the logarithmic LDA score."
#
# But that's just four columns.  What about the fifth?
# From run_lefse.py:
#
# 97     outres = {}
# 98     outres['lda_res_th'] = lda_res_th
# 99     outres['lda_res'] = lda_res
# 100     outres['cls_means'] = cls_means
# 101     outres['cls_means_kord'] = kord
# 102     outres['wilcox_res'] = wilcoxon_res
lefse_load_res <- function(fp) {
  if (file.size(fp) == 0) {
    return(NULL)
  }
  data <- read.table(fp, sep="\t", header = FALSE, stringsAsFactors = FALSE, na.strings = c("", "-"))
  colnames(data) <- c("Feature", "LogHighestMean", "ClassHighestMean", "LogLDAScore", "WilcoxRes")
  data$Feature <- lefse_unmangle_features(data$Feature)
  data
}

# specific corrections for our feature names.
lefse_unmangle_features <- function(vec) {
  if (any(grepl("^path_", vec))) {
    vec <- sub("^path_", "path:", vec)
  } else if (any(grepl("^md_", vec))) {
    vec <- sub("^md_", "md:", vec)
  } else if (any(grepl("^ec_", vec))) {
    vec <- gsub("_", ".", sub("^ec_", "ec:", vec), fixed = TRUE)
  } else {
    warning("unrecognized category/feature name")
  }
  vec
}

# create a LEfSe-style grouped bar graph with scores per feature.
lefse_plot_bars <- function(data) {
  data <- subset(data, ! is.na(ClassHighestMean))
  data <- data[order(data$ClassHighestMean, data$LogLDAScore), ]
  data$pos <- 1:nrow(data)
  plt <- ggplot2::ggplot(data) +
    ggplot2::geom_bar(aes(x = pos,
                          y = LogLDAScore,
                          fill = ClassHighestMean),
                      stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::labs(x=NULL, y = "LDA SCORE (log 10)") +
    ggplot2::scale_x_discrete(limits = data$Feature) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
  plt
}

# Take a TSV file with a matrix of observations across samples, a CSV file with
# sample metadata, and a column name to use from the sample metadata, and create
# a combined TSV file that LEfSe/format_input.py will understand.  Currently
# supports just one class, no subclass.
# If all the class values end up being the same, an empty file is written, to
# signal that lefse should be skipped for this case.
# weights_fp: path to weighted matrix file
# metadata_fp: path to CSV sample attributes file
# column_name: name of column in metadata_fp to use as LEfSe class.  If NULL,
# samples with NA for their metadata value are removed.
# column_name_txt: what should NA values become in the column named above?
# out_path: output file path
prepare_lefse_input <- function(weights_fp,
                                metadata_fp,
                                column_name,
                                column_na_txt = NA,
                                out_path = NULL,
                                sample_id_name = "sampleID") {
  weights <- load_weights_tsv(weights_fp)
  s_attrs <- load_sample_attrs(metadata_fp)
  rowidx <- match(rownames(weights), s_attrs[[sample_id_name]])
  if (any(is.na(rowidx))) {
    stop("sample ID mismatch")
  }
  lefse_metadata <- data.frame(
     s_attrs[[sample_id_name]][rowidx],
     s_attrs[[column_name]][rowidx],
     stringsAsFactors = FALSE,
     check.names = FALSE
  )
  # Set column names
  colnames(lefse_metadata) <- c(sample_id_name, column_name)
  # Substitute NA values in the grouping column, or remove NA entries entirely
  idxl <- is.na(lefse_metadata[[column_name]])
  if (is.null(column_na_txt)) {
    # remove entire rows (both metadata and weights) that had NA values for
    # column_name
    lefse_metadata <- lefse_metadata[! idxl, ]
    weights <- weights[! idxl, ]
    combo <- rbind(t(lefse_metadata), t(weights))
  } else {
    lefse_metadata[[column_name]][idxl] <- column_na_txt
    # Combine and transpose all the information, so metadata is the first few rows
    # and measurements are the rows that follow
    combo <- rbind(t(lefse_metadata), t(weights))
  }
  if (! is.null(out_path)) {
    if (length(unique(lefse_metadata[[column_name]])) == 1) {
      cat("", file = out_path)
    } else {
      save_lefse_input(combo, out_path)
    }
  }
  return(combo)
}

save_lefse_input <- function(data, fp) {
  write.table(x = data, file = fp, quote = FALSE, sep = "\t",
              row.names = TRUE, col.names = FALSE)
}


# File load/save helpers --------------------------------------------------


# Load a matrix of numeric values from a TSV file.
# input data:
# first column is name of the pathways/enzymes/modules values.  All subsequent
# columns are samples.  Each row is a single measurement.  Each cell is a value
# for a given sample+measurement.  Zeros are stored as NAs.
# output matrix:
# columns are measurements (within pathways/enzymes/modules)
# rows are samples
# zeros are zeros
# all cells are weight values
load_weights_tsv <- function(fp) {
  data <- data.table::fread(fp, header = TRUE, check.names = FALSE)
  na_rows <- which(apply(data, 1, function(row) all(is.na(row))))
  if (length(na_rows)) {
    warning("NA rows found, removing")
    data <- data[-na_rows, ]
  }
  rn <- sub("_1(\\.dedup)?", "", colnames(data)[-1])
  cn <- unlist(data[, 1])
  data <- data[, -1]
  data <- t(as.matrix(data))
  # (yes, it will actually set named vectors as row/column names if you let it)
  rownames(data) <- unname(rn)
  colnames(data) <- unname(cn)
  data[is.na(data)] <- 0
  data
}

# Load sample attributes from a CSV file.
load_sample_attrs <- function(fp) {
  # first column is just a row number.
  data <- data.table::fread(fp,
                            drop = 1,
                            colClasses = c(sampleID = "factor"))
  data <- as.data.frame(data)
  rownames(data) <- data$sampleID
  data
}

# load a simple list of terms from a text file into a character vector.
load_txt <- function(fp) {
  data <- read.csv(fp,
                   header = FALSE,
                   stringsAsFactors = FALSE)[, 1]
  names(data) <- data
  data
}

load_var_names <- function(fp) {
  # Handle empty-file case automatically
  if (file.size(fp) == 0) {
    data <- NULL
  } else {
    data <- data.table::fread(fp, col.names = c("VarID", "VarName"))
    data <- as.data.frame(data)
    rownames(data) <- data$VarID
  }
  data
}

load_config <- function(fp) {
  yaml::read_yaml(fp)
}


# Other Functions ---------------------------------------------------------


# plot measurement values (normalized to abundances per sample and z-scores per 
# measurement) present in the given LEfSe .res data, and order samples on rows
# by group and measurements on columns by LEfSe enrichment reported.  To produce
# a summary heatap with just "interesting" measurements, filter the lefse_data
# beforehand, for example on LogLDAScore.
# data: matrix of numeric values to plot.  rows are samples (should match 
#       s_attrs$sampleID), columns observations.
# lefse_data: a .res data frame from running LEfSe on the given data.
# s_attrs: sample attributes data frame
# md_var: column name in sample attributes that was used for LEfSe
# grp_colors: optional mapping of values to color strings for md_var
# var_names: optional data frame of substitutions for specific values in
# colnames of data and lefse_data$Feature
plot_heatmap_with_lefse <- function(data, lefse_data, s_attrs, md_var,
                                    grp_colors = NULL, var_names = NULL) {
  if (is.null(lefse_data) || nrow(lefse_data) == 0) {
    warning("No LEfSe data given for heatmap filtering")
    plot.new()
    plot.window(c(-1, 1), c(-1, 1))
    text(0, 0, "No LEfSe data given for heatmap filtering")
    return(plot.new())
  }

  # If custom variable names were given, subsitute those in for both the Feature
  # column of lefse_data and column names in the data matrix.
  if (! is.null(var_names)) {
    idx <- match(lefse_data$Feature, var_names$VarID)
    lefse_data$Feature <- ifelse(
      is.na(idx),
      lefse_data$Feature,
      var_names$VarName[idx])
    idx <- match(colnames(data), var_names$VarID)
    colnames(data) <- ifelse(
      is.na(idx),
      colnames(data),
      var_names$VarName[idx])
  }

  # Annotation for columns: which variable is enriched in which group?
  anno_col <- data.frame(
    EnrichedIn = factor(lefse_data$ClassHighestMean, levels = levels(factor(s_attrs[[md_var]]))),
    row.names = lefse_data$Feature,
    stringsAsFactors = FALSE
  )
  anno_col <- anno_col[order(anno_col$EnrichedIn), , drop = FALSE]
  anno_col <- anno_col[rownames(anno_col) %in% colnames(data), , drop = FALSE]
  
  # Annotation for rows: which sample belongs in which group?
  anno_row <- data.frame(x = factor(s_attrs[[md_var]]))
  colnames(anno_row) <- md_var
  rownames(anno_row) <- s_attrs$sampleID
  anno_row <- anno_row[order(anno_row[[md_var]]), , drop=FALSE]
  # Only keep annotation rows that are relevant for the given data
  anno_row <- anno_row[rownames(anno_row) %in% rownames(data), , drop = FALSE]
  # exclude NAs for md_var as well
  anno_row <- anno_row[! is.na(anno_row[[md_var]]), , drop = FALSE]
  
  # Standardize the levels to a consistent set but with no empties.
  lvls <- unique(c(as.character(anno_row[, md_var]),
                   as.character(anno_col[, "EnrichedIn"])))
  anno_row[[md_var]] <- factor(anno_row[[md_var]], levels = lvls)
  anno_col[["EnrichedIn"]] <- factor(anno_col[["EnrichedIn"]], levels = lvls)
  
  # Normalize per sample to unit sum, and scale per value column to z-score.
  data_norm <- t(apply(data, 1, function(row) row/sum(row)))
  data_norm <- scale(data_norm)
  
  # Take only the selected rows and columns
  data_notable <- data_norm[rownames(anno_row), rownames(anno_col), drop = FALSE]
  
  # Now, plot values
  
  # Colors: define a set of colors for each factor level in the grouping 
  # metadata variable.  Duplicate the set so we use the same ones for the rows
  # and columns.

  # If a color mapping list was given, see if the md_var given is present, and
  # if so, use the named colors.
  if (! is.null(grp_colors)) {
    colors <- grp_colors[levels(anno_row[[md_var]])]
    if (any(is.na(colors))) {
      idxl <- ! levels(anno_row[[md_var]]) %in% names(grp_colors)
      missings <- levels(anno_row[[md_var]])[idxl]
      stop(paste("missing color definitions for",
                 md_var,
                 ":", paste(missings, sep = ", ")))
    }
  } else {
    colors <- 1 + seq_along(levels(anno_row[[md_var]]))
  }
  names(colors) <- levels(anno_row[[md_var]])
  annotation_colors <- list(colors, colors)
  names(annotation_colors) <- c(md_var, "EnrichedIn")
  
  args <- list(mat = data_notable,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_row = anno_row,
               annotation_col = anno_col,
               annotation_colors = annotation_colors)
  
  do.call(pheatmap::pheatmap, args)
}

# Wrapper for plot_heatmap_with_lefse using file paths and performing filtering
# via LEfSe results.
# weights_fp: path to TSV matrix of observations
# res_fp: path to TSV LEfSe results file
# metadata_fp: path to CSV file with sample metadata
# column_name: sample metadata column to use for groupings
# out_path: path to output PDF file
# lda_score_min: minimum LogLDAScore in LEfSe results for variables to keep
# config_fp: path to YAML configuration file; used for grouping colors
# var_names_fp: path to TSV mapping of variable IDs to long names
# column_na_txt: value to filter out for ClassHighestMean in LEfSe results
plot_heatmap_with_lefse_files <- function(weights_fp, res_fp, metadata_fp,
                                          column_name, out_path,
                                          lda_score_min = 3,
                                          config_fp = NULL,
                                          var_names_fp = NULL,
                                          column_na_txt = "NA") {
  # If config was given, try to find a color mapping to use
  if (! is.null(config_fp)) {
    config <- load_config(config_fp)
    grp_colors <- config[["level_colors"]][[column_name]]
    # turn into named vector instead of list
    vecnames <- names(grp_colors)
    grp_colors <- as.character(grp_colors)
    names(grp_colors) <- vecnames
    if (length(grp_colors) == 0) {
      grp_colors <- NULL
    }
  } else {
    grp_colors <- NULL
  }
  # If var_names were given, apply the mapping to use more descriptive names.
  if (! is.null(var_names_fp)) {
    var_names <- load_var_names(var_names_fp)
  } else {
    var_names <- NULL
  }
  s_attrs <- load_sample_attrs(metadata_fp)
  weights <- load_weights_tsv(weights_fp)
  lefse_data <- lefse_load_res(res_fp)
  if (is.null(lefse_data)) {
    notables <- NULL
  } else {
    notables <- subset(lefse_data,
                       LogLDAScore > lda_score_min &
                         ! ClassHighestMean %in% c(column_na_txt))
  }
  # Found by eye that this height in inches seems to work well
  pdf(file = out_path, width = 11, height = nrow(weights) / 6)
  plot_heatmap_with_lefse(data = weights,
                          lefse_data = notables,
                          s_attrs = s_attrs,
                          md_var = column_name,
                          grp_colors = grp_colors,
                          var_names = var_names)
  dev.off()
}

# expand out a set of categories (pathways/modules/enzymes) and metadata
# variables, load a .res file for each combination, and reverse the expected
# modifications LEfSe has made to our category names.
lefse_load_res_all <- function(category_names, md_vars) {
  res_fields <- expand.grid(Prefix = "weights",
                            Category = category_names,
                            Group = md_vars,
                            Suffix = "res",
                            stringsAsFactors = FALSE)
  res_fields$Path <- file.path("lefse-results",
                               do.call(paste,
                                       c(res_fields, list(sep = "."))))
  res_fields <- subset(res_fields, select = -c(Prefix, Suffix))
  res_fields$Case <- with(res_fields, paste(Category, Group, sep = ":"))
  #res <- lapply(res_fields$Path, lefse_load_res)
  res <- apply(res_fields, 1, function(row) {
    lefse_load_res(row[["Path"]])
  })
  names(res) <- res_fields$Case
  list(res_fields = res_fields, res = res)
}
