#
# These are functions that help manage input/output for LEfSe.
#

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
  data <- read.table(fp, sep="\t", header = FALSE, stringsAsFactors = FALSE, na.strings = c("", "-"))
  colnames(data) <- c("Feature", "LogHighestMean", "ClassHighestMean", "LogLDAScore", "WilcoxRes")
  data
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
# weights_fp: path to weighted matrix file
# metadata_fp: path to CSV sample attributes file
# column_name: name of column in metadata_fp to use as LEfSe class
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
  # Substitute NA values in the grouping column
  idxl <- is.na(lefse_metadata[[column_name]])
  lefse_metadata[[column_name]][idxl] <- column_na_txt
  # Combine and transpose all the information, so metadata is the first few rows
  # and measurements are the rows that follow
  combo <- rbind(t(lefse_metadata), t(weights))
  if (! is.null(out_path)) {
    save_lefse_input(combo, out_path)
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


# Other Functions ---------------------------------------------------------


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
    data <- lefse_load_res(row[["Path"]])
    if (row[["Category"]] == "pathways") {
      data$Feature <- sub("^path_", "path:", data$Feature)
    } else if (row[["Category"]] == "modules") {
      data$Feature <- sub("^md_", "md:", data$Feature)
    } else if (row[["Category"]] == "enzymes") {
      data$Feature <- gsub("_", ".", sub("^ec_", "ec:", data$Feature), fixed = TRUE)
    } else {
      warning("unrecognized category/feature name")
    }
    data
  })
  names(res) <- res_fields$Case
  list(res_fields = res_fields, res = res)
}
