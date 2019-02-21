#
# This is currently a big messy glob of examples and attempts and dead ends.
# Take with a grain of salt!
#


# Custom Functions --------------------------------------------------------


# LDA.
# data: Assumes row names of matrix are sample identifiers, column names 
# variable identifiers, all values are observations for a single sample and
# variable.
# s_attrs: data frame of sample atrributes
# sample.id: character name of sample attributes column storing the sample names
# to match to data's rownames.
# group.id: character name of sample attributes column storing the grouping
# variable to use.
# keep.na: what should we do with NA values for the grouping column?
#  TRUE: make a new category named "NA"
#  FALSE: filter out those samples
#  default, NULL: make a new category named "NA" but complain about it
lda <- function(data, s_attrs, sample.id="SunbeamID", group.id="StudyGroup", keep.na=NULL) {
  idx <- match(rownames(data), s_attrs[[sample.id]])
  grp <- s_attrs[[group.id]][idx]
  # NA handling for values in the grouping column
  if (any(is.na(grp))) {
    if (is.null(keep.na)) {
      warning("NAs in grouping variable")
      grp[is.na(grp)] <- "NA"
    } else if (keep.na) {
      grp[is.na(grp)] <- "NA"
    } else {
      idx <- ! is.na(grp)
      data <- data[idx, ]
      grp <- grp[idx]
    }
  }
  grp <- factor(grp)
  grp <- droplevels(grp)
  if (length(levels(grp)) < 2) {
    stop("All samples belong to only one group.")
  }
  # Filter out any variables that aren't doing anything for us, that is to say,
  # those that have only one unique value or one unique value per group.
  # Otherwise lda() fails.
  
  # See also caret::nearZeroVar
  #keeps <- apply(data, 2, function(x) length(unique(x)) > 1 )
  keeps <- apply(data, 2, function(vec) ncol(table(grp, vec)) > length(levels(grp)))
  if (sum(keeps) < ncol(data)) {
    warning("Dropping variables with too low variance")
  }
  data <- data[, keeps]
  
  
  # Scale each measurement to mean 1 and sd 0.
  # This makes the scaling multipliers directly comparable.
  data <- scale(data)
  lda_out <- MASS::lda(x = data, grouping = grp)
  # Attach sample groupings and modified dataunity matrix
  lda_out$groups <- grp
  lda_out$data <- data
  # use the LDA scalings to convert data matrix from N-dimensional
  # measurement-space to (num groups - 1)-dimensional LD space
  lda_out$prediction <- predict(lda_out, data)
  lda_out$group.id <- group.id
  lda_out
}

# Create a data frame of linear discriminant scaling values and note which group
# had the highest mean for each measurement.
build_lda_table <- function(lda_data) {
  # which group has the highest mean across each measurement?
  groups_enriched <- unlist(apply(lda_data$means, 2, function(vec) {
    unname(which(max(vec) == vec)[1])
  }))
  enriched <- rownames(lda_data$means)[groups_enriched]
  names(enriched) <- names(groups_enriched)
  
  lda_effect <- with(lda_data,
                     data.frame(Measurement = rownames(scaling),
                                EnrichedIn = factor(enriched,
                                                    levels = levels(lda_data$groups)),
                                stringsAsFactors = FALSE))
  
  
  # What is the maximum absolute scaling value for each measurement across LDs?
  lda_effect$LDMaxMag <- apply(lda_data$scaling, 1, function(vec) max(abs(vec)))
  
  # Bind all individual LD scaling columns
  lda_effect <- cbind(lda_effect, lda_data$scaling)
  
  # Order by group and then by magnitude of scaling.
  # lda_effect <- lda_effect[order(lda_effect$EnrichedIn,
  #                                abs(lda_effect$LDMaxMag), decreasing = TRUE), ]
  # Order by magnitude of scaling.
  lda_effect <- lda_effect[order(lda_effect$LDMaxMag, decreasing = TRUE), ]
  lda_effect
}

# gg scatterplot of the samples projected onto the top two linear discriminants,
# color-coded by the grouping factor that was used.
plot_lda_scatter <- function(lds) {
  if (ncol(lds$prediction$x) >= 2) {
    data <- data.frame(
      LD1 = lds$prediction$x[, 1],
      LD2 = lds$prediction$x[, 2],
      Group = lds$groups
    )
    txt_title <- "LDA"
    if (! is.na(lds$group.id)) {
      txt_title <- paste("LDA by", lds$group.id)
    }
    ggplot(data) +
      geom_point(aes(x = LD1, y = LD2, col = Group)) +
      ggtitle(txt_title)  
  } else {
    data <- data.frame(
      LD1 = lds$prediction$x[, 1],
      Group = lds$groups
    )
    ggplot(data) + 
      geom_boxplot(aes(x = Group, y = LD1))
  }
}

plot_lda_bars <- function(ld_scalings, n=20, is_per_grp=FALSE) {
  if (is_per_grp) {
    top_scalings <- do.call(rbind,
                            lapply(split(ld_scalings,
                                         ld_scalings$EnrichedIn),
                                   function(chunk) {
      head(chunk, n=n)
    }))
  } else {
    idx <- order(ld_scalings$LDMaxMag, decreasing = TRUE)
    ld_scalings <- ld_scalings[idx, ]
    top_scalings <- head(ld_scalings, n=n)  
  }
  
  # assume all positive, but if we only have two groups, put the second
  # negative.
  lvls <- levels(ld_scalings$EnrichedIn)
  top_scalings$DirScaling <- abs(top_scalings$LDMaxMag)
  if (length(lvls) == 2) {
    top_scalings$DirScaling <- with(top_scalings, abs(LDMaxMag)*((EnrichedIn == lvls[2])*2 - 1))
  }
  
  idx <- with(top_scalings, order(EnrichedIn, DirScaling))
  top_scalings <- top_scalings[idx, ]
  top_scalings$pos <- 1:nrow(top_scalings)
  
  ylim <- max(abs(top_scalings$LDMaxMag))*c(-1, 1)*2.2
  
  # TODO this is relative to the axis but shouldn't be.  fix this.
  s <- rep(2, nrow(top_scalings))
  nudge_y <- 0
  if (length(lvls) == 2) {
    s <- (top_scalings$EnrichedIn == lvls[2]) + 1
    nudge_y <- c(-1,1)[s]*nudge_y
  }

  ggplot(top_scalings,
         aes(x = pos,
             y = DirScaling,
             fill = EnrichedIn)) +
    geom_col() +
    scale_y_continuous(limits = ylim, labels = abs) +
    coord_flip() +
    geom_text(label = top_scalings$Measurement,
              nudge_y = nudge_y,
              hjust = c("right", "left")[s]) +
    xlab("Measurement") +
    ylab("Magnitude of LD Scaling") +
    labs(title = "Most discriminating measurements between groups",
         fill = "Enriched In") +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom"
    )
}

columns_summary <- function(data) {
  na_cts <- sapply(data, function(x) sum(is.na(x)))
  uniques <- sapply(data, function(x) length(unique(x)))
  uniques_no_na <- sapply(data, function(x) length(unique(x[!is.na(x)])))
  cts <- data.frame(
    Variable = colnames(data),
    NACount = na_cts,
    Uniques = uniques,
    UniquesNoNA = uniques_no_na)
  cts <- cts[order(cts[, "NACount"]), ]
  cts
}


# Analysis ----------------------------------------------------------------


library(magrittr)
library(dplyr)
library(viridis)
library(pheatmap)
library(ggplot2)
source("lefse_helpers.R")

s_attrs <- load_sample_attrs("shotgun_metadata_2519.csv")
weights <- list(
  pathways = load_weights_tsv("pathways/all_weighted_pathways_matrix.tsv"),
  enzymes = load_weights_tsv("enzymes/all_weighted_enzymes_matrix.tsv"),
  modules = load_weights_tsv("modules/all_weighted_modules_matrix.tsv")
)

# These are metadata variables Meagan is particularly interested in evaluating
# with LEfSe.
metadata_vars <- load_txt("shotgun_metadata_cols_of_interest.txt")

# A list of two things: a data frame of info on what we ran LEfSe on, and a big
# ol' list of output data frames for each case.  I've kept these separate for
# now since the values across cases should never really be presented together, 
# But for ease of management it could all be combined into one big data frame.
res_all <- lefse_load_res_all(names(weights), metadata_vars)


# Analysis - Custom -------------------------------------------------------


# First off, make sure we're only working with columns that we have useful
# grouping data for.

s_attrs[rownames(weights$pathways), ] %>%
  columns_summary() %>%
  filter(Variable %in% metadata_vars) %>%
  arrange(match(Variable, metadata_vars)) ->
  metadata_vars_stats

metadata_vars_stats %>%
  filter(UniquesNoNA < 2) %>%
  pull(Variable) %>%
  as.character() ->
  metadata_vars_stub

if (length(metadata_vars_stub)) {
  msg <- c("Some metadata variables have fewer than",
               "two unique groups for the data given.",
               "Removing:", metadata_vars_stub)
  warning(do.call(paste, as.list(msg)))
  metadata_vars <- metadata_vars[! metadata_vars %in% metadata_vars_stub]
}

# For each matrix of measurements across samples (pathways/modules/enzymes), Run
# LDA across each of the chosen metadata columns. The output here is a list 
# across pathways/modules/enzymes, each containing a list across metadata 
# variables, each containing a list bundling all LDA results for that 
# combination.
datasets <- names(weights)
names(datasets) <- datasets
lda_all <- lapply(datasets, function(d) {
  lapply(metadata_vars, function(v) {
    within(list(), {
      cat("LDA", d, v, end = "\n", file = 2)
      dataset <- d
      variable <- v
      lds <- lda(data = weights[[dataset]],
                 s_attrs = s_attrs,
                 sample.id = "sampleID",
                 group.id = v,
                 keep.na = FALSE)
      ld_scalings <- build_lda_table(lds)
    })
  })
})

# Add in bar plots for each one.
lda_all <- lapply(lda_all, function(lda_set) {
  lapply(lda_set, function(data) {
    within(data, {
      plt_bars <- plot_lda_bars(ld_scalings, n = 5, is_per_grp=TRUE) +
        ggtitle(paste0("Most discriminating measurements between groups (", variable, ")"))
    })  
  })
})

# Save bar plots to disk.
for (set in lda_all) {
  for (entry in set) {
    cat("LDA Save Plot: ", entry$dataset, "/", entry$variable, end = "\n", file = 2)
    filename <- paste0("lda_bars_", entry$dataset, "_", entry$variable, ".pdf")
    ggsave(filename = filename,
           width = unit(8.5, "in"),
           height = unit(11, "in"),
           plot = entry$plt_bars)
  }
}
 

# An example
grp <- "High_Pos_Giardia"
w <- weights[["pathways"]]
wvar_meas <- "path:map00944"
data <- data.frame(
  sampleID = rownames(w),
  Group = s_attrs[rownames(w), grp],
  Measurement = w[, var_meas])
data <- filter(data, ! is.na(Group))
ggplot(data) + geom_boxplot(aes(x=Group, y=Measurement))
