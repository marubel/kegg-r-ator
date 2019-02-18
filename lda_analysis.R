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

# by: group and annotate samples by this s_attrs column.
plot_heatmap <- function(data, s_attrs, by, use_log=TRUE, logstep=0.1,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         anno_col=NULL, annotation_colors=NULL) {
  # sort and annotate rows
  anno_row <- data.frame(x = s_attrs[[by]], stringsAsFactors = FALSE)
  colnames(anno_row) <- by
  rownames(anno_row) <- s_attrs$sampleID
  anno_row <- anno_row[rownames(data), , drop = FALSE]
  anno_row <- anno_row[order(anno_row[[by]]), , drop=FALSE]
  data <- data[order(match(rownames(data), rownames(anno_row))), ]
  
  if (use_log) {
    data[data==0] <- NA
    data <- log10(data)
    x_min <- ceiling(min(data, na.rm=TRUE))
    x_max <- floor(max(data, na.rm=TRUE))
    breaks <- seq(x_min, x_max, logstep)
    colors <- viridis::viridis(length(breaks)-1)  
  }
  
  args <- list(mat = data,
               cluster_rows = cluster_rows,
               cluster_cols = cluster_cols,
               annotation_row = anno_row)

  if (use_log) {
    args <- c(args, list(
      breaks = breaks,
      color = colors
    ))
  }
  
  if (!is.null(anno_col)) {
    args <- c(args, list(
      annotation_col = anno_col
      ))
  }
  if (!is.null(annotation_colors)) {
    args <- c(args, list(
      annotation_colors = annotation_colors
    ))
  }

  do.call(pheatmap::pheatmap, args)
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

# group duplicate columns
group_duplicate_columns <- function(data) {
  # group together all column indexes by identical column values.
  sets <- list()
  for (i in 1:(ncol(data)-1)) {
    # i already counted, skip
    if (i %in% unlist(sets)) {
      next
    }
    # starting new set with i
    set <- i
    for (j in (i+1):ncol(data)) {
      if (identical(data[, i], data[, j])) {
        set <- c(set, j)
      }
    }
    sets <- c(sets, list(set))
  }
  sets <- lapply(sets, function(idx) colnames(data)[idx])
  sets
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



# Example: results for pathways <-> subsistence as barcharts
plt <- lefse_plot_bars(res_all$res[["pathways:Subsistence"]])
#ggsave("lda.pdf", plt, width = unit(8.5, "in"), height=unit(30, "in"))




# Based on LEfSe results, how do our raw values for pathways look?
# Take just columns for the notable features, and rows where Subsistence is set
notables <- subset(res_all$res[["pathways:Subsistence"]],
                   LogLDAScore > 3 & ! ClassHighestMean %in% c("NA", "0"))
columns <- notables$Feature
rows <- which(rownames(weights$pathways) %in% as.character(subset(s_attrs, ! is.na(Subsistence))$sampleID))


# Normalize per sample to unit sum, and scale per value column to z-score.
pathways_norm <- t(apply(weights$pathways, 1, function(row) row/sum(row)))
pathways_notable <- pathways_norm[rows, columns]
pathways_notable <- scale(pathways_notable)

# Now, plot values
anno_col <- data.frame(
  EnrichedIn = notables$ClassHighestMean,
  row.names = notables$Feature,
  stringsAsFactors = FALSE
)
anno_col <- anno_col[order(anno_col$EnrichedIn), , drop = FALSE]
colors <- c(Agropastoralist = 2, "Hunter-Gatherer" = 3, Pastoralist = 4)
annotation_colors <- list(Subsistence = colors, EnrichedIn = colors)
pathways_notable <- pathways_notable[, rownames(anno_col)]
plot_heatmap(pathways_notable,
             s_attrs,
             "Subsistence",
             anno_col = anno_col,
             annotation_colors = annotation_colors,
             use_log = FALSE)


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


# Try Heatmap -------------------------------------------------------------

# I'm still curious to see if I can leverage this to get a more concise and
# interesting heatmap visualization of the raw values.

# whoa, too many variables!
set <- "pathways"
grouping <- "Subsistence"
idxl <- !is.na(s_attrs[rownames(weights[[set]]), grouping])
data <- weights[[set]][idxl, ]
plot_heatmap(data, s_attrs, by = grouping)
#plot_heatmap(data, s_attrs, by = grouping, cluster_rows = TRUE, cluster_cols = TRUE)

# What if we use just the top variables by LDA scaling?

keeps <- head(lda_all[[set]][[grouping]]$ld_scalings, n = 10)
keeps <- keeps[order(keeps$EnrichedIn, keeps$LDMaxMag), ]
data <- data[, keeps$Measurement]

plot_heatmap(data, s_attrs, by = grouping)
plot_heatmap(data, s_attrs, by = grouping, cluster_rows = TRUE, cluster_cols = TRUE)

anno_col <- data.frame(EnrichedIn = keeps$EnrichedIn)
rownames(anno_col) <- colnames(data)

lvls <- levels(anno_col$EnrichedIn)
colors <- (1:length(lvls)) + 1
names(colors) <- lvls
plot_heatmap(data, s_attrs, by = grouping, anno_col = anno_col,
             annotation_colors = list(EnrichedIn = colors, Subsistence = colors))
plot_heatmap(data, s_attrs, by = grouping,
             cluster_rows = TRUE, cluster_cols = TRUE,
             anno_col = anno_col,
             annotation_colors = list(EnrichedIn = colors, Subsistence = colors))


# Scaling Z Scores --------------------------------------------------------

# looks pretty but didn't help much.

breaks <- seq(-10, 10, 0.5)

zscores <- do.call(rbind, lapply(lda_all$pathways, function(lda_for_var) {
  z <- scale(lda_for_var$lds$scaling[, 1])
  hist(z, breaks = breaks, plot=FALSE)$counts
}))

nonzeros <- zscores != 0
fromleft <- apply(nonzeros, 1, function(x) match(TRUE, x))
fromright <- apply(nonzeros, 1, function(x) match(TRUE, rev(x)))
outliers <- sort(apply(cbind(fromleft, fromright), 1, min))
zscores <- zscores[names(outliers), ]
zscores[zscores==0] <- NA
zscores <- log10(zscores)

labels <- breaks
idxl <- ((1:length(labels))-1)%%5 != 0
labels[idxl] <- ""
pheatmap(zscores,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = viridis(10),
         border_color = NA,
         labels_col = labels)

# Try summary -------------------------------------------------------------



# how many NA entries do we have for each metadata column?  sorted by lease to
# most.
columns_summary(s_attrs) %>%
  filter(Variable != "sampleID") ->
  metadata_columns_summary


# Take the top ten variables that have fewer than 10 unique values
metadata_columns_summary %>%
  filter(Uniques < 10) %>%
  pull(Variable) %>%
  head(n=30) %>%
  as.character() ->
  x
names(x) <- x

# Try an LDA run for each of these variables.
ld_results <- lapply(x, function(v) {
  within(list(), {
    lds <- lda(data = weights$pathways,
        s_attrs = s_attrs,
        sample.id = "sampleID",
        group.id = v,
        keep.na = FALSE)
    ld_scalings <- build_lda_table(lds)  
  })
})

# Now, if we look *across* these metadata variables, do any particular ones 
# stand out as having just one or a few pathways that have significant scaling
# multipliers?
# Let's pool all scaling values together for each case, and scale them to zero 
# mean and unit standard deviation.  Then we can easily glance across them.
ld_z <- lapply(ld_results, function(data) {
  scale(as.vector(data$lds$scaling))
})

# which ones have values more than N SD increments away from the mean?
sapply(ld_z, function(z) {
  sum(abs(z)> 8)
})


# Example - Pathways ------------------------------------------------------


# First off, we see pathways separate samples very nicely for multiple possible
# types of grouping variables.

lds <- lda(data = weights$pathways,
    s_attrs = s_attrs,
    sample.id = "sampleID",
    group.id = "Population",
    keep.na=TRUE)
plot_lda_scatter(lds) + ggtitle("LDA by Population (All Pathways)")
table(lds$groups, lds$prediction$class)

# or, collapsed down onto one LD axis:
data <- data.frame(
  sampleID <- rownames(lds$prediction$x),
  LD1 = lds$prediction$x[, 1],
  Population = s_attrs[rownames(lds$prediction$x), "Population"]
)
ggplot(data) + geom_boxplot(aes(x = Population, y = LD1)) + ggtitle("LD1 by Population (All Pathways)")


### What pathways give the most separating power for Population?

# get a data frame of sorted scaling values per pathways variable
ld_scalings <- build_lda_table(lds)
hist(ld_scalings$LDMaxMag)


## hmm...
for (fraction in seq(1, 0.1, -0.1)) {
  nc <- ceiling(fraction*nrow(ld_scalings))
  columns <- head(ld_scalings$Measurement, nc)
  data_fraction <- weights$pathways[, columns]
  lds_fraction <- lda(data = data_fraction, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)
  plot_lda_scatter(lds_fraction) +
    ggtitle(paste0("LDA by Population (", fraction*100, "% of pathways)")) ->
    #coord_cartesian(xlim = c(-300, 300), ylim=c(-100, 100)) ->
    p
  plot(p)
  #ggsave(paste0("plot", formatC(nrow(ld_scalings) - nc, digits=3, flag="0"), ".png"), p)
}

# post-scaling, what columns are identical to other columns?
dups <- group_duplicate_columns(lds$data)

# what happens if we trim these out?  shouldn't matter which in terms of the math.
data <- weights$pathways[, unlist(lapply(dups, `[`, 1))]
lds <- lda(data = data[, 1:176], s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)


# take just the top N entries per group.
# ld_scalings %>%
#   group_by(EnrichedIn) %>%
#   top_n(50, abs(Scaling)) ->
#   ld_scalings_top



# just take the top portion of the most-separating measurements.
size <- ceiling(nrow(ld_scalings)/2)
ld_scalings_top <- head(ld_scalings, size)

pathways_subset <- weights$pathways[, ld_scalings_top$Measurement]
pathways_rand <- weights$pathways[, sample(1:ncol(weights$pathways), size)]


lds_top <- lda(data = pathways_subset, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)
lds_rand <- lda(data = pathways_rand, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)

plot_lda_scatter(lds) + ggtitle("LDA by Population")
plot_lda_scatter(lds_top) + ggtitle("LDA by Population, portion of measurements")
plot_lda_scatter(lds_rand) + ggtitle("LDA by Population, portion of measurements (random)")




plot_heatmap(pathways_subset, s_attrs, "Population")


# huh, that doesn't look like much at first glance.  let's sanity check a few directly.

idx <- which(max(abs(ld_scalings$Scaling)) == abs(ld_scalings$Scaling))
path <- ld_scalings$Measurement[idx]
vec <- weights$pathways[, path]
grp <- s_attrs$Population[match(names(vec), s_attrs$sampleID)]
boxplot(vec ~ grp)
# ok, that's something, but it's not amazing.  looks like maybe no single
# vaiable from pathways is quite enough.


# Examples - Enzymes ------------------------------------------------------


lds <- lda(data = weights$enzymes, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)
plot_lda_scatter(lds)
# Those NA samples squish the rest of the display
lds <- lda(data = weights$enzymes, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=FALSE)
plot_lda_scatter(lds)

# ld_scalings <- build_lda_table(lds)
# ld_scalings %>%
#   group_by(EnrichedIn) %>%
#   top_n(10, abs(Scaling)) ->
#   ld_scalings_top
# enzymes_subset <- weights$enzymes[, ld_scalings_top$Measurement]
# plot_heatmap(enzymes_subset, s_attrs, "Population")


# Example - Module --------------------------------------------------------


lds <- lda(data = weights$modules, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)
plot_lda_scatter(lds)


# Example - Combined ------------------------------------------------------


# No reason we have to keep these variables separate.  They'll get turned into z
# scores for the analysis anyway.

combo <- weights$pathways
combo <- cbind(combo, weights$enzymes)
combo <- cbind(combo, weights$modules)

######## Wait, what?  how did that make it worse?
lds <- lda(data = combo, s_attrs = s_attrs, sample.id = "sampleID", group.id = "Population", keep.na=TRUE)
plot_lda_scatter(lds)
table(lds$groups, lds$prediction$class)
