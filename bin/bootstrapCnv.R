#!/usr/bin/env Rscript
# R script for exomedepth data (requires a readcount data object created by readCount.R)

# Install necessary packages if not already installed
install_packages_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
}

install_bioc_packages_if_missing <- function(packages) {
  if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(packages)
}

# Check for required packages and install if missing
needed_packages <- c("ExomeDepth", "GenomicRanges", "Rsamtools", "stringr", "dplyr", "moments")
install_packages_if_missing(needed_packages)

# Check for required Bioconductor packages and install if missing
bioconductor_packages <- c("Biostrings", "Rsamtools", "GenomicRanges", "GenomicAlignments")
install_bioc_packages_if_missing(bioconductor_packages)

# Load the necessary libraries
sapply(c("GenomicRanges", "ExomeDepth", "Rsamtools", "stringr", "dplyr", "moments"), library, character.only = TRUE)

# Functions
'%nin%' <- function(x, y) !('%in%'(x, y))

entropy <- function(target) {
  freq <- table(target) / length(target)
  vec <- as.data.frame(freq)[,2]
  vec <- vec[vec > 0]
  -sum(vec * log2(vec))
}

read_numeric_arg <- function(arg, default_value) {
  numeric_value <- as.numeric(arg)
  if (is.na(numeric_value)) {
    return(default_value)
  }
  return(numeric_value)
}

read_bed_to_granges <- function(bed) {
  bed_data <- read.table(bed, header = FALSE, fill = TRUE)
  return(GRanges(seqnames = bed_data$V1, IRanges(start = bed_data$V2 + 1, end = bed_data$V3)))
}

# Set regex for sex matching of sample
sex <- regex('_[MF]_')

# Main script
args <- commandArgs(trailingOnly = TRUE)
setwd(dirname(args[1]))
load(args[1])

num_tries <- read_numeric_arg(args[2], 100)
sample_size <- read_numeric_arg(args[3], 30)
range_from <- read_numeric_arg(args[4], 1)
range_to <- read_numeric_arg(args[5], length(testsamplenames))

bed_files <- args[which(endsWith(args, '.bed'))]

if (length(bed_files) > 0) {
  message('Getting Segment annotations...')
  ranges <- lapply(bed_files, read_bed_to_granges)
  annotations <- unlist(as(ranges, "GRangesList"))
}

# Function to analyse scores
analyze_scores <- function(scores) {
  scores %>%
    filter(is.na(common) | !common) %>%
    group_by(exon) %>%
    summarise(
      entropy = entropy(prop),
      skewness = skewness(prop),
      n = n(),
      prevalence = n() / length(testsamplenames)
    ) -> exon.entropy
  
  common_exons <- exon.entropy %>%
    filter(entropy > 2 & prevalence > 0.4) %>%
    pull(exon)
  
  scores %>%
    filter(!(exon %in% common_exons) & prop > 0.5) %>%
    select(sample)
}

# Function to aggregate CNVs
aggregate_cnvs <- function(samplecnvs, exons, annotations) {
  dups <- table(unlist(apply(samplecnvs[samplecnvs$reads.ratio > 1, ], 1, function(x) x[1]:x[2]))) / tries
  dels <- table(unlist(apply(samplecnvs[samplecnvs$reads.ratio < 1, ], 1, function(x) x[1]:x[2]))) / tries
  
  rois <- unique(c(names(dups), names(dels)))
  results <- lapply(rois, function(roi) {
    exon.name <- exons$names[as.numeric(roi)]
    cnvs.ovp <- ifelse(is.null(annotations), NA, overlapsAny(exons[as.numeric(roi)], annotations))
    
    data <- data.frame(name = testsample, exon = exon.name, prop = NA, common = cnvs.ovp)
    if (!is.na(dups[roi])) {
      data$type <- "duplication"
      data$prop <- dups[roi]
    } else {
      data$type <- "deletion"
      data$prop <- dels[roi]
    }
    data
  })
  
  do.call(rbind, results)
}

# Check if there are scores files
score_files <- args[which(endsWith(args, '.scores'))]
if (length(score_files) > 0) {
  proceed <- "RUN" %in% args
  
  message('Analysing previous runs...')
  scores <- do.call(rbind, lapply(score_files, function(file) read.table(file, header = FALSE)))
  colnames(scores) <- c('sample', 'exon', 'type', 'prop', 'common')
  
  exclusions <- analyze_scores(scores)
  exclusions <- unique(exclusions$sample)
  
  if (length(exclusions) == 0) {
    stop("Good set of normals already")
  } else {
    if (!proceed) stop('To proceed, run same command and append RUN as an argument')
    testsamplenames <- setdiff(testsamplenames, exclusions)
    refsamplenames <- testsamplenames
  }
}

from <- ifelse(!is.na(as.numeric(args[4])), as.numeric(args[4]), 1)
to <- ifelse(!is.na(as.numeric(args[5])), as.numeric(args[5]), length(testsamplenames))
message(paste('Read Count Data:', args[1]))

if (!identical(refsamplenames, testsamplenames)) {
  stop("You must prepare data in batch normalisation mode.")
}

timestamp <- format(Sys.time(), "%y%m%d%H%M%S")
cnv.file <- sub("[.][^.]*$", paste0(".", timestamp, ".scores"), args[1])
data.file <- sub("[.][^.]*$", paste0(".", timestamp, ".RData"), args[1])

for (testsample in testsamplenames[from:to]) {
  cat(sprintf("%d/%d\n", i + from, to))
  i <- i + 1

  samplecnvs <- data.frame()
  results <- data.frame(name = character(),
                        exon = character(),
                        type = character(),
                        prop = numeric(),
                        common = logical())

  refsamples <- refsamplenames[refsamplenames != testsample]
  cat(sprintf("Reference samples: %s\n", toString(refsamplenames)))

  hasXY <- any(rois[, 1] == "X") || any(rois[, 1] == "Y")
  tssx <- str_extract(testsample, sex)
  rcsx <- str_extract(refsamples, sex)
  hasSex <- !(is.na(tssx) || any(is.na(rcsx)))

  sampling_msg <- if (hasXY && hasSex) {
    refsamples <- refsamples[rcsx == tssx]
    sprintf("==> SEX_MATCHED SAMPLING %d of %d", samplesize, length(refsamples))
  } else {
    sprintf("==> SAMPLING %d of %d", samplesize, length(refsamples))
  }
  cat(sampling_msg, "\n")

  for (t in 1:tries) {
    cat(sprintf("** %d of %d\n", t, tries))

    refsamples.sample <- sample(refsamples, samplesize, replace = TRUE)
    refcounts <- as.matrix(counts[, refsamples.sample])
    refset <- suppressWarnings(select.reference.set(
      test.counts = counts[, testsample],
      reference.counts = refcounts,
      bin.length = counts$end - counts$start
    ))

    reference.selected <- apply(refcounts[, refset$reference.choice, drop = FALSE], 1, sum)

    message("*** Creating ExomeDepth object...")
    suppressWarnings(ED <- new("ExomeDepth",
                               test = counts[, testsample],
                               reference = reference.selected,
                               formula = "cbind(test, reference) ~ 1"))

    message("*** Calling CNVs...")
    result <- CallCNVs(x = ED,
                       transition.probability = 10^-4,
                       chromosome = counts$chromosome,
                       start = counts$start,
                       end = counts$end,
                       name = counts$exon)

    if (length(result@CNV.calls) > 0) {
      cnvs <- result@CNV.calls[result@CNV.calls[, 9] > 0, ]
      if (nrow(cnvs) > 0) {
        samplecnvs <- rbind(samplecnvs, cnvs[, c(1, 2, 9, 10, 11, 12)])
      }
    }
  }

  cat("************\n", samplecnvs, "\n************\n")

  if (nrow(samplecnvs) > 0) {
    dups <- table(unlist(lapply(samplecnvs[samplecnvs$reads.ratio > 1, ], function(x) x[1]:x[2]))) / tries
    dels <- table(unlist(lapply(samplecnvs[samplecnvs$reads.ratio < 1, ], function(x) x[1]:x[2]))) / tries

    for (roi in unique(c(names(dups), names(dels)))) {
      exon.name <- exons$names[as.numeric(roi)]
      cnvs.ovp <- if (is.null(annotations)) NA else overlapsAny(exons[as.numeric(roi)], annotations)

      cat(sprintf("** aggregating %s\n", roi))

      if (!is.na(dups[roi])) {
        results <- rbind(results, data.frame(name = testsample, exon = exon.name, type = "duplication", prop = dups[roi], common = cnvs.ovp))
      }
      if (!is.na(dels[roi])) {
        results <- rbind(results, data.frame(name = testsample, exon = exon.name, type = "deletion", prop = dels[roi], common = cnvs.ovp))
      }
    }

    write.table(results, file = cnv.file, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE, append = TRUE)
    save.image(data.file)
  }
}