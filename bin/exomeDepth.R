#!/usr/bin/env Rscript
# Rscript for ExomeDepth analysis with specific package versions

# Define packages and their required versions
packages <- list(
  GenomicRanges = "1.50.2",
  ExomeDepth = "1.1.16",
  xtable = "1.8.4",
  knitr = "1.45",
  kableExtra = "1.3.4",
  randomForest = "4.7.1.1",
  Rsamtools = "2.14.0",
  stringr = "1.5.1",
  dplyr = "1.1.4"
)

# Function to check and load packages
checkAndLoadPackage <- function(package, version) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(paste(package, "is not installed."))
  } else if (packageVersion(package) != version) {
    stop(paste("Version mismatch for", package, "- required:", version, "but installed:", packageVersion(package)))
  }
  library(package, character.only = TRUE)
}

# Apply the function to each package
lapply(names(packages), function(pkg) checkAndLoadPackage(pkg, packages[[pkg]]))

# Source the ExomeDepth script
source("/root/ed2vcf.R")

# Read arguments
cmd <- commandArgs(trailingOnly = FALSE)
runningScript <- unlist(strsplit(cmd[which(substr(cmd,1,7)=='--file=')], '='))[2]
scriptDirectory <- normalizePath(dirname(runningScript))
args <- commandArgs(trailingOnly = TRUE)

# Variables to configure pon and batch normnalisation plots
if (length(args) >= 7) {
  plot_pon <- args[6]
  plot_batch <- args[7]
} else {
  plot_pon <- "true"
  plot_batch <- "true"
}

setwd(dirname(args[2])) # set working directory to target directory

# Set other variables
panel <- unlist(strsplit(args[3], ':'))  # Covered regions name
testsample <- basename(unlist(strsplit(args[5], ':'))[1])  # RGID
samplename <- as.character(unlist(strsplit(args[5], ':'))[2]) # RGSM
threshold <- as.numeric(unlist(strsplit(args[5], ':'))[3]) # Threshold set in ngs_config

# Make sure ngs_config variables are set if missed out
if(is.na(threshold)) {threshold <- 10^-3}
if(is.na(plot_pon)) {plot_pon <- "true"}
if(is.na(plot_batch)) {plot_batch <- "true"}

# Version report header
argversion <- args[1]
versionfile <- paste(scriptDirectory, 'VERSION', sep = '/')
pipeversion <- ifelse(test = file.exists(versionfile),
                    paste(argversion, readChar(versionfile, 7), sep = "-"),
                    argversion)
message(paste('Running', pipeversion))

# Loads counts, samplenames, rois (extended exons.hg19)
load(args[4])

# Currently not implemented
extras <- NA

# Report run configuration
message(paste('        Version:', args[1]))
message(paste('         Output:', args[2]))
message(paste('            ROI:', panel[1]))
message(paste('     Panel Name:', panel[2]))
message(paste('    Read Counts:', args[4]))
message(paste('   Original BAM:', testsample))
message(paste('     SampleName:', samplename))
message(paste('  Normalisation:', names(refsamplenames)))
message(paste('         Extras:', extras))

# Check if testsample in refsamples
if (!testsample%in%names(refsets)) {
  stop(paste("The requested sample not available in", paste(names(refsets), collapse = ', ')))
}

results <- list()
qc <- data.frame() # Empty df to save qc table into

# Run if enough samples
if (length(names(refsamplenames)) > 0) {

  # Read exons/ROI and create subset
  message('Getting ROIs...')
  counts <- counts[which(counts$exon%in%rois$name), ] # Reduce to target regions (from readCount)
  exons <- GRanges(seqnames = rois$chromosome, IRanges(start = rois$start, end = rois$end), names = rois$name)
  covered <- with(read.table(panel[1], header = FALSE),
                  GRanges(seqnames = V1, IRanges(start = V2 + 1, end = V3, names = V4), names = V4)) # BED FILE 0-based
  coveredexons <- subsetByOverlaps(exons, covered)
  ce <- data.frame(seqnames = seqnames(coveredexons), starts = start(coveredexons)-1, ends = end(coveredexons), name = mcols(coveredexons)$names)
  ce$name <- na.omit(as.character(ce$name)) # Ensure 'ce$name' is a character vector and omit NAs
  selected.genes <- unique(unlist(lapply(strsplit(ce$name, '_'), function(x) x[1])))
  
  # Set defaults and load QC limits (override)
  limits <- list(
    medcor = c(NA, 0.90),    # Median correlation within batch
    maxcor = c(0.95, 0.90),  # Max correlation within batch
    refcor = c(0.95, 0.90),  # Reference set correlation
    refcount = c(3, 1),      # Refernce set size (selected reference samples)
    coeffvar = c(NA, NA),    # Coefficient of variation
    cvzscore = c(1, 4),      # Z-Score of the coefficient of variation
    coverage = c(100),       # Minimum exon depth (read count)
    expectedbf = c(5.0, NA), # Expected BF
    minrefs = c(2, Inf)      # Minimum reference set size
  )
  annotations <- NA  # CNV annotations
  
  # Load from file (including annotations)
  if (!is.na(extras)) {
    extras <- ifelse(startsWith(extras, '/'), extras, paste(scriptDirectory, extras, sep = '/'))
    load(extras)
  }

  # Define low coverage exons
  message('Identifying low coverage exons...')
  exonnames <- coveredexons@elementMetadata@listData$names
  coverage.df <- data.frame()
  for (rs in names(refsamplenames)) {
    coverage.df <- rbind(coverage.df, data.frame(
      refset = rs, 
      exon = counts$exon, 
      gc = counts$GC, 
      coverage.min = apply(counts[, refsamplenames[[rs]]], 1, min),
      coverage.median = apply(counts[, refsamplenames[[rs]]], 1, median),
      coverage.max = apply(counts[, refsamplenames[[rs]]], 1, max)
    ))
  }
  coverage.table <- coverage.df[which(coverage.df$coverage.median<limits$coverage & coverage.df$exon%in%exonnames), ]

  # Check if reference sample choice has mismatched sex
  message('Checking sex match...')
  sexMismatch <- list()
  for (rs in names(refsamplenames)) {
	sexMismatch[[rs]] <- rep(NA, length(refsets[[testsample]][[rs]]$reference.choice))
    if (any(rois[, 1] == "X") || any(rois[, 1] == "Y")) {
      sex <- regex('_[FM]_')
      # Get sex of testsample (ignore undefined)
      tssx <- str_extract(testsample, sex)
      rcsx <- str_extract(refsets[[testsample]][[rs]]$reference.choice, sex)
      # If sex information available
      if (!is.na(tssx) && !all(is.na(rcsx))) {
        sexMismatch[[rs]] <- !(rcsx == tssx)
      }
    }
  }

  # Prepare reference (sum reference choice)
  message('Preparing reference set...')
  reference.selected <- list()
  ref.correlation <- list()
  for (rs in names(refsamplenames)) {
    reference.selected[[rs]] <- as.vector(apply(
      X = as.matrix(counts[, refsets[[testsample]][[rs]]$reference.choice, drop = FALSE]),
      MAR = 1,
      FUN = sum
    ))
    ref.correlation[[rs]] <- cor(cbind(rpkm[, testsample], calcRPKM(reference.selected[[rs]], (counts$end-counts$start+1))))[1, 2]
  }

  # Build QC table
  message('Building QC summary...')
  decide <- function(v, t) {
    result <- list(
      value = round(v, 3), 
      threshold = 'No threshold set',
      warning = NA,
      fail = NA,
      status = "PASS"
    )
    if (!is.null(t)) {
      cmp <- ifelse(any(is.na(t)) || t[1]>t[2], function(m, n) m >= n, function(m, n) m < n)
      threshold <- ifelse(any(is.na(t)) || t[1]>t[2], "equal or greater than", "less than" )
      status <- ifelse(!is.na(t[2]) && !cmp(v, t[2]), "FAIL", 
                     ifelse(!is.na(t[1]) && !cmp(v, t[1]), "CAUTION", "PASS"))
      result <- list(
        value = round(v, 3),
        threshold = threshold,
        warning = t[1],
        fail = t[2],
        status = status
      )
    }
    result
  }
  
  save.image()
  message('•••SAVED•••')

  # Create a QC table for each refset
  qc <- list()
  for (rs in unique(stats$refset)) {
    stats.refset <- stats[which(stats$refset ==  rs), 
                        which(colnames(stats) != "refset")]

    # Calc CV Z-score for each batch
    batches <- as.factor(unlist(lapply(strsplit(rownames(stats.refset), '_'), '[[', 1)))
    batches.mean <- by(stats.refset$coeff.var, batches, mean)
    batches.sd <- by(stats.refset$coeff.var, batches, sd)

    # Get/calc testsample data
    testsample.cv <- stats.refset[which(stats.refset$sample == testsample), 'coeff.var']
    testsample.batch <- batches[which(stats.refset$sample == testsample)]
    batch.mean <- batches.mean[[testsample.batch]]
    batch.sd <- batches.sd[[testsample.batch]]
    testsample.cvz <- (testsample.cv - batch.mean)/batch.sd

    qc[[rs]] <- rbind(
      decide(stats.refset[testsample, "batch.mediancor"], limits$medcor), 
      decide(stats.refset[testsample, "batch.maxcor"], limits$maxcor), 
      decide(stats.refset[testsample, "coeff.var"], limits$coeffvar), 
      decide(testsample.cvz, limits$cvzscore), 
      decide(ref.correlation[[rs]], limits$refcor), 
      decide(stats.refset[testsample, "refsamples"], limits$refcount), 
      decide(stats.refset[testsample, "min.refs"], limits$minrefs), 
      decide(stats.refset[testsample, "expected.BF"], limits$expectedbf)
    )
    rownames(qc[[rs]]) = c(
      "Median correlation in batch",
      "Maximum correlation in batch",
      "Coefficient of variation (CV)",
      "CV Z-score within batch",
      "Correlation with reference",
      "Size of reference set",
      "Forced minimum reference set size",
      "Expected BF with reference"
    )
    qc[[rs]] <- as.data.frame(qc[[rs]])
    print(qc[[rs]])
  }

  # Run ExomeDepth
  cnvs_all <- list()
  for (rs in names(reference.selected)) {
    message('*** Creating ExomeDepth object...')
    suppressWarnings(ED <- new('ExomeDepth',
                               test = counts[, testsample],
                               reference = reference.selected[[rs]],
                               formula = 'cbind(test, reference) ~ 1'))

    # Call CNV
    message('*** Calling CNVs...')
    cnvs <- CallCNVs(x = ED,
                   transition.probability = threshold,
                   chromosome = counts$chromosome,
                   start = counts$start,
                   end = counts$end,
                   name = counts$exon)
    print(paste('Raw CNV count:', length(cnvs@CNV.calls)))

    # Annotate results
    message("Annotating CNVs...")
    exon_overlap_frac <- 0.000000001  # 1bp/1Gb, basically any overlap
    extra_overlap_frac <- 0.1  # 1bp/10bp (only valid for if annotation is large known SegDups)
    if (length(cnvs@CNV.calls)>0) {
      # Add exon numbers (from subset)
      cnvs.annotated <- AnnotateExtra(x = cnvs, 
                                      reference.annotation = coveredexons,
                                      min.overlap = exon_overlap_frac,
                                      column.name = "exons.hg19")

      # Add extra annotation
      if (all(is.na(annotations))) {
        message("No CNV annotaions available")
      } else {
        message("Adding extra annotations...")
        cnvs.annotated <- AnnotateExtra(x = cnvs.annotated, 
                                        reference.annotation = annotations,
                                        min.overlap = extra_overlap_frac,
                                        column.name = "annotation")
      }
      cnvs.annotated@annotations$name <- as.factor(
        sapply(strsplit(
          as.character(cnvs.annotated@annotations$name), '_'), "[[", 1))
      results[[rs]] <- cnvs.annotated
    } else {
      message(paste('No CNVs called with refset', toupper(rs)))
      results[[rs]] <- NA  # No CNVs called
    }

    # Write BED file of CNVs
    message(paste("Writing BED file for refset", toupper(rs)))
    getBed <- function(x) {
      if (class(x) ==  "ExomeDepth") {
        beddata <- x@CNV.calls[which(x@CNV.calls[, 9] > 0), c(7, 5, 6, 3, 9)]
        beddata[, 2] <- beddata[, 2] - 1
        return(beddata)
      }
      NULL
    }
    bed.filename <- sub("[.][^.]*$", paste0(".", rs, ".bed"), args[2],
                        perl = TRUE)
    write.table(
      getBed(results[[rs]]), 
      file = bed.filename, sep = "\t", 
      col.names = FALSE, quote = FALSE, row.names = FALSE)

    # Write VCF file
    message(paste("Writing VCF file for refset", toupper(rs)))
    if (file.exists(referenceFasta)) {
      vcf.filename <- sub("[.][^.]*$", paste0(".", rs, ".vcf"), args[2],
                          perl = TRUE)
      ed2vcf(results[[rs]], vcf.filename, referenceFasta, rois, samplename)
    }
  }
} else {
  message("Skipping ExomeDepth. No suitable reference set available")
}

if (plot_pon == "false") {
  results["pon"] <- NA
}

if (plot_batch == "false") {
  results["batch"] <- NA
}

# Knit report (using refsets, results)
knitrScript <- paste(scriptDirectory, "exomeDepth.Rnw", sep = "/")
if (sub(".*[.]", "", args[2], perl = TRUE) == "pdf") {
  knit2pdf(knitrScript, output = sub("[.][^.]*$", ".tex", args[2], perl = TRUE))
}

save.image(sub("[.][^.]*$", ".RData", args[2], perl = TRUE))