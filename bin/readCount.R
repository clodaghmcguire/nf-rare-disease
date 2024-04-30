#!/usr/bin/env Rscript

# Rscript for exomedepth data

#### Install packages ####

# BiocManager and ExomeDepth
options(menu.graphics = FALSE)

# Required packages with specified versions
required_packages <- list(
  plyr = "1.8.9",
  stringr = "1.5.1",
  ExomeDepth = "1.1.16",
  Biostrings = "2.64.1",
  IRanges = "2.30.1",
  Rsamtools = "2.14.0",
  GenomicRanges = "1.50.2",
  aod = "1.3.2",
  VGAM = "1.1.9",
  GenomicAlignments = "1.32.1",
  dplyr = "1.1.4"
)

installed_pkgs <- installed.packages()[, "Package"]
installed_versions <- installed.packages()[, "Version"]

for (pkg in names(required_packages)) {
  if (pkg %in% installed_pkgs && installed_versions[installed_pkgs == pkg] == required_packages[[pkg]]) {
    library(pkg, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    cat(paste("Loaded", pkg, "version", required_packages[[pkg]], "\n"))
  } else {
    cat(paste("Required version of", pkg, "(", required_packages[[pkg]], ") is not installed.\n"))
  }
}

#### Setup ####

# Configs #
normalprefix <- 'NORMAL'
cmp.cols <- 3
sex <- regex('_[MF]_')
min.refs <- 2
mincov.pon <- 0

# Read stdin arguments
args <- commandArgs(trailingOnly = TRUE)

# Print what will be printed
message(paste('     Output:', args[1]))
message(paste('   Reference:', args[2]))
referenceFasta <- args[2]

# Read ROI
targetfiles <- vector()
for (roi in unlist(strsplit(args[3], ','))) {
  targetfiles <- append(targetfiles, roi)
  message(paste('      ROI:', roi))
}

# Load/combine exon/ROI data and make unique
rois <- NULL
for (tf in targetfiles) rois <- rbind(rois, read.table(tf, header = FALSE)[, 1:4])
rois <- unique(rois)
colnames(rois) <- c('chromosome', 'start', 'end', 'name')

# Read tracks (BAM files, precomputed normals)
bam <- vector()
pcn <- vector()
for (bamfile in args[4:length(args)]) {
  if (endsWith(unlist(strsplit(bamfile, ':'))[1], '.bam')) {
    bam <- append(bam, unlist(strsplit(bamfile, ':'))[1])
    message(paste('   Alignments:', bamfile))
  } else if (endsWith(bamfile, '.RData')) {
    pcn <- append(pcn, bamfile)
    message(paste('     Counts:', bamfile))
  } else {
    message(paste('    Skipped:', bamfile))
  }
}

# Omit normals that are in supplied counts file
names.imported <- vector()
for (c in pcn) {
  attach(c)
  names.imported <- append(names.imported,
  colnames(counts)[which(startsWith(colnames(counts), normalprefix))])
  detach()
}

# Remove imported normals from bam list
bams.to.count <- which(!as.vector(sapply(bam, basename))%in%names.imported)
bam <- bam[bams.to.count]

# Generate read counts from BAM files for all exons
counts.computed <- NA
if (length(bam) > 0) {
    suppressWarnings(counts.computed <- getBamCounts(
        bed.frame = rois,
        bam.files = bam,
        min.mapq = 10,
        include.chr = FALSE,  # Chrom start with chr prefix
        referenceFasta = referenceFasta))
    counts.computed <- as(counts.computed, 'data.frame')
}

for(column in 6:length(names(counts.computed))){
  if(nchar(names(counts.computed)[column])==20){
    names(counts.computed)[column] <- sub('X', '', names(counts.computed)[column])
  }
}

# Add precomputed counts (normals only)
for (c in pcn) {
  # Import precomputed normals
  attach(c)
  targets.imported <- counts[, c(1:cmp.cols)]
  counts.imported <- counts
  detach()
  # If compatible normals, amend counts table
  if (length(bam)>0) {
    targets.same <- isTRUE(all.equal(
      counts.computed[, c(1:cmp.cols)], targets.imported))
    if (targets.same && ncol(counts.imported) > 0) {
      counts.computed <- cbind(counts.computed, counts.imported)
    } else {
      message(paste('ERROR: Supplied count data not compatible', c))
    }
  } else {
    message('No computed read depth data, using supplied counts')
    # If no computed counts -> just use imported from RData file
    counts.computed <- counts.imported
  }
}

counts <- counts.computed

# Select normal samples (if 3+ specified)
bams <- colnames(counts)[which(endsWith(colnames(counts), '.bam'))]
normals <- which(substr(bams, 1, nchar(normalprefix)) == normalprefix)
tests <- which(substr(bams, 1, nchar(normalprefix)) != normalprefix)

# Test samples are same regardless of normalisation method
testsamplenames <- bams[tests]

# Get refsamplenames
refsamplenames <- list()
if (length(tests) >= 3) {
    message('-> Batch normalisation')
    refsamplenames[['batch']] <- bams[tests]
}

if (length(normals) >= 2) {
  message(paste("-> Panel of Normals [", length(normals), "]..."))
    refsamplenames[['pon']] <- bams[normals]
}

# Create PoN if no test samples
if (length(testsamplenames) == 0) {
  message('Creating Panel of Normals (PoN)...')

    # Remove all samples from PoN that do not match coverage requirement
    bad.normals <- names(which(apply(counts[, bams], 2, min)<mincov.pon))
    if (length(bad.normals)>0) {
        message("Removing normals that have insufficient ROI coverage")
        print(bad.normals)
        counts <- counts[, -which(colnames(counts) %in% bad.normals)]
    }

  # Save normal counts
  save(list = c("refsamplenames",  # supplied normals
        "counts"),     # normal counts
     file = args[1])
    quit(status = 0, save = 'no')
}

# Find reference sets
message('Building reference sets...')

# Pick reference sample set from batch and pon
selectReferenceSet <- function(testsample) {
    selected.refs <- list()
    hasXY <- any(rois[,1] == "X") || any(rois[,1] == "Y")
    for (rs in names(refsamplenames)) {
        message(paste('Picking', toupper(rs), 'reference for', testsample))

        # exclude test sample
        refsamples <- refsamplenames[[rs]][which(
          refsamplenames[[rs]] != testsample)]

        # if PoN && XY && Sex -> sex match
        hasPon <- rs == 'pon'
        tssx <- str_extract(testsample, sex)
        rcsx <- str_extract(refsamples, sex)
        hasSex <- !(is.na(tssx) || any(is.na(rcsx)))
        if (hasPon && hasXY && hasSex) {
            message(paste('Enforcing sex match of', testsample, 'to reference'))
            refsamples <- refsamples[which(rcsx == tssx)]
        } else if (!hasPon && hasXY && !is.na(tssx)) {
          # If no PoN, XY targets and testsample has sex, do sex
          # match if at least 2 matched references exist within batch
            refsamples_candidates <- refsamples[which(rcsx == tssx)]
            if (length(refsamples_candidates)>2) {
                message(paste(
                  'Enforcing sex match of', testsample, 'within batch'
                  ))
                refsamples <- refsamples_candidates
            }
        }

        # Pick reference 
        suppressWarnings(selected <- select.reference.set(
            test.counts = counts[, testsample], 
            reference.counts = as.matrix(counts[, refsamples]), 
            bin.length = (counts$end - counts$start)
        ))

        # If not enough reference samples, enlarge reference set
        ref.count <- which(selected$summary.stats$selected)
        if (ref.count < min.refs) {
            # Create new df
            new.summary.stats <- selected$summary.stats
            new.summary.stats$selected <- FALSE
            # Get row with highest expected.BF after min.refs rows
            BF.ranks <- order(new.summary.stats$expected.BF, decreasing = T, na.last = T)
            best.BF.rank <- which.min(BF.ranks[min.refs:length(BF.ranks)])+min.refs-1
            # Mark selected row 
            new.summary.stats$selected[best.BF.rank] <- TRUE
            # Build new selected reference
            new.selected <- list(
                reference.choice = as.character(new.summary.stats$ref.samples[1:best.BF.rank]), 
                summary.stats = new.summary.stats, 
                min.refs = min.refs)
            message(paste('Enlarged reference set:', ref.count, ' = >', best.BF.rank))
            selected <- new.selected
        } else {
            message(paste('Reference set size:', ref.count))
            selected$min.refs <- 1
        }

        # Save reference set to list
        selected.refs[[rs]] <- selected
    }
    selected.refs
}

refsets <- lapply(testsamplenames, selectReferenceSet)
names(refsets) <- testsamplenames

#### Statistics output ####

# Calculate RPKM
calcRPKM <- function(c, l) c/(l*sum(c)/10^6)
counts.len <- counts$end-counts$start+1
rpkm <- apply(counts[, testsamplenames], 2, function(x) calcRPKM(x, counts.len))

# Calculate batch statistics (all samples excluding normals if any)
batch.cv <- apply(rpkm, 2, function(r) sd(r)/mean(r)*100)
batch.cor <- cor(rpkm)
diag(batch.cor) <- NA
stats <- data.frame()

for (rs in names(refsamplenames)) {
    for (testsample in testsamplenames) {
        message(paste('Collecting stats for', testsample, 'with', toupper(rs), 'reference set'))
        d <- cbind(
                 refset = rs, 
                 sample = testsample, 
                 min.refs = refsets[[testsample]][[rs]]$min.refs, 
                 refsamples = which(refsets[[testsample]][[rs]]$summary.stats$selected), 
                 refsets[[testsample]][[rs]]$summary.stats[which(refsets[[testsample]][[rs]]$summary.stats$selected), 
                                        c('correlations', 'expected.BF', 'phi', 'RatioSd', 'mean.p', 'median.depth')], 
                 batch.maxcor = max(batch.cor[testsample, ], na.rm = TRUE), 
                 batch.mediancor = median(batch.cor[testsample, ], na.rm = TRUE), 
                 coeff.var = batch.cv[testsample]
                )
        rownames(d) <- testsample
        stats <- rbind(stats, d)
    }
}

print(stats)

#### Output ####

# Save table
write.table(stats, file = sub("[.][^.]*$", ".csv", args[1], perl = TRUE), sep = '\t', quote = FALSE, row.names = FALSE)

# Save read count table as .RData file
save(list = c(
  "referenceFasta",   # reference fasta
  "refsamplenames",   # selected reference samples (normals if provided)
  "testsamplenames",  # supplied test samples
  "counts",           # all read counts
  "rois",             # regions of interest (full)
  "batch.cv",         # Coefficient of Variation for each sample in batch
  "batch.cor",        # Correlation within batch
  "refsets",          # picked reference sets
  "rpkm",             # RPKM calculation
  "calcRPKM",         # function to calculate RPKM
  "stats"),           # model and batch statistics
  file = args[1])

message('DONE')
