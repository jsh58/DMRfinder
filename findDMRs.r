#!/usr/bin/Rscript

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Dec. 2016

# Finding differentially methylated regions from a table
#   of counts produced by bisulfite sequencing data.
# Underlying statistics are performed using the R package 'DSS'
#   (https://bioconductor.org/packages/release/bioc/html/DSS.html).

version <- '0.1'
copyright <- 'Copyright (C) 2016 John M. Gaspar (jsh58@wildcats.unh.edu)'

printVersion <- function() {
  cat('findDMRs.r from DMRfinder, version', version, '\n')
  cat(copyright, '\n')
  q()
}

usage <- function() {
  cat('Usage: Rscript findDMRs.r  [options]  -i <input>  -o <output>  \\
             <groupList1>  <groupList2>  [...]
    -i <input>    File listing genomic regions and methylation counts
    -o <output>   Output file listing methylation results
    <groupList>   Comma-separated list of sample names (at least two
                    such lists must be provided)
  Options:
    -n <str>      Comma-separated list of group names
    -k <str>      Column names of <input> to copy to <output> (comma-
                    separated; def. "chr, start, end, CpG")
    -s <str>      Column names of DSS output to include in <output>
                    (comma-separated; def. "mu, diff, pval")
    -c <int>      Min. number of CpGs in a region (def. 3)
    -d <float>    Min. methylation difference between sample groups
                    ([0-1]; def. 0.10)
    -p <float>    Max. p-value ([0-1]; def. 0.05)
    -q <float>    Max. q-value ([0-1]; def. 1)
    -up           Report only regions hypermethylated in later group
    -down         Report only regions hypomethylated in later group
    -t <int>      Report regions with at least <int> comparisons
                    that are significant (def. 1)
')
  q()
}

# default args/parameters
infile <- outfile <- groups <- NULL
names <- list()        # list of sample names
minCpG <- 3            # min. number of CpGs
minDiff <- 0.10        # min. methylation difference
maxPval <- 0.05        # max. p-value
maxQval <- 1           # max. q-value (fdr)
up <- down <- F        # report NA/hyper-/hypo- methylated results
tCount <- 1            # min. number of significant comparisons
verbose <- F           # verbose option
keep <- c('chr', 'start', 'end', 'CpG')  # columns of input to keep
dss <- c('chr', 'pos', 'mu', 'diff', 'pval') # columns of DSS output to keep

# get CL args
args <- commandArgs(trailingOnly=T)
i <- 1
while (i <= length(args)) {
  if (substr(args[i], 1, 1) == '-') {
    if (args[i] == '-h' || args[i] == '--help') {
      usage()
    } else if (args[i] == '--version') {
      printVersion()
    } else if (args[i] == '-v') {
      verbose <- T
    } else if (args[i] == '-up') {
      up <- T
    } else if (args[i] == '-down') {
      down <- T
    } else if (i < length(args)) {
      if (args[i] == '-i') {
        infile <- args[i + 1]
      } else if (args[i] == '-o') {
        outfile <- args[i + 1]
      } else if (args[i] == '-n') {
        temp <- strsplit(args[i + 1], '[ ,]')[[1]]
        groups <- temp[nchar(temp) > 0]
      } else if (args[i] == '-k') {
        temp <- strsplit(args[i + 1], '[ ,]')[[1]]
        keep <- c(keep, temp[nchar(temp) > 0])
      } else if (args[i] == '-s') {
        temp <- strsplit(args[i + 1], '[ ,]')[[1]]
        dss <- c(dss, temp[nchar(temp) > 0])
      } else if (args[i] == '-c') {
        minCpG <- as.integer(args[i + 1])
      } else if (args[i] == '-d') {
        minDiff <- as.double(args[i + 1])
      } else if (args[i] == '-p') {
        maxPval <- as.double(args[i + 1])
      } else if (args[i] == '-q') {
        maxQval <- as.double(args[i + 1])
        dss <- c(dss, 'fdr')
      } else if (args[i] == '-t') {
        tCount <- as.integer(args[i + 1])
      } else {
        cat('Error! Unknown parameter:', args[i], '\n')
        usage()
      }
      i <- i + 1
    } else {
      cat('Error! Unknown parameter with no arg:', args[i], '\n')
      usage()
    }
  } else {
    temp <- strsplit(args[i], '[ ,]')[[1]]
    names <- c(names, list(temp[nchar(temp) > 0]))
  }
  i <- i + 1
}

# check for I/O errors, repeated columns
if (is.null(infile) || is.null(outfile)) {
  cat('Error! Must specify input and output files\n')
  usage()
}
input <- tryCatch( file(infile, 'r'), warning=NULL,
  error=function(e) stop('Cannot read input file ',
  infile, '\n', call.=F) )
output <- tryCatch( file(outfile, 'w'), warning=NULL,
  error=function(e) stop('Cannot write to output file ',
  outfile, '\n', call.=F) )
if (length(names) < 2) {
  cat('Error! Must specify at least two groups of samples\n')
  usage()
}
if (any(duplicated(unlist(names)))) {
  stop('Sample(s) repeated in different groups: ',
    paste(unique(unlist(names)[duplicated(unlist(names))]),
    collapse=', '), '\n')
}
keep <- unique(keep)
dss <- unique(dss)

# group samples into a named list
samples <- list()
for (i in 1:length(names)) {
  if (! is.null(groups) && i <= length(groups)) {
    group <- groups[i]
  } else {
    group <- paste(names[i][[1]], collapse='_')
  }
  if (group %in% names(samples)) {
    stop('Duplicated group name: ', group)
  }
  samples[[ group ]] <- names[i][[1]]
}

# load DSS
if (verbose) {
  cat('Loading DSS package\n')
}
suppressMessages(library(DSS))

# load data, check for errors
if (verbose) {
  cat('Loading methylation data from', infile, '\n')
}
data <- read.csv(input, sep='\t', header=T, check.names=F)
if (any( ! keep %in% colnames(data))) {
  stop('Missing column(s) in input file ', infile, ': ',
    paste(keep[! keep %in% colnames(data)], collapse=', '), '\n')
}
if (colnames(data)[1] != 'chr') {
  stop('Improperly formatted input file ', infile, ':\n',
    '  Must have "chr" as first column\n')
}

# determine columns for samples
idx <- list()
idx[[ 'N' ]] <- list()
idx[[ 'X' ]] <- list()
for (i in names(samples)) {
  idx[[ 'N' ]][[ i ]] <- rep(NA, length(samples[[ i ]]))
  idx[[ 'X' ]][[ i ]] <- rep(NA, length(samples[[ i ]]))
}
for (i in names(samples)) {
  for (j in 1:length(samples[[ i ]])) {
    for (k in 1:ncol(data)) {
      spl <- strsplit(colnames(data)[k], '-')[[1]]
      if (length(spl) < 2) { next }
      if (spl[-length(spl)] == samples[[ i ]][ j ]) {
        if (spl[length(spl)] == 'N') {
          idx[[ 'N' ]][[ i ]][ j ] <- k
        } else if (spl[length(spl)] == 'X') {
          idx[[ 'X' ]][[ i ]][ j ] <- k
        }
      }
    }
    if ( is.na(idx[[ 'N' ]][[ i ]][ j ])
        || is.na(idx[[ 'X' ]][[ i ]][ j ]) ) {
      stop('Missing information from input file ', infile, ':\n',
        '  For sample "', samples[[ i ]][ j ], '", need both "',
        samples[[ i ]][ j ], '-N" and "',
        samples[[ i ]][ j ], '-X" columns')
    }
  }
}

# for each sample, create data frames to meet DSS requirements
frames <- list()
for (i in names(samples)) {
  for (j in 1:length(samples[[ i ]])) {
    tab <- data.frame('chr'=data$chr, 'pos'=data$start,
      'N'=data[, idx[[ 'N' ]][[ i ]][ j ] ],
      'X'=data[, idx[[ 'X' ]][[ i ]][ j ] ])
    frames[[ samples[[ i ]][ j ] ]] <- tab
  }
}

# perform DML pairwise tests using DSS
bsdata <- makeBSseqData(frames, names(frames))
res <- data[, keep]  # results table
mat <- matrix(nrow=nrow(res), ncol=length(samples)*(length(samples)-1)/2)
  # matrix of booleans: does region meet threshold(s) for each comparison
comps <- c()  # group comparison strings
for (i in 1:(length(samples)-1)) {
  for (j in (i+1):length(samples)) {

    # perform DML test
    comp <- paste(names(samples)[i], names(samples)[j], sep='->')
    comps <- c(comps, comp)
    if (verbose) {
      cat('Comparing group "', names(samples)[i],
        '" to group "', names(samples)[j], '"\n  ', sep='')
    }
    if (length(samples[[ i ]]) < 2 || length(samples[[ j ]]) < 2) {
      # without replicates, must set equal.disp=T
      if (verbose) {
        dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]],
          equal.disp=T)
      } else {
        sink('/dev/null')
        dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]],
          equal.disp=T)
        sink()
      }
    } else {
      if (verbose) {
        dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]])
      } else {
        sink('/dev/null')
        dml <- DMLtest(bsdata, group1=samples[[ i ]], group2=samples[[ j ]])
        sink()
      }
    }

    # make sure necessary columns are present, remove extraneous
    col <- colnames(dml)
    if (any( ! dss %in% col & ! paste(dss, '1', sep='') %in% col)) {
      stop('Missing column(s) from DSS result: ',
        paste(dss[! dss %in% col & ! paste(dss, '1', sep='') %in% col],
        collapse=', '), '\n')
    }
    dml[, ! col %in% dss & ! substr(col, 1, nchar(col)-1) %in% dss] <- NULL

    # add results to res table
    start <- ncol(res) + 1
    res <- suppressWarnings( merge(res, dml,
      by.x=c('chr', 'start'), by.y=c('chr', 'pos'), all.x=T) )

    # determine if rows meet threshold(s)
    if (maxQval < 1) {
      mat[, length(comps)] <- ! (is.na(res[, 'diff'])
        | abs(res[, 'diff']) < minDiff
        | (up & res[, 'diff'] > 0) | (down & res[, 'diff'] < 0)
        | is.na(res[, 'pval']) | res[, 'pval'] > maxPval
        | is.na(res[, 'fdr']) | res[, 'fdr'] > maxQval)
    } else {
      mat[, length(comps)] <- ! (is.na(res[, 'diff'])
        | abs(res[, 'diff']) < minDiff
        | (up & res[, 'diff'] > 0) | (down & res[, 'diff'] < 0)
        | is.na(res[, 'pval']) | res[, 'pval'] > maxPval)
    }

    # add groups to column names
    for (k in start:ncol(res)) {
      col <- colnames(res)[k]
      if (substr(col, nchar(col), nchar(col)) == '1') {
        colnames(res)[k] <- paste(names(samples)[i],
          substr(col, 1, nchar(col)-1), sep=':')
      } else if (substr(col, nchar(col), nchar(col)) == '2') {
        colnames(res)[k] <- paste(names(samples)[j],
          substr(col, 1, nchar(col)-1), sep=':')
      } else {
        colnames(res)[k] <- paste(comp, col, sep=':')
      }
    }

  }
}

# filter regions based on CpG sites and mat matrix
if (verbose) {
  cat('Producing output file', outfile, '\n')
}
res <- res[res[, 'CpG'] >= minCpG & rowSums(mat) >= tCount, ]

# for repeated columns, average the values
repCols <- c()
sampleCols <- c()
groupCols <- c()
for (i in 1:ncol(res)) {
  if (i %in% repCols) { next }

  # find duplicated columns
  repNow <- c(i)
  if (i < ncol(res)) {
    for (j in (i + 1):ncol(res)) {
      if (colnames(res)[i] == colnames(res)[j]) {
        repNow <- c(repNow, j)
      }
    }
  }

  # average duplicates
  if (length(repNow) > 1) {
    res[, i] <- rowMeans(res[, repNow], na.rm=T)
    repCols <- c(repCols, repNow)
    sampleCols <- c(sampleCols, colnames(res)[i])
  } else if (! colnames(res)[i] %in% keep) {
    groupCols <- c(groupCols, colnames(res)[i])
  }
}
res <- res[, c(keep, sampleCols, groupCols)]

# limit results to 7 digits; reverse sign on diffs
options(scipen=999)
for (col in c(sampleCols, groupCols)) {
  spl <- strsplit(col, ':')[[1]]
  if (spl[length(spl)] == 'diff') {
    res[, col] <- -round(res[, col], digits=7)
  } else {
    res[, col] <- round(res[, col], digits=7)
  }
}

# sort chromosome names by number/letter
level <- levels(res$chr)
intChr <- strChr <- intLev <- strLev <- c()
for (i in 1:length(level)) {
  if (substr(level[i], 1, 3) == 'chr') {
    sub <- substr(level[i], 4, nchar(level[i]))
    if (!is.na(suppressWarnings(as.integer(sub)))) {
      intChr <- c(intChr, as.numeric(sub))
    } else {
      strChr <- c(strChr, sub)
    }
  } else {
    sub <- level[i]
    if (!is.na(suppressWarnings(as.integer(sub)))) {
      intLev <- c(intLev, as.numeric(sub))
    } else {
      strLev <- c(strLev, sub)
    }
  }
}
# put numeric chroms first, then strings
chrOrder <- c(paste('chr', levels(factor(intChr)), sep=''),
  levels(factor(intLev)),
  paste('chr', levels(factor(strChr)), sep=''),
  levels(factor(strLev)))

# write output results
write.table(res[order(match(res$chr, chrOrder), res$start), ],
  output, sep='\t', quote=F, row.names=F)
if (verbose) {
  cat('Regions reported:', nrow(res), '\n')
}
