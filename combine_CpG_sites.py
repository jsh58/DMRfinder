#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Dec. 2016

# Combining multiple samples' single-base resolution
#   methylation data into a set of genomic regions.

import sys
import os.path
import gzip
import math
version = '0.4_dev'
copyright = 'Copyright (C) 2016 John M. Gaspar (jsh58@wildcats.unh.edu)'

def printVersion():
  sys.stderr.write('combine_CpG_sites.py from DMRfinder, version' \
    + ' %s\n' % version)
  sys.stderr.write(copyright + '\n')
  sys.exit(-1)

def usage():
  sys.stderr.write('''Usage: python combine_CpG_sites.py  [options]  -o <output>  [<input>]+
    [<input>]+    One or more files, each listing methylation counts
                    for a particular sample
    -o <output>   Output file listing genomic regions and combined
                    methylation counts for each sample
  Options:
    To consider a particular CpG:
      -r <int>    Min. number of counts at a position (def. 3)
      -s <int>    Min. number of samples with -r counts (def. 1)
    To analyze a region of CpGs:
      -d <int>    Max. distance between CpG sites (def. 100)
      -c <int>    Min. number of CpGs in a region (def. 3)
      -x <int>    Max. length of a region (def. 500)
    To report a particular result:
      -m <int>    Min. total counts in a region (def. 20)
    Other:
      -f          Report methylation fraction for each sample
      -y <file>   Report clusters of valid CpG sites;
                    if <file> exists, use these clusters
      -b          Analyze one chromosome at a time (memory-saving)
      -e <file>   File listing ordered chromosome names (comma-
                    separated; used only with -b option)
''')
  sys.exit(-1)

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot read input file %s\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
  '''
  if filename == '-':
    return sys.stdout
  try:
    f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot write to output file %s\n' % filename)
    sys.exit(-1)
  return f

def getInt(arg):
  '''
  Convert given argument to int.
  '''
  try:
    val = int(arg)
  except ValueError:
    try:
      val = int(float(arg))
    except ValueError:
      sys.stderr.write('Error! Cannot convert %s to int\n' % arg)
      sys.exit(-1)
  return val

def openFiles(infiles):
  '''
  Open files and save file names.
  '''
  files = []
  samples = []
  for infile in infiles:
    files.append(openRead(infile))

    # save sample name (basename of file)
    sample = infile.split('/')[-1].split('.')[0]
    while sample in samples:
      sample += '-'
    if sample[0] == '-':
      sample = '_' + sample[1:]  # change leading '-'
    samples.append(sample)

  return files, samples

def writeHeader(fOut, samples, fraction):
  '''
  Write the header for the output file.
  '''
  fOut.write('\t'.join(['chr', 'start', 'end', 'CpG']))
  if fraction:
    fOut.write('\t' + '\t'.join(samples))
  else:
    for sample in samples:
      for letter in ['N', 'X']:
        fOut.write('\t' + sample + '-' + letter)
  fOut.write('\n')

def loadChrOrder(infiles, files):
  '''
  Return a list of ordered chromosome names,
    constructed from input count files.
  '''
  # load chroms from input files
  order = [[] for i in range(len(files))]
  for i in range(len(files)):
    refChrom = ''
    for line in files[i]:
      chrom = line.split('\t')[0]
      if chrom != refChrom:
        if chrom in order[i]:
          sys.stderr.write('Error! Unsorted input file:' \
            + ' %s\n' % infiles[i])
          sys.exit(-1)
        order[i].append(chrom)
        refChrom = chrom
    files[i].seek(0)  # rewind file

  # create master list (chrOrder)
  #   (not optimal, but this is an SCS problem)
  order.sort(key=len, reverse=True)
  chrOrder = order[0]
  for i in range(1, len(order)):

    prevIdx = -1
    j = 0
    while j < len(order[i]):

      # find missing samples, insert into chrOrder
      miss = []
      while order[i][j] not in chrOrder:
        miss.append(order[i][j])
        j += 1
      if miss:
        prevIdx = chrOrder.index(order[i][j])
        for m in miss[::-1]:
          chrOrder.insert(prevIdx, m)
        continue

      # make sure orders match
      idx = chrOrder.index(order[i][j])
      if idx < prevIdx:
        sys.stderr.write('Error! Input files not in same order\n')
        sys.stderr.write('  (if you think this is incorrect, please specify\n')
        sys.stderr.write('  the chromosome order via the -e argument)\n')
        sys.exit(-1)

      prevIdx = idx
      j += 1

  return chrOrder

def saveChrOrder(chrOrderFile, infiles, files, verbose):
  '''
  Determine order of chromosomes to analyze.
  '''
  chrOrder = []  # order of chromosomes to process
  if chrOrderFile != None:
    # load chromosome order from given file
    f = openRead(chrOrderFile)
    for line in f:
      chrOrder.extend(line.rstrip().split(','))
    if f != sys.stdin:
      f.close()

  else:
    # need to construct order from input count files
    if verbose:
      sys.stderr.write('Constructing chromosome order\n')
    chrOrder = loadChrOrder(infiles, files)

  return chrOrder

def splitRegion(chrom, reg, count, minCpG, minReg, \
    maxLen, samples, fraction, fOut, fClus):
  '''
  Split a CpG region that is too large and process
    each subregion via processRegion().
  '''
  # determine number of subregions and length
  subReg = math.ceil((reg[-1] - reg[0]) / float(maxLen))
  while len(reg) / subReg < minCpG:
    subReg -= 1  # too few CpGs: decrease number
  lengthReg = (reg[-1] - reg[0]) / subReg
  subReg = int(subReg)

  # create subregions based on length
  start = 0     # index of beginning of subregion
  prev = reg[0] # genomic position of beginning of subregion
  ends = []
  for j in range(len(reg)):
    if reg[j] > prev + lengthReg and j - start >= minCpG:
      ends.append(j)
      if len(ends) == subReg - 1: break
      start = j
      prev += lengthReg
  while len(ends) < subReg:
    ends.append(len(reg))

  # make sure each region has at least minCpG
  j = len(ends) - 1
  while j and ends[j] - ends[j - 1] < minCpG:
    ends[j - 1] = ends[j] - minCpG
    j -= 1

  # process subregions
  start = 0     # index of beginning of subregion
  total = 0     # number of regions printed
  for end in ends:
    # pass to processRegion()
    total += processRegion(chrom, reg[start:end], count, \
      minCpG, minReg, float('inf'), samples, fraction, \
      fOut, fClus)
    start = end

  return total

def processRegion(chrom, reg, count, minCpG, minReg, \
    maxLen, samples, fraction, fOut, fClus):
  '''
  Produce output for a given region of CpGs: a line
    containing chromosome name, start and end
    coordinates, number of CpGs, and methylation
    data for each sample, all tab-delimited.
  To print a line, the region must have at least
    <minCpG> sites, and at least one sample must have
    at least <minReg> counts.
  Any region longer than <maxLen> will be split via
    splitRegion(), as long as <minCpG> is still
    maintained by the subregions.
  Any sample that does not have <minReg> counts gets
    an 'NA' designation.
  '''
  if len(reg) < minCpG:
    return 0

  # split region larger than maxLen
  if reg[-1] - reg[0] > maxLen:
    return splitRegion(chrom, reg, count, minCpG, minReg, \
      maxLen, samples, fraction, fOut, fClus)

  flag = False  # boolean for printing line
  res = '%s\t%d\t%d\t%d' % (chrom, reg[0], reg[-1], len(reg))
  for sample in samples:
    meth = unmeth = 0

    # sum methylated/unmeth bases at each position in region
    for r in reg:
      pos = str(r)
      if sample in count[pos]:
        meth += count[pos][sample][0]
        unmeth += count[pos][sample][1]

    if meth + unmeth < minReg:
      # less than minimum number of counts
      res += '\tNA'
      if not fraction:
        res += '\tNA'
    else:
      if fraction:
        # compute methylated fraction
        res += '\t%f' % (meth / float(meth + unmeth))
      else:
        res += '\t%d\t%d' % (meth + unmeth, meth)  # actual counts
      flag = True  # at least 1 valid sample

  if flag:
    fOut.write(res + '\n')

    # record CpG sites in cluster
    if fClus:
      fClus.write(chrom + '\t' + str(reg[0]))
      for r in reg[1:]:
        fClus.write(',' + str(r))
      fClus.write('\n')
    return 1

  return 0

def combineRegions(count, total, chrom, minSamples, maxDist, \
    minCpG, minReg, maxLen, samples, fraction, fOut, fClus):
  '''
  Combine data from CpG positions that are close to each
    other (a modified single-linkage clustering, with
    distance parameter maxDist). Process combined regions
    on the fly (via processRegion() function).
  '''
  printed = 0  # count of printed regions
  reg = []  # for saving connected positions
  pos3 = 0
  for pos in sorted(total, key=int):

    # require a min. number of samples
    if total[pos] >= minSamples:
      loc = getInt(pos)

      # if next position is more than maxDist away,
      #   process previous genomic region
      if pos3 and loc - pos3 > maxDist:
        printed += processRegion(chrom, reg, count, minCpG, \
          minReg, maxLen, samples, fraction, fOut, fClus)
        reg = []  # reset list
      reg.append(loc)
      pos3 = loc

  # process last genomic region for this chromosome
  printed += processRegion(chrom, reg, count, minCpG, \
    minReg, maxLen, samples, fraction, fOut, fClus)
  return printed

def loadCounts(f, minReads, count, total, order, sample):
  '''
  Load the methylated/unmethylated counts for a file.
  '''
  # load counts from file
  for line in f:
    try:
      chrom, pos, end, pct, meth, unmeth \
        = line.rstrip().split('\t')
    except ValueError:
      sys.stderr.write('Error! Poorly formatted record:\n%s' % line)
      sys.exit(-1)
    meth = getInt(meth)
    unmeth = getInt(unmeth)

    # save counts and total
    if chrom not in count:
      count[chrom] = {}
      total[chrom] = {}
      order.append(chrom)
    if pos not in count[chrom]:
      count[chrom][pos] = {}
    count[chrom][pos][sample] = [meth, unmeth]
    # save to 'total' dict. only if sufficient coverage
    if meth + unmeth >= minReads:
      total[chrom][pos] = total[chrom].get(pos, 0) + 1

def processFiles(infiles, files, samples, minReads, \
    minSamples, maxDist, minCpG, minReg, maxLen, fraction, \
    fOut, fClus, verbose):
  '''
  Load all methylation counts from the inputs,
    cluster, and output results.
  '''
  # load methylation information for each sample
  count = {}    # for methylated, unmethylated counts
  total = {}    # for number of samples with min. coverage
  order = []    # for ordered chromosome names
  if verbose:
    sys.stderr.write('Loading methylation information\n')
  for i in range(len(infiles)):
    if verbose:
      sys.stderr.write('  file: %s\n' % infiles[i])
    loadCounts(files[i], minReads, count, total, order, samples[i])

  # cluster and produce output
  if verbose:
    sys.stderr.write('Combining regions and producing output\n')
  printed = 0
  for chrom in order:
    printed += combineRegions(count[chrom], total[chrom], \
      chrom, minSamples, maxDist, minCpG, minReg, maxLen, \
      samples, fraction, fOut, fClus)
  return printed

def loadChromCounts(files, samples, lines, refChrom, \
    minReads, count, total):
  '''
  Load the methylated/unmethylated counts for a
    chromosome from all the input files.
  '''
  # load counts from each file
  for i in range(len(files)):
    if not lines[i]: continue
    try:
      chrom, pos, end, pct, meth, unmeth \
        = lines[i].rstrip().split('\t')
    except ValueError:
      sys.stderr.write('Error! Poorly formatted record:' \
        + '\n%s' % lines[i])
      sys.exit(-1)
    if chrom != refChrom: continue

    # load counts as long as chrom matches
    while chrom == refChrom:
      meth = getInt(meth)
      unmeth = getInt(unmeth)

      if pos not in count:
        count[pos] = {}
      count[pos][samples[i]] = [meth, unmeth]
      # save to 'total' dict. only if sufficient coverage
      if meth + unmeth >= minReads:
        total[pos] = total.get(pos, 0) + 1

      line = files[i].readline()
      if not line:
        break
      try:
        chrom, pos, end, pct, meth, unmeth \
          = line.rstrip().split('\t')
      except ValueError:
        sys.stderr.write('Error! Poorly formatted record:' \
          + '\n%s' % line)
        sys.exit(-1)

    # reset next line
    lines[i] = line

def processChrom(infiles, files, samples, minReads, \
    minSamples, maxDist, minCpG, minReg, maxLen, \
    fraction, fOut, chrOrder, fClus, verbose):
  '''
  Process the input files, one chromosome at a time.
  '''
  # load first lines of files
  lines = []
  for f in files:
    lines.append(f.readline())

  # process each chromosome separately
  if verbose:
    sys.stderr.write('Loading methylation information,\n')
    sys.stderr.write('  combining regions, and producing output\n')
  printed = 0
  for refChrom in chrOrder:
    if verbose:
      sys.stderr.write('  chromosome: %s\n' % refChrom)

    # load counts from each file
    count = {}    # for methylated, unmethylated counts
    total = {}    # for number of samples with min. coverage
    loadChromCounts(files, samples, lines, refChrom, \
      minReads, count, total)

    # cluster and print output
    printed += combineRegions(count, total, refChrom, \
      minSamples, maxDist, minCpG, minReg, maxLen, \
      samples, fraction, fOut, fClus)

  # check for unprocessed records
  unProc = False
  for i in range(len(lines)):
    if lines[i]:
      if not unProc:
        sys.stderr.write('Error! Unprocessed records in inputs:\n')
      sys.stderr.write('File %s:\n%s' % (infiles[i], lines[i]))
      unProc = True
  if unProc:
    sys.exit(-1)

  return printed

def loadClusters(clusFile, verbose):
  '''
  Load info for all clusters from external file.
  '''
  if verbose:
    sys.stderr.write('Loading cluster information\n')
  f = openRead(clusFile)
  total = sites = 0
  idx = dict()
  clus = list()
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 2:
      sys.stderr.write('Error! Poorly formatted cluster record:' \
        + '\n%s' % line)
      sys.exit(-1)
    pos = spl[1].split(',')

    # save sites and cluster header
    for p in pos:
      idx[spl[0] + ' ' + p] = total
      sites += 1
    clus.append('\t'.join([spl[0], pos[0], pos[-1], \
      str(len(pos))]))
    total += 1

  if f != sys.stdin:
    f.close()
  if verbose:
    sys.stderr.write('  Clusters loaded: %d\n' % total)
    sys.stderr.write('  CpG sites: %d\n' % sites)
  return clus, idx

def loadCountsClus(f, count, idx, sample):
  '''
  Load the methylated/unmethylated counts for a file.
  '''
  # load counts from file
  for line in f:
    try:
      chrom, pos, end, pct, meth, unmeth \
        = line.rstrip().split('\t')
    except ValueError:
      sys.stderr.write('Error! Poorly formatted record:\n%s' % line)
      sys.exit(-1)

    # save counts
    k = chrom + ' ' + pos
    if k in idx:
      if sample in count[idx[k]]:
        count[idx[k]][sample][0] += getInt(meth)
        count[idx[k]][sample][1] += getInt(unmeth)
      else:
        count[idx[k]][sample] = [getInt(meth), getInt(unmeth)]

def printClus(infiles, files, samples, minReg, \
    fraction, fOut, count, clus):
  '''
  Print counts for each cluster.
  '''
  printed = 0
  for i in range(len(clus)):
    res = clus[i]
    for sample in samples:

      if sample in count[i]:
        meth = count[i][sample][0]
        unmeth = count[i][sample][1]
        if meth + unmeth < minReg:
          # less than minimum number of counts
          res += '\tNA'
          if not fraction:
            res += '\tNA'
        else:
          if fraction:
            # compute methylated fraction
            res += '\t%f' % (meth / float(meth + unmeth))
          else:
            res += '\t%d\t%d' % (meth + unmeth, meth)  # actual counts
      else:
        res += '\tNA'
        if not fraction:
          res += '\tNA'

    fOut.write(res + '\n')
    printed += 1

  return printed

def processFilesClus(infiles, files, samples, minReg, \
    fraction, fOut, clusFile, verbose):
  '''
  Control processing of files when given file of
    cluster information.
  '''
  # load cluster information from file
  clus, idx = loadClusters(clusFile, verbose)

  # load methylation information for each sample
  count = [{} for i in range(len(clus))]
  if verbose:
    sys.stderr.write('Loading methylation information\n')
  for i in range(len(infiles)):
    if verbose:
      sys.stderr.write('  file: %s\n' % infiles[i])
    loadCountsClus(files[i], count, idx, samples[i])

  # print output for each cluster
  if verbose:
    sys.stderr.write('Printing output\n')
  return printClus(infiles, files, samples, minReg, \
    fraction, fOut, count, clus)

def loadSites(f, line, refChrom, idx, clus):
  '''
  Load cluster information for one chrom.
  '''
  total = cpg = 0
  while line:
    spl = line.split('\t')
    if len(spl) < 2:
      sys.stderr.write('Error! Poorly formatted cluster record:' \
        + '\n%s' % line)
      sys.exit(-1)
    if spl[0] != refChrom:
      break
    pos = spl[1].split(',')

    # save sites and cluster header
    for p in pos:
      idx[spl[0] + ' ' + p] = total
      cpg += 1
    clus.append('\t'.join([spl[0], pos[0], pos[-1], \
      str(len(pos))]))
    total += 1

    line = f.readline().rstrip()

  return line, cpg

def loadCountsClusChrom(files, samples, lines, refChrom, \
      count, idx):
  '''
  Load the methylated/unmethylated counts for a set
    of files for one chrom (pre-clustered option).
  '''
  # load counts from each file
  for i in range(len(files)):
    if not lines[i]: continue
    try:
      chrom, pos, end, pct, meth, unmeth \
        = lines[i].rstrip().split('\t')
    except ValueError:
      sys.stderr.write('Error! Poorly formatted record:' \
        + '\n%s' % lines[i])
      sys.exit(-1)
    if chrom != refChrom: continue

    # save counts as long as chrom matches
    while chrom == refChrom:
      k = chrom + ' ' + pos
      if k in idx:
        if samples[i] in count[idx[k]]:
          count[idx[k]][samples[i]][0] += getInt(meth)
          count[idx[k]][samples[i]][1] += getInt(unmeth)
        else:
          count[idx[k]][samples[i]] = [getInt(meth), getInt(unmeth)]

      # load next line
      line = files[i].readline()
      if not line:
        break
      try:
        chrom, pos, end, pct, meth, unmeth \
          = line.rstrip().split('\t')
      except ValueError:
        sys.stderr.write('Error! Poorly formatted record:' \
          + '\n%s' % line)
        sys.exit(-1)

    # reset next line
    lines[i] = line

def processChromClus(infiles, files, samples, minReg, \
    fraction, fOut, chrOrder, clusFile, verbose):
  '''
  Control processing of files when given file of
    cluster information, one chrom at a time.
  '''
  if verbose:
    sys.stderr.write('Loading cluster information and ' \
      'methylation counts,\n  and producing output\n')

  # load first lines of methylation files
  lines = []
  for f in files:
    lines.append(f.readline())

  # load first line of cluster info
  f = openRead(clusFile)
  line = f.readline().rstrip()

  total = sites = 0
  printed = 0

  # process each chrom
  for refChrom in chrOrder:
    if verbose:
      sys.stderr.write('  chromosome: %s\n' % refChrom)

    # load cluster sites
    idx = dict()
    clus = list()
    line, cpg = loadSites(f, line, refChrom, idx, clus)
    total += len(clus)
    sites += cpg

    # load counts
    count = [{} for i in range(len(clus))]
    loadCountsClusChrom(files, samples, lines, refChrom, \
      count, idx)

    # print output
    printed += printClus(infiles, files, samples, minReg, \
      fraction, fOut, count, clus)

  # check for unprocessed records
  unProc = False
  for i in range(len(lines)):
    if lines[i]:
      if not unProc:
        sys.stderr.write('Error! Unprocessed records in inputs:\n')
      sys.stderr.write('File %s:\n%s' % (infiles[i], lines[i]))
      unProc = True
  if unProc:
    sys.exit(-1)

  if f != sys.stdin:
    f.close()
  if verbose:
    sys.stderr.write('  Clusters loaded: %d\n' % total)
    sys.stderr.write('  CpG sites: %d\n' % sites)
  return printed

def main():
  '''
  Main.
  '''
  # Default parameters
  infiles = []        # list of input files
  outfile = None      # output file
  minReads = 3        # min. reads in a sample at a position
  minSamples = 1      # min. samples with min. reads at a position
  maxDist = 100       # max. distance between CpGs
  minCpG = 3          # min. CpGs in a region
  minReg = 20         # min. reads in a sample for a region
  maxLen = 500        # max. length of a combined region
  fraction = False    # report methylated fractions option
  clusFile = None     # file listing clusters of valid CpGs
  byChrom = False     # process by chromosome (memory saving)
  chrOrderFile = None # file listing order of chromosomes to process
  verbose = False     # verbose option

  # Get command-line args
  args = sys.argv[1:]
  i = 0
  while i < len(args):
    if args[i] == '-h' or args[i] == '--help':
      usage()
    elif args[i] == '--version':
      printVersion()
    elif args[i] == '-v':
      verbose = True
    elif args[i] == '-b':
      byChrom = True
    elif args[i] == '-f':
      fraction = True
    elif args[i][0] == '-' and i < len(args) - 1:
      if args[i] == '-r':
        minReads = getInt(args[i+1])
      elif args[i] == '-s':
        minSamples = getInt(args[i+1])
      elif args[i] == '-d':
        maxDist = getInt(args[i+1])
      elif args[i] == '-c':
        minCpG = getInt(args[i+1])
      elif args[i] == '-m':
        minReg = getInt(args[i+1])
      elif args[i] == '-x':
        maxLen = getInt(args[i+1])
      elif args[i] == '-e':
        chrOrderFile = args[i+1]
      elif args[i] == '-o':
        outfile = args[i+1]
      elif args[i] == '-y':
        clusFile = args[i+1]
      else:
        sys.stderr.write('Error! Unknown parameter: %s\n' % args[i])
        usage()
      i += 1
    elif os.path.isfile(args[i]):
      infiles.append(args[i])
    else:
      sys.stderr.write('Error! Unknown parameter with no arg: ' \
        + '%s\n' % args[i])
      usage()
    i += 1

  # check for I/O errors
  if outfile == None:
    sys.stderr.write('Error! Must specify an output file\n')
    usage()
  if len(infiles) == 0:
    sys.stderr.write('Error! Must specify one or more input files\n')
    usage()
  fOut = openWrite(outfile)
  files, samples = openFiles(infiles)
  writeHeader(fOut, samples, fraction)

  if byChrom:
    chrOrder = saveChrOrder(chrOrderFile, infiles, files, \
      verbose)

  # if provided clusters, process directly
  if clusFile and os.path.isfile(clusFile):
    if byChrom:
      printed = processChromClus(infiles, files, samples, \
        minReg, fraction, fOut, chrOrder, clusFile, verbose)
    else:
      printed = processFilesClus(infiles, files, samples, \
        minReg, fraction, fOut, clusFile, verbose)

  # default analysis: cluster samples and produce output
  else:
    # open output file for cluster info
    fClus = None
    if clusFile:
      fClus = openWrite(clusFile)

    # cluster and produce output
    if byChrom:
      printed = processChrom(infiles, files, samples, minReads, \
        minSamples, maxDist, minCpG, minReg, maxLen, fraction, \
        fOut, chrOrder, fClus, verbose)
    else:
      printed = processFiles(infiles, files, samples, minReads, \
        minSamples, maxDist, minCpG, minReg, maxLen, fraction, \
        fOut, fClus, verbose)

    if fClus and fClus != sys.stdout:
      fClus.close()

  # finish up
  if fOut != sys.stdout:
    fOut.close()
  if verbose:
    sys.stderr.write('Genomic regions printed: %d\n' % printed)
    sys.stderr.write('Valid sample names:')
    for sample in samples:
      sys.stderr.write(' ' + sample)
    sys.stderr.write('\n')

if __name__ == '__main__':
  main()
