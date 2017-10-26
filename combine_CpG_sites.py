#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Dec. 2016

# Combining multiple samples' single-base resolution
#   methylation data into a set of genomic regions.

import sys
import os.path
import gzip
import math
version = '0.3'
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
      -b          Memory-saving option (may take longer)
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

def splitRegion(chrom, reg, count, minCpG, minReg, \
    maxLen, samples, fraction, fOut):
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
      minCpG, minReg, float('inf'), samples, fraction, fOut)
    start = end

  return total

def processRegion(chrom, reg, count, minCpG, minReg, \
    maxLen, samples, fraction, fOut):
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
      maxLen, samples, fraction, fOut)

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
      flag = True
  if flag:
    fOut.write(res + '\n')
    return 1
  return 0

def combineRegions(count, total, chrom, minSamples, maxDist, \
    minCpG, minReg, maxLen, samples, fraction, fOut):
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
          minReg, maxLen, samples, fraction, fOut)
        reg = []  # reset list
      reg.append(loc)
      pos3 = loc
  # process last genomic region for this chromosome
  printed += processRegion(chrom, reg, count, minCpG, \
    minReg, maxLen, samples, fraction, fOut)
  return printed

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

def processFiles(infiles, files, samples, minReads,
    minSamples, maxDist, minCpG, minReg, maxLen, fraction,
    fOut, verbose):
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
      samples, fraction, fOut)
  return printed

def loadChrOrder(infiles, files):
  '''
  Return a list of ordered chromosome names.
  '''
  # load chroms from input files
  order = [[] for i in range(len(files))]
  for i in range(len(files)):
    refChrom = ''
    for line in files[i]:
      chrom = line.split('\t')[0]
      if chrom != refChrom:
        if chrom in order[i]:
          sys.stderr.write('Error! Unsorted input file:' + \
            ' %s\n' % infiles[i])
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

def loadChromCounts(files, samples, lines, refChrom,
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
      sys.stderr.write('Error! Poorly formatted record:' + \
        '\n%s' % lines[i])
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
        sys.stderr.write('Error! Poorly formatted record:' + \
          '\n%s' % line)
        sys.exit(-1)

    # reset next line
    lines[i] = line

def processChrom(infiles, files, samples, minReads,
    minSamples, maxDist, minCpG, minReg, maxLen,
    fraction, fOut, chrOrder, verbose):
  '''
  Process the input files, one chromosome at a time.
  '''
  # load ordered chromosome names (if not specified)
  if not chrOrder:
    if verbose:
      sys.stderr.write('Loading chromosome order\n')
    chrOrder = loadChrOrder(infiles, files)

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
    loadChromCounts(files, samples, lines, refChrom,
      minReads, count, total)

    # cluster and print output
    printed += combineRegions(count, total, refChrom, \
      minSamples, maxDist, minCpG, minReg, maxLen, \
      samples, fraction, fOut)

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

  # run program
  if byChrom:
    chrOrder = []   # order of chromosomes to process
    if chrOrderFile != None:
      # load chromosome order from given file
      f = openRead(chrOrderFile)
      for line in f:
        chrOrder.extend(line.rstrip().split(','))
      if f != sys.stdin:
        f.close()
    printed = processChrom(infiles, files, samples, minReads, \
      minSamples, maxDist, minCpG, minReg, maxLen, fraction, \
      fOut, chrOrder, verbose)
  else:
    printed = processFiles(infiles, files, samples, minReads, \
      minSamples, maxDist, minCpG, minReg, maxLen, fraction, \
      fOut, verbose)

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
