#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Dec. 2016

# This script produces a sorted table of methylated/unmethylated
#   counts directly from a SAM file made by Bismark.

import sys
import os.path
import re
import gzip
version = '0.1'
copyright = 'Copyright (C) 2016 John M. Gaspar (jsh58@wildcats.unh.edu)'

def printVersion():
  sys.stderr.write('extract_CpG_data.py from DMRfinder, version' \
    + ' %s\n' % version)
  sys.stderr.write(copyright + '\n')
  sys.exit(-1)

def usage():
  sys.stderr.write('''Usage: python extract_CpG_data.py  [options]  -i <input>  -o <output>
    -i <input>    SAM alignment file produced by Bismark (must have
                    a header and 'XM' methylation strings; can use
                    '-' for stdin)
    -o <output>   Output file listing counts of methylated and
                    unmethylated CpGs, merged and sorted
  Options:
    -m <int>      Minimum coverage (methylation counts) to report a
                    CpG site (def. 1 [all sites reported])
    -s            Report strand in third column of output
    -b <file>     BED file listing regions for which to collect
                    linked methylation data
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
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
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

def loadBed(bedFile, bedRegions, bedSites):
  '''
  Load BED regions from file. Open output file.
  '''
  # load BED regions, do not allow repeated names
  f = openRead(bedFile)
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 4:
      sys.stderr.write('Error! Poorly formatted BED record:\n%s' % line)
      sys.exit(-1)
    if spl[3] in bedRegions:
      sys.stderr.write('Error! Repeated BED region name: %s\n' % spl[3])
      sys.exit(-1)
    # save chrom, start, and end
    bedRegions[spl[3]] = (spl[0], getInt(spl[1]), getInt(spl[2]))
    bedSites[spl[3]] = [] # CpG sites will be loaded in parseSAM()
  f.close()

  # open output file (strip suffix and add '_linked.txt')
  fname = '.'.join(bedFile.split('.')[:-1])
  while os.path.isfile(fname + '_linked.txt'):
    fname += '-'
  return openWrite(fname + '_linked.txt')

def printBed(bedOut, bedRegions, bedSites, linkedMeth):
  '''
  Print linked-methylation information for designated regions.
  '''
  for reg in sorted(bedRegions):
    # print header
    chrom = bedRegions[reg][0]
    bedOut.write('Region: %s, %s:%d-%d\nSites: ' % (reg, \
      chrom, bedRegions[reg][1], bedRegions[reg][2]))
    count = 0
    # print list of sorted CpG sites
    for pos in sorted(bedSites[reg]):
      if count:
        bedOut.write(', %d' % pos)
      else:
        bedOut.write('%d' % pos)
      count += 1
    if count == 0:
      bedOut.write('<none>\n\n')  # no CpG sites
      continue
    bedOut.write('\n')

    # compile methylation results for each read
    for head in linkedMeth[reg]:
      res = ''  # result string -- meth data for each read
      for pos in sorted(bedSites[reg]):
        if pos not in linkedMeth[reg][head][chrom]:
          res += '-'  # no data: labeled '-'
          continue
        unmeth, meth = linkedMeth[reg][head][chrom][pos]
        if unmeth == 1:
          res += '0'  # unmethylated: labeled '0'
        elif meth == 1:
          res += '1'  # methylated: labeled '1'
        else:
          sys.stderr.write('Error! Problem parsing linked '
            + 'methylation information for read %s\n' % head)
          sys.exit(-1)
      bedOut.write('%s\t%s\n' % (res, head))
    bedOut.write('\n')

def printOutput(fOut, genome, meth, minCov, strand, verbose):
  '''
  Print the sorted output -- location and methylation counts
    for each CpG with at least minCov data values.
  '''
  if verbose:
    sys.stderr.write('Printing the output file\n')
  printed = 0
  for chrom in genome:
    for loc in sorted(meth[chrom]):
      total = meth[chrom][loc][0] + meth[chrom][loc][1]
      if total < minCov:
        continue  # fails to meet minimum coverage
      if strand:
        # 3rd column is strand ('+')
        fOut.write('%s\t%d\t+\t%.6f\t%d\t%d\n' % (chrom, loc,
          100.0 * meth[chrom][loc][1] / total,
          meth[chrom][loc][1], meth[chrom][loc][0]))
      else:
        # 3rd column is end (loc+1)
        fOut.write('%s\t%d\t%d\t%.6f\t%d\t%d\n' % (chrom, loc, loc+1,
          100.0 * meth[chrom][loc][1] / total,
          meth[chrom][loc][1], meth[chrom][loc][0]))
      printed += 1
  return printed

def parseCigar(cigar):
  '''
  Return string representation of CIGAR.
  '''
  ops = re.findall(r'(\d+)([IDM])', cigar)
  cigar = ''
  for op in ops:
    cigar += int(op[0]) * op[1]
  return cigar

def getTag(lis, tag):
  '''
  Get optional tag from a SAM record.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  sys.stderr.write('Error! Cannot find %s in SAM record\n' % tag)
  sys.exit(-1)

def saveMeth(d, chrom, loc, meth):
  '''
  Save the methylation info for a given genomic position
    to the given dict (d).
  '''
  if chrom not in d:
    d[chrom] = {}
  if loc not in d[chrom]:
    d[chrom][loc] = [0, 0]
  if meth == 'z':
    d[chrom][loc][0] += 1  # unmethylated: index 0
  else:
    d[chrom][loc][1] += 1  # methylated: index 1

def loadMeth(cigar, strXM, chrom, pos, rc, meth, ins, dup,
    peMeth, head, bedRegions, bedSites, linkedMeth):
  '''
  Load methylation info using a methylation string (strXM).
  '''
  methCount = count = 0  # counting variables
  offset = 0  # in/del offset
  cigPos = 0  # position in cigar -- mirrors i (position in meth)
  for i in range(len(strXM)):

    # CpG methylation calls are 'z' (unmethylated) or 'Z' (methylated)
    while strXM[i] in ['z', 'Z']:

      # determine genomic location (C of 'CG' on the forward strand)
      loc = pos + i - rc + offset

      # for "novel" CpG sites created by a deletion,
      #   adjust location to the 5' end of the 'D's
      if rc and cigar[cigPos-1] == 'D':
        j = cigPos - 1
        while j > -1 and cigar[j] != 'M':
          j -= 1
        loc -= cigPos - j - 1

      # skip if position has been counted in a previous alignment
      if dup == 2 and chrom in peMeth and loc in peMeth[chrom]:
        break

      # for "novel" CpG sites created by an insertion,
      #   save to 'ins' dictionary
      if (not rc and cigar[cigPos] == 'I') or \
          (rc and cigPos > 0 and cigar[cigPos-1] == 'I'):
        saveMeth(ins, chrom, loc, strXM[i])
      # otherwise, save to regular 'meth' dict
      else:
        saveMeth(meth, chrom, loc, strXM[i])

      # if p-e alignment, also save to 'peMeth' dict
      if dup == 1:
        saveMeth(peMeth, chrom, loc, strXM[i])

      # save methylation data if it falls within a BED region
      for reg in bedRegions:
        if bedRegions[reg][0] == chrom and bedRegions[reg][1] <= loc \
            and bedRegions[reg][2] > loc:
          # save meth data to linkedMeth dict (using read header)
          if reg not in linkedMeth:
            linkedMeth[reg] = {}
          if head not in linkedMeth[reg]:
            linkedMeth[reg][head] = {}
          saveMeth(linkedMeth[reg][head], chrom, loc, strXM[i])
          # save location to bedSites dict
          if loc not in bedSites[reg]:
            bedSites[reg].append(loc)

      # update counts
      if strXM[i] == 'Z':
        methCount += 1
      count += 1
      break

    # change in/del offset
    cigPos += 1
    while cigPos < len(cigar) and cigar[cigPos] == 'D':
      offset += 1
      cigPos += 1
    if cigPos < len(cigar) and cigar[cigPos] == 'I':
      offset -= 1

  return methCount, count

def parseSAM(f, meth, bedRegions, bedSites, linkedMeth, verbose):
  '''
  Parse the SAM file. Save methylation data.
  '''
  if verbose:
    sys.stderr.write('Processing the SAM file\n')
  genome = []  # for chromosome names, ordered in SAM header
  ins = {}     # for novel CpGs caused by insertion
  peMeth = {}  # for checking overlapping of paired-end alignments
  total = mapped = 0     # counting variables for reads
  methCount = count = 0  # counting variables for methylation data
  for line in f:

    # save chromosome names
    if line[0] == '@':
      spl = line.rstrip().split('\t')[1].split(':')
      if line[1:3] == 'SQ' and spl[0] == 'SN':
        meth[spl[1]] = {}
        genome.append(spl[1])
      continue

    # load SAM data
    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record\n' + line)
      sys.exit(-1)

    # get alignment info
    flag = getInt(spl[1])
    if flag & 0x900:
      continue  # skip secondary/supplementary
    total += 1
    if verbose and total % 1000000 == 0:
      sys.stderr.write('Reads processed so far: %d\n' % total)
    if flag & 0x4:
      continue  # skip unmapped
    mapped += 1
    rc = 0
    if getTag(spl[11:], 'XG') == 'GA':
      rc = 1  # alignment is to G->A converted genome,
              #   so methylation data is on G of 'CG'

    # determine if read has multiple segments
    dup = 0  # 0 -> single-end alignment
             # 1 -> paired-end alignment, not seen before
             # 2 -> paired-end alignment, seen before
    if flag & 0x1:
      if spl[0] in peMeth:
        dup = 2
      else:
        # first segment -- initialize dict
        peMeth[spl[0]] = {}
        dup = 1

    # load chrom, position
    chrom = spl[2]
    if chrom not in meth:
      sys.stderr.write('Error! Cannot find chromosome %s' % chrom \
        + ' in genome\n  (make sure input SAM has a header)\n')
      sys.exit(-1)
    pos = getInt(spl[3])

    # load CIGAR, methylation string
    cigar = parseCigar(spl[5])
    strXM = getTag(spl[11:], 'XM')  # methylation string from Bismark

    # load CpG methylation info
    count1, count2 = loadMeth(cigar, strXM, chrom, pos, rc, \
      meth, ins, dup, peMeth.get(spl[0], None), \
      spl[0], bedRegions, bedSites, linkedMeth)
    methCount += count1
    count += count2

  # warn about novel inserted CpGs
  if verbose and ins:
    sys.stderr.write('Warning! Novel CpG(s) caused by ' \
      + 'insertion(s) -- will be ignored.\n')

  return genome, total, mapped, methCount, count

def main():
  '''
  Main.
  '''
  # Default parameters
  infile = None     # input file
  outfile = None    # output file
  minCov = 1        # min. coverage to report a CpG site
  strand = False    # report strand in output
  bedFile = None    # (optional) BED file for linked meth. data
  verbose = False   # verbose option

  # get command-line args
  args = sys.argv[1:]
  i = 0
  while i < len(args):
    if args[i] == '-h' or args[i] == '--help':
      usage()
    elif args[i] == '--version':
      printVersion()
    elif args[i] == '-v':
      verbose = True
    elif args[i] == '-s':
      strand = True
    elif i < len(args) - 1:
      if args[i] == '-i':
        infile = openRead(args[i+1])
      elif args[i] == '-o':
        outfile = openWrite(args[i+1])
      elif args[i] == '-m':
        minCov = max(getInt(args[i+1]), 1)
      elif args[i] == '-b':
        bedFile = args[i+1]
      else:
        sys.stderr.write('Error! Unknown parameter: %s\n' % args[i])
        usage()
      i += 1
    else:
      sys.stderr.write('Error! Unknown parameter with no arg: ' \
        + '%s\n' % args[i])
      usage()
    i += 1

  # check for I/O errors
  if infile == None or outfile == None:
    sys.stderr.write('Error! Must specify input and output files\n')
    usage()

  # load BED file regions (optional)
  bedRegions = {}  # defining genomic regions for each BED record
  bedSites = {}    # CpG sites for each region (loaded in parseSAM())
  if bedFile != None:
    bedOut = loadBed(bedFile, bedRegions, bedSites)
    if verbose:
      sys.stderr.write('BED regions loaded: %d\n' % len(bedRegions))

  # process file
  meth = {}        # dict for methylation counts
  linkedMeth = {}  # dict for linked methylation data
  genome, total, mapped, methCount, count = parseSAM(infile, meth,
    bedRegions, bedSites, linkedMeth, verbose)
  if infile != sys.stdin:
    infile.close()

  # print summary counts
  if verbose:
    sys.stderr.write('Reads analyzed: %d\n' % total \
      + '  Mapped: %d\n' % mapped \
      + '  Total CpG methylation values in the reads: %d\n' % count \
      + '    Methylated: %d\n' % methCount \
      + '    Unmethylated: %d\n' % (count - methCount))
    if count:
      sys.stderr.write('    Percent methylated: %.1f%%\n' % \
        (100.0 * methCount / count))
    else:
      sys.stderr.write('    Percent methylated: n/a\n')

  # print output file
  printed = printOutput(outfile, genome, meth, minCov, strand, \
    verbose)
  if outfile != sys.stdout:
    outfile.close()
  if bedFile != None:
    printBed(bedOut, bedRegions, bedSites, linkedMeth)
    bedOut.close()
  if verbose:
    sys.stderr.write('Genomic CpG sites printed: %d\n' % printed)

if __name__ == '__main__':
  main()
