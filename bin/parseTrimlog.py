#!/usr/bin/env python2
__doc__='''
##########################################################################
## parseTrimlog.py | Simple read count statistics from trimmomatic log  ##
##########################################################################
'''

__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT"
__version__ = "4.0.7"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
import re
from collections import Counter, defaultdict, OrderedDict
import datetime
import numpy as np

# SE
# 3JTQU:01170:00546 234 12 246 7

# PE (illumina)
# M03623:13:000000000-AEUCW:1:2115:9806:16147 2:N:0:2 288 0 288 13

# read name (and group)
# surviving sequence length
# location of the first surviving base, aka. the amount trimmed from the start
# location of the last surviving base in the original read
# amount trimmed from the end

# miseq/nextseq PE
# SRA
# SE
parsers = [
    re.compile('([^:]+:[^:]+:[^:]+:[^:]+:(\d+):\d+:\d+)\s+(?:\S+\s+)?(\d+)\s+\d+\s+\d+\s+(\d+)$')
]

# parse trimming stats
re_trimlog = None
success = Counter()
lengths = defaultdict(list)
tree = lambda: defaultdict(tree)
tiletrim = tree()
passed = Counter()  # counts nonzero pairs
mates = []
mates_rlen = []

for i, line in enumerate(sys.stdin):
    if not re_trimlog:
        ## find best parser
        for parser in parsers:
            print >> sys.stderr, "INFO: testing pattern", parser.pattern
            m = parser.match(line)
            if m:
                re_trimlog = parser
                print >> sys.stderr, 'INFO: parsed', m.groups()
                break
        if not re_trimlog:
            print >> sys.stderr, "WARNING: No suitable parser, no trimming stats will be calculated"
            print >> sys.stdout, "## custom.read.metrics"
            print >> sys.stdout, "# parseTrimlog.py STDIN"
            print >> sys.stdout, "## custom.read.metrics"
            print >> sys.stdout, "# Started on:", datetime.datetime.now().isoformat()
            print >> sys.stdout, "\n# Sorry, no suitable parser available"
            print >> sys.stdout  # add empty line
            sys.exit(0)

    try:
        m = re_trimlog.match(line)
        assert m
    except AssertionError:
        print >> sys.stderr, "WARNING: parser error, skipping line", i
        continue

    # parse line
    fragment, tile, rlen, tlen = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
    # set group (count consecutive read names) and mates read length (for fragment pass stats)
    if mates and mates[-1] == fragment:
        mates.append(fragment)
        mates_rlen.append(rlen)
    else:
        mates = [ fragment ]
        # update tile stats for mates
        passed[tile] += 1 if all(mates_rlen) else 0
        mates_rlen = [ rlen ]
    group = len(mates)

    # count success (read is kept)
    if rlen != 0:
        success[group] += 1

    # store read length
    lengths[group].append(rlen)
    try:
        tiletrim[tile][group][rlen] += 1
    except TypeError:
        tiletrim[tile][group][rlen] = 1
    except:
        raise

## count last mate group
passed[tile] += 1 if all(mates_rlen) else 0


# calculcate tile statistics
tilestatlines = []
rlen_tiles = defaultdict(list)
for tile in sorted(tiletrim.keys()):
    median_read_len = []
    for grp in sorted(tiletrim[tile].keys()):
        rsum = sum(tiletrim[tile][grp].values())
        previous, s = 0,0  # s is count of all reads
        for rlen in sorted(tiletrim[tile][grp].keys()):
            s += tiletrim[tile][grp][rlen]
            if s > rsum/2:  # passed median
                tile_median = (previous+rlen)/2 if previous and rlen-previous > 1 else rlen
                break
            previous = rlen
        median_read_len.append(tile_median)
        rlen_tiles[grp].append(tile_median)
    tilestatlines.append('\t'.join(map(str,[tile, ','.join(map(str,median_read_len)), passed[tile]])))

# summary statistics
stats = OrderedDict()
stats['TOTAL_READS'] = ','.join([str(len(lengths[group])) for group in sorted(lengths.keys())])
stats['PASSED_READS'] = ','.join([str(success[group]) for group in sorted(lengths.keys())])
stats['PASSED_FRAC'] = ','.join(['{:.3f}'.format(success[group]/float(len(lengths[group]))) for group in sorted(lengths.keys())])
stats['LENGTH_MEDIAN'] = ','.join([str(np.median(lengths[group])) for group in  sorted(lengths.keys())])
stats['LENGTH_MEAN'] = ','.join([str(np.mean(lengths[group])) for group in  sorted(lengths.keys())])
stats['LENGTH_STD'] = ','.join([str(np.std(lengths[group])) for group in  sorted(lengths.keys())])
stats['TILE_MEDIAN_LENGTH_STD'] = ','.join([ '{:.3f}'.format(np.std(rlen_tiles[r])) for r in sorted(rlen_tiles.keys()) ])
stats['TILE_MEDIAN_LENGTH_MIN'] = ','.join([ str(min(rlen_tiles[r])) for r in sorted(rlen_tiles.keys())])
stats['TILE_MEDIAN_LENGTH_MAX'] = ','.join([ str(max(rlen_tiles[r])) for r in sorted(rlen_tiles.keys())])


# print to metrics file
print >> sys.stdout, "## custom.read.metrics"
print >> sys.stdout, "# parseTrimlog.py STDIN"
print >> sys.stdout, "## custom.read.metrics"
print >> sys.stdout, "# Started on:", datetime.datetime.now().isoformat()
print >> sys.stdout, "\n# BASIC FILTERING METRICS"
print >> sys.stdout, '\t'.join(stats.keys())
print >> sys.stdout, '\t'.join([ stats[k] for k in stats.keys() ])
# summarize trimming lengths per tile
print >> sys.stdout, "\n# MEDIAN READ LENGTH PER TILE"
print >> sys.stdout, '\t'.join(['TILE','MEDIAN_LENGTH','NONZERO'])
print >> sys.stdout, '\n'.join(tilestatlines)


print >> sys.stdout  # add empty line
