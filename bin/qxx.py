#!/usr/bin/env python2
__doc__='''
################################################################
## q30.py | Calculate Q30 metrics from picard quality by base ##
################################################################
'''

__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT"
__version__ = "4.0.7"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys, os
import time, datetime
from collections import OrderedDict, defaultdict, Counter
ts = time.time()

'''list statisitics'''
def QXX(data, thresholds):
    results = Counter()
    for t in thresholds:
        count, total = 0, 0
        for q, c in data.items():
            try:
                qual = int(q)
                reads = int(c)
            except:
                pass
            else:
                if qual>=t:
                    count += reads
                total += reads
        if total:
            results[t] = '{:.3f}'.format(float(count)/total)
    return results

if __name__=="__main__":
    breaks = [10,15,20,25,30,35,40]
    stats = OrderedDict()
    for fi in sys.argv[1:]:
        data = {}
        started = False
        # parse input
        with open(fi) as fh:
            for line in fh:
                if started:
                    try:
                        f = line.split()
                        data[int(f[0])] = int(f[1])
                    except:
                        pass
                elif line.startswith('QUALITY'):
                    started = True
        # Q30 statistics
        stats[fi] = QXX(data,breaks)

    # print stats to stdout
    print >> sys.stdout, "## custom.qualitydistribution.metrics"
    print >> sys.stdout, "#", ' '.join(sys.argv)
    print >> sys.stdout, "## custom.qualitydistribution.metrics"
    print >> sys.stdout, "# Started on:", datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print >> sys.stdout, "\n# QUALITY DISTRIBUTION METRICS"
    print >> sys.stdout, '\t'.join(['FILENAME'] + [ 'Q'+x for x in map(str,breaks) ] )
    for fi, result in stats.items():
        print >> sys.stdout, '\t'.join([os.path.basename(fi)] + map(str,[ result[b] for b in breaks ]))
    print >> sys.stdout  # add empty line
