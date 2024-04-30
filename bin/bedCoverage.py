#!/usr/bin/env python2

__doc__=="""
calculates coverage over exons from an overlap as such:
bedtools intersect -wao -a EXONFILE.BED6 -b BEDGRAPH (eg. from bedtools genomecov)
Returns 'picard-style' metrics
  1. EXON_LENGTH, COVERED_BASES, COVERED_PCT, COVERED_XXFOLD_PCT (as specified by steps)
  2. not covered intervals (zero coverage)
  3. coverage by gene GENE covered, covered_thresholded
"""
__author__ = "David Brawand"
__credits__ = "David Brawand"
__license__ = "MIT"
__version__ = "4.0.7"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
import datetime, time
from collections import Counter, OrderedDict

ts = time.time()

## chr     start           end             name            score   strand  chr     start           end             score   length
## 7       116422041       116422151       NM_000245       0       +       7       116422034       116422044       458     3
## 7       116422041       116422151       NM_000245       0       +       7       116422044       116422045       457     1
## 7       116422041       116422151       NM_000245       0       +       7       116422045       116422052       458     7
## 7       116422041       116422151       NM_000245       0       +       7       116422052       116422053       457     1
## 7       55268008        55268106        NM_005228       0       +       7       55268104        55268106        2686    2
## 7       39893941        39899991        FAKE            0       +       .       -1              -1              .       0


def main(lowcoverthreshold,bins):
    # read intersection data
    exons = {}
    for line in sys.stdin:
        f = line.split()
        # check if correctly formatted
        formatcheck(f)
        # add
        exonID = tuple(f[:4])
        try:
            exons[exonID].extend(f)
        except KeyError:
            exons[exonID] = Exon(f)
        # # skip 0 coverage
        # try:
        #     cov = int(f[7])
        #     assert cov != 0
        # except:
        #     continue

    # get coverage lengths
    notCovered = []
    coverage = Counter()
    total = 0
    for e in exons.values():
        total += len(e)
        c = e.getCoverage()
        for c,l in c.items():
            coverage[c] += l

    # generate coveragebins (covered equal or more)
    coveragebins = Counter()
    coveredbases = 0
    for k,v in coverage.items():
        if k>0:
            coveredbases += v
        for c in bins:
            if k>=c:
                coveragebins[c] += v

    # generate stats dict
    stats = OrderedDict()
    stats['EXON_LENGTH'] = sum([len(e) for e in exons.values()])
    stats['COVERED_BASES'] = coveredbases
    for c in sorted(bins):
        stats['COVERED_'+str(c)+'X_PCT'] = '{:.3f}'.format(coveragebins[c]/float(stats['EXON_LENGTH']))

    # print stats
    print >> sys.stdout, "## custom.coverage.metrics"
    print >> sys.stdout, "# Started on:", datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print >> sys.stdout, "## custom.coverage.metrics"
    print >> sys.stdout, "\n# COVERAGE BINS"
    print >> sys.stdout, '\t'.join(stats.keys())
    print >> sys.stdout, '\t'.join(map(str,stats.values()))
    ## no coverage (0X)
    print >> sys.stdout, "\n# NOT COVERED"
    print >> sys.stdout, "chr\tstart\tend\tname\tcov"
    for e in sorted(exons.values()):
        for c in e.notCovered(0):
            print >> sys.stdout, '\t'.join([e.chr,str(c[0]),str(c[1]),e.name,str(c[2])])
    ## low coverage (below 20X)
    print >> sys.stdout, "\n# LOW COVERAGE (<"+str(lowcoverthreshold)+"X)"
    print >> sys.stdout, "chr\tstart\tend\tname\tavgcov"
    for e in sorted(exons.values()):
        for c in e.notCovered(lowcoverthreshold-1):  # <20x by default
            print >> sys.stdout, '\t'.join([e.chr,str(c[0]),str(c[1]),e.name,str(c[2])])
    ## percentage covered (lowcoverage threshold and bins)
    coveragebins = sorted(list(set([lowcoverthreshold] + bins)))
    cov = {}
    print >> sys.stdout, "\n# PERCENT COVERED"
    print >> sys.stdout, "name\t" + '\t'.join(["PCT"+str(x)+"X" for x in coveragebins])
    for e in sorted(exons.values()):
        try:
            cov[e.name]['length'] += len(e)
        except KeyError:
            cov[e.name] = { 'length': len(e), 'covered': OrderedDict(zip(coveragebins,[0]*len(coveragebins))) }
        except:
            raise
        # add coverage
        for c in cov[e.name]['covered'].keys():
            cov[e.name]['covered'][c] += e.coveredBases(c)
    for gene in sorted(cov.keys()):
        # caluculate coverage
        for c in cov[gene]['covered'].keys():
           cov[gene]['covered'][c] /= float(cov[gene]['length'])
        # print result
        print >> sys.stdout, '\t'.join([gene]+['{:.3f}'.format(cov[gene]['covered'][c]) for c in cov[gene]['covered'].keys()])
    # empty line to terminate output
    print >> sys.stdout


def formatcheck(f):
    try:
        assert len(f) == 9
        int(f[1])  #FROM
        int(f[2])  #TO
        int(f[5])  #FROM
        int(f[6])  #TO
        assert f[7]=='.' or int(f[7])>=0  #SCORE
        int(f[8])  #LENGTH
    except (ValueError, AssertionError):
        sys.stderr.write("please input with 'bedtools intersect -wao -a EXONFILE.BED4 -b BEDGRAPH'\n")
        sys.stderr.write(">"+str(f)+"<\n")
        raise Exception('FormatError')


class Exon(object):
    def __init__(self,f):
        self.exonID = tuple(f[:4])
        self.chr = f[0]
        self.coord = tuple(map(int,f[1:3]))
        try:
            assert self.coord[0] <= self.coord[1]
        except:
            raise
        self.name = f[3]
        self.coverage = []
        self.extend(f)  # add coverage segment
        return

    def __len__(self):
        return self.coord[1] - self.coord[0]

    def __lt__(self,other):
        return (self.chr, self.coord[0], self.coord[1]) < (other.chr, other.coord[0], other.coord[1])

    def extend(self,f):
        if f[4] != '.':
            cc = sorted(list(self.coord) + map(int,f[5:7]))
            assert len(cc)==4
            covered = cc[1:3]
            assert covered[1]-covered[0]==int(f[8])
            self.coverage.append(covered + [int(f[7])])
        return

    '''returns length of limit<=coverage'''
    def coveredBases(self, limit=1):
        return sum([x[1]-x[0] for x in self.coverage if x[2]>=limit])

    '''returns dict of coverage including not covered regions'''
    def getCoverage(self):
        self.coverage.sort()
        cov = Counter()
        cov[0] = sum([x[1]-x[0] for x in self.notCovered(0)])
        for x in self.coverage:
            cov[x[2]] += x[1] - x[0]
        assert sum([b for b in cov.values()]) == len(self)  # paranoia
        return cov

    '''returns tuples of uncovered regions'''
    def notCovered(self, limit=0):
        self.coverage.sort()
        uncovered = []  # start,end,coverage*len
        if self.coverage:  # has coverage
            # add uncovered start and first coverage segment
            if self.coverage[0][0] > self.coord[0]:  # leading empty segment (no coverage)
                uncovered.append([self.coord[0], self.coverage[0][0], 0])
            # first segement
            if self.coverage[0][2] <= limit:
                uncovered.append([self.coverage[0][0], self.coverage[0][1], (self.coverage[0][1]-self.coverage[0][0])*self.coverage[0][2] ])
            # extend with remaining segment
            if len(self.coverage) > 1:
                for s in [ (self.coverage[x-1], self.coverage[x]) for x in range(1,len(self.coverage)) ]:
                    i,j = s[0],s[1]
                    if s[0][1] != s[1][0]: # gap between covered segments
                        uncovered.append([s[0][1], s[1][0], 0])
                    if s[1][2] <= limit: # if le limit
                        uncovered.append([s[1][0], s[1][1], (s[1][1]-s[1][0])*s[1][2]])
            # add uncovered end
            if self.coverage[-1][1] < self.coord[1]:  # trailing segment (if no coverage)
                uncovered.append([self.coverage[-1][1],self.coord[1],0])
        else:  # no coverage
            uncovered.append([self.coord[0],self.coord[1],0])
        # average segment coverage if limit is not 0
        if limit > 0:
            averaged = uncovered[:1]
            for i in range(1,len(uncovered)):
                if uncovered[i][0] == averaged[-1][1]:  #extend
                    averaged[-1][1] = uncovered[i][1]  # set new end
                    averaged[-1][2] += uncovered[i][2]  # add coverage*len
                else:  #append
                    averaged.append(uncovered[i])
            for av in averaged:
                av[2] = float('{:.2f}'.format(av[2]/float(av[1]-av[0])))
            uncovered = averaged
        return sorted(uncovered)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--threshold', metavar='INT', help='coverage threshold', type=int, default=20)
    parser.add_argument('-s','--steps', nargs='+', help='Coverage Steps', default=[1,10,20,40,80,160,320,640])
    args = parser.parse_args()

    main(args.threshold,map(int,args.steps))
