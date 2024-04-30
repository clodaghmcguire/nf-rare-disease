#!/usr/bin/env python2

import sys
import HGVS

__doc__=="""BED to HGVS intervals"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "4.0.7"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Development"

def readBed(fh):
    for line in fh:
        if not line.startswith('\n'):
            yield line.split()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile', metavar='FILE', help='BED input', type=str)
    parser.add_argument('-g','--genome', metavar='FILE', help='Genome (FASTA)', type=str, required=True)
    parser.add_argument('-a','--annotation', metavar='FILE', help='Genome Annotation (GenePred)', type=str,required=True)
    parser.add_argument('-e','--exclude', metavar='HGVS_PREFIX', help='Exclude comma seperated HGVS prefixes [None]', default=None)
    parser.add_argument('-c','--cleanup', action="store_true", help='Remove non-matching genes', default=False)

    args = parser.parse_args()

    babel = HGVS.Babelfish(args.genome, args.annotation)

    excludehgvs = set(args.exclude.split(',')) if args.exclude else set([])

    # read BED file
    with open(args.bedfile) as fh:
        for iv in readBed(fh):
            original_hgvs_result = hgvs_results = babel.interval2HGVS(iv[0],int(iv[1]),int(iv[2]))
            # remove unwanted genes
            if len(iv)>3 and args.cleanup:
                resultgenes = set([ r[0] for r in hgvs_results ])
                if set(iv[3].split(',')).issubset(resultgenes):  # all query genes have to be in subset
                    hgvs_results = [ r for r in hgvs_results if r[0] in iv[3].split(',') ]
            # print result table
            try:
                assert hgvs_results
            except AssertionError:
                hgvs_results = original_hgvs_result  # fallback
                #raise Exception('EmptyHGVSResult')
            except:
                raise
            else:
                for r in hgvs_results:
                    prefixes = set([r[2][:1], r[3][:1]])
                    if not prefixes.issubset(excludehgvs):
                        print '\t'.join(iv + list(r))
