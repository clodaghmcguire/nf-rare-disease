#!/usr/local/bin/env python

import sys
import re
import warnings
import pyhgvs as hgvs
#import pyhgvs.utils as hgvs_utils
from pyhgvs.utils import make_transcript, read_refgene
from pyhgvs.models import Exon,Position, Transcript
from pygr.seqdb import SequenceFileDB
from intervaltree import Interval, IntervalTree
from collections import defaultdict
from itertools import imap

__doc__=="""interconversion of Genomic and HGVS notations"""
__author__ = "David Brawand"
__license__ = "MIT"
__version__ = "4.0.7"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

# interval tree
class IT(dict):
    def __init__(self,fi=None,chrom=0,chromStart=1,chromEnd=2,name=3):
        if fi:
            with open(fi,'r') as fh:
                for line in fh:
                    f = line.split()
                    self.add(f[chrom],f[chromStart],f[chromEnd],f[name])
        return

    def add(self, chrom, chromStart, chromEnd, name=None):
        if chrom not in self.keys():
            self[chrom] = IntervalTree()
        self[chrom].add(Interval(int(chromStart),int(chromEnd),name))

    def get(self,chrom,start,end):
        try:
            self[chrom]
        except:
            return []
        else:
            return [ x.data for x in self[chrom][start:end] ]

    def __contains__(self,f):
        try:
            return bool(self[f[0]][int(f[1]):int(f[2])])
        except KeyError:
            return False
        except:
            raise

    def __str__(self):
        return '\n'.join([ k+":"+str(v) for k,v in self.items()])

    def __repr__(self):
        return '\n'.join([ k+":"+str(v) for k,v in self.items()])


class Babelfish(object):
    def __init__(self,genome,transcripts,extend=500):
        # read genome
        print >> sys.stderr, 'Reading Genome...',
        self.genome = SequenceFileDB(genome)
        # read transcripts
        print >> sys.stderr, '\rReading Annotation...',
        try:
            infile = open(transcripts)
        except:
            raise
        else:  # read transcripts which are in reference
            self.transcripts = [ t for t in list(imap(make_transcript, read_refgene(infile))) \
                if t.tx_position.chrom in self.genome.keys() ]
        # INDEX TRANSCRIPTS AND CREATE INTERVALTREE
        print >> sys.stderr, '\rIndexing Transcripts...',
        self.transcript_index = defaultdict(list)
        self.intervaltree = IT()
        for t in self.transcripts:
            self.transcript_index[t.name].append(t)
            self.transcript_index[t.full_name].append(t)
            pos = t.tx_position
            self.intervaltree.add(pos.chrom, pos.chrom_start-extend, pos.chrom_stop+extend, t)
        print >> sys.stderr, '\rHGVS Babelfish is ready!'

    def genomic2HGVS(self,chrom,offset,ref,alt,transcripts=None,transcriptNames=[]):
        transcripts = [] if transcripts is None else transcripts
        # get transcripts by name
        if transcriptNames:
            for t in transcriptNames:
                transcripts += self.transcript_index.get(t,[])
            transcripts = list(set(transcripts))
        # get transcripts by interval
        if not transcripts:
            transcripts = self.intervaltree.get(chrom,int(offset),int(offset)+max([len(ref),len(alt)]))
        # calc HGVS or each transcript
        results = defaultdict(lambda: 'NA')
        for transcript in transcripts:
            results[transcript.full_name] = hgvs.format_hgvs_name(chrom, int(offset), ref, alt, self.genome, \
                transcript, use_prefix=False, use_gene=False, max_allele_length=50)
        return results

    def interval2HGVS(self,chrom,start,end,transcriptNames=None):
        # grab/find transcripts
        if not transcriptNames:
            transcripts = self.intervaltree.get(chrom,int(start),int(end))
        else:
            transcripts = []
            for t in transcriptNames:
                transcripts += self.transcript_index.get(t)
            transcripts = list(set(transcripts))
        # get transcripts
        results = []
        if not transcripts:  # print genomic_to_cdna_coord
            prefix = 'g.'
            results.append(('genomic',chrom,prefix+str(start+1),prefix+str(end)))
        else:  #
            for tx in transcripts:
                prefix = 'c.' if tx.is_coding else 'n.'
                # check overlap with query
                try:
                    assert tx.tx_position.chrom == chrom
                except:
                    continue
                if tx.strand == '-':
                    txStart = hgvs.genomic_to_cdna_coord(tx,end)  # [from,to] closed interval
                    txEnd = hgvs.genomic_to_cdna_coord(tx,start+1)
                else:
                    txStart = hgvs.genomic_to_cdna_coord(tx,start+1)
                    txEnd = hgvs.genomic_to_cdna_coord(tx,end)  # [from,to] closed interval
                results.append((tx.gene.name,tx.full_name,prefix+str(txStart),prefix+str(txEnd)))
        return results

    def HGVS2genomic(self,query):
        # validate query (snpEff 4.2)
        try:
            assert query.count(':') == 1
        except AssertionError:
            q = query.split(':')
            query = ':'.join([q[0],q[-1]])
        except:
            raise
        try:
            i,j = query.split(':')
            assert i.count('.') < 2  # only one version suffix
        except AssertionError:
            m = re.match('(\w{2}_\d+\.\d+)',i)
            if m:
                query = ':'.join([m.group(1),j])
        except:
            raise

        hgvs_name = hgvs.HGVSName(query)
        # get transcript (or gene)
        transcripts = self.transcript_index.get(hgvs_name.transcript,[]) if hgvs_name.transcript \
            else self.transcript_index.get(hgvs_name.gene,[])
        # get genomic locations
        genomicLoci = []
        for t in transcripts:
            try:
                gl = hgvs.parse_hgvs_name(query, self.genome, transcript=t)
            except KeyError:
                if t.tx_position.chrom not in self.genome.keys():
                    print >> sys.stderr, "WARNING:", t.tx_position.chrom, "not in reference genome ("+t.full_name+")"
                    continue
                raise
            except:
                raise
            else:
                genomicLoci.append(gl)
        return genomicLoci


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('query', nargs='?',help='HGVS/COORDINATE QUERY', default=None)
    parser.add_argument('-g','--genome', metavar='FASTA', help='Genome', default='/home/vagrant/work/genome/human_g1k_v37.fasta')
    parser.add_argument('-a','--annotation', metavar='GENEPRED', help='Annotation', default='/home/vagrant/work/genome/refGene')
    args = parser.parse_args()

    babel = Babelfish(args.genome,args.annotation)

    if args.query:
        import re
        # guess format
        if re.match(r'\w+:\d+-\d+',args.query):  # genomic interval
            m = re.match(r'(\w+):(\d+)-(\d+)',args.query)
            print babel.interval2HGVS(m.group(1),int(m.group(2)),int(m.group(3)))
        elif re.match(r'\w+:\d+:\w+:\w+',args.query):  # variant
            q = args.query.split(':')
            print [ (k,v) for k, v in babel.genomic2HGVS(q[0],int(q[1]),q[2],q[3]).items() ]
        elif re.match(r'.._\d+(\.\d+)?:[cgn]\.',args.query):  # HGVS
            print babel.HGVS2genomic(args.query)
        else:
            print >> sys.stderr, "Could not recognise format, valid arguments are:"
            print >> sys.stderr, "\tCHROM:START:END"
            print >> sys.stderr, "\tCHROM:POS:REF:ALT"
            print >> sys.stderr, "\tTRANSCRIPT:HGVS"
            sys.exit(1)
    else:
        # genomic->HGVS->genomic
        #q = ('11', 17496508, 'T', 'C')
        q = ('4', 144918712, 'C', 'G')
        print "QUERY", q
        hgvs_results = babel.genomic2HGVS(*q)
        for r in hgvs_results.items():
            genomic = babel.HGVS2genomic(':'.join(map(str,r)))
            print '\tRESULT', r
            for g in genomic:
                print '\t\tMAPS', g, q == g

        # genomic interval -> HGVS interval
        for iv in [('16', 226900, 227500), ('11', 5247800, 5247900)]:
            print "QUERY", iv
            for r in babel.interval2HGVS(*iv):
                print "\tRESULT", r
