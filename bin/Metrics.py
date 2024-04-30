#!/usr/bin/env python2

__doc__ = '''
class for parsing of metrics files
default: picard style
other: kraken
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
import datetime
import json
from copy import deepcopy

## import from package or define withing module if called in non-package context
#try:
 #   from . import MetricsFileError
#except ValueError:
 #   class MetricsFileError(Exception):
  #      def __init__(self, message):
   #         super(MetricsFileError, self).__init__(message)

## specific parsers
class MetricsFile(list):
    def __init__(self, fi, noguess=False):
        self.filename = fi
        self.fileheader = []  # file header, if any
        ## get file type from extension
        if fi.endswith('kraken'):  # kraken metrics
            krakenline = re.compile(r"\s*(\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(.*)$")
            self.source = fi[fi.rfind('.')+1:]
            self.append({
                'nrows':6,
                'idxname':{},
                'rownames':['frac','rootedreads','reads','rank','taxonid','name'][::-1],
                'cols':{} })  # initialize
            with open(fi,'r') as fh:
                for i, line in enumerate(fh):
                    if line.startswith('#'):
                        pass  # comment line
                    else:
                        m = krakenline.match(line.rstrip())
                        if m:
                            taxonkey = "_".join([ m.group(4) ] + m.group(6).strip().split()).lower()
                            self[0]['cols'][taxonkey] = list(m.groups())[::-1]  # reverse labels that frac is last row (for intervalmetrics)
                            self[0]['cols'][taxonkey][-1] = str(float(self[0]['cols'][taxonkey][-1])/100)  # convert percentage to fraction
                            self[0]['idxname'][i+1] = taxonkey  # assume first colums (index0) was the rownames
                        else:
                            pass
        else:  # Picard Style Metrics (or simple table with commented header)
            fileheaderbuffer, descriptorbuffer = [], []
            self.source = fi[fi.rfind('.')+1:]
            with open(fi,'r') as fh:
                # parse blocks
                headerfound = False
                for line in fh:
                    if line.startswith('#'):
                        if self.fileheader:  # Fileheader set -> is block description (Picard)
                            descriptorbuffer.append(line.rstrip())
                        elif fileheaderbuffer or line.startswith('##'):
                            fileheaderbuffer.append(line.rstrip())  # store fileheader
                        else:  # assume simple table header
                            try:
                                h = line.split()
                                assert len(h)>1  # more than one column
                            except:
                                raise MetricsFileError('FormatError')
                            else:
                                self.append({'cols': {}, 'idxname': {}, 'nrows':0, 'descriptor': ''})
                                for i, heading in enumerate(h):
                                    heading = heading.strip('#')
                                    self[-1]['cols'][heading] = []
                                    self[-1]['idxname'][i] = heading
                                headerfound = True
                    elif len(line)>1:
                        if not headerfound:  # block header
                            # get header and setup, store descriptor
                            self.append({'cols': {}, 'idxname': {}, 'nrows':0, 'descriptor': deepcopy(descriptorbuffer)})
                            del descriptorbuffer[:]  # empty descriptorbuffer
                            for i, heading in enumerate(line.split()):
                                heading = heading.strip()
                                self[-1]['cols'][heading] = []
                                self[-1]['idxname'][i] = heading
                            headerfound = True
                        else:  # data row
                            self[-1]['nrows'] += 1
                            for i, cell in enumerate(line.split()):
                                self[-1]['cols'][self[-1]['idxname'][i]].append(cell.strip())
                    else: # separator (empty line) increments block
                        headerfound = False  # reset block header
                        if descriptorbuffer:  # store block even if it is empty
                            self.append({'cols': {}, 'idxname': {}, 'nrows':0, 'descriptor': deepcopy(descriptorbuffer)})
                            del descriptorbuffer[:]
                        if fileheaderbuffer:
                            self.fileheader += deepcopy(fileheaderbuffer)
                            del fileheaderbuffer[:]

                # try finding row labels (first column with unique non numeric fields) and default to metric of not found
                for i in range(len(self)):
                    if not noguess:
                        for colnumber in sorted(self[i]['idxname'].keys()):
                            colname = self[i]['idxname'][colnumber]
                            try:
                                multirow = self[i]['nrows'] > 1
                                distinct = len(set(self[i]['cols'][colname])) == len(self[i]['cols'][colname])
                                named = colname.startswith('name')
                                isnumeric = all([ re.match("^(\d+)?\.?(\d+)?$", x) for x in self[i]['cols'][colname][0].split(',')])
                                assert multirow and distinct and (named or not isnumeric)
                            except IndexError:  # no row (empty block)
                                continue
                            except AssertionError:  # not distinct or name or all digit
                                continue
                            except:
                                raise
                            else:
                                self[i]['rownames'] = self[i]['cols'][colname]
                                del self[i]['cols'][colname]  # remove from numerical self
                                del self[i]['idxname'][colnumber]  # remove from col name index
                                break

                    # set default if no labels (only print length of rows as in first column)
                    if 'rownames' not in self[i].keys():
                        # define rowname (whatever precedes metrics in filename or 'metric')
                        m = re.search(r'\.([^\.]+)\.metrics\.',fi)
                        if m:
                            rn = m.group(1)
                        else:
                            rn = 'metric'
                        self[i]['rownames'] = [rn] * self[i]['nrows']
        return

    '''return list of block of metric file'''
    def getTable(self,block,header=False,rownames=False):
        prettylines = []
        try:
            colheaders = [ self[block]['idxname'][c] for c in sorted(self[block]['idxname'].keys()) ]
        except IndexError:
            pass  # dont print anything
        except:
            raise
        else:
            if header and colheaders:
                if rownames:
                    prettylines.append("\t".join(['name']+colheaders))  ## HEADER
                else:
                    prettylines.append("\t".join(colheaders))  ## HEADER
            for r, row in enumerate(self[block]['rownames']):
                if rownames:
                    prettylines.append("\t".join([row] + [ self[block]['cols'][col][r] if col in self[block]['cols'].keys() and r < len(self[block]['cols'][col]) else "NULL" for col in colheaders ]))
                else:
                    prettylines.append("\t".join([ self[block]['cols'][col][r] if col in self[block]['cols'].keys() and r < len(self[block]['cols'][col]) else "NULL" for col in colheaders ]))
        return prettylines

    '''JSON output (for SQVD)'''
    def json(self,block=0,**kwargs):
        source = self.fileheader
        # get metrics dictionary
        data = { k: dict(zip(self[block]['rownames'],v)) for k,v in self[block]['cols'].items() } \
            if self[block]['nrows'] > 1 else { k: v[0] if v else '' for k,v in self[block]['cols'].items() }
        # cast floats if possible
        for k,v in data.items():
            try:
                data[k] = float(v)
            except:
                data[k] = v
        # create data structure (crimson-like)
        datum = {
            "type": self.source,
            "source": self.filename,
            "timestamp": datetime.datetime.utcnow().isoformat(),
            "data": {
                "header": self.fileheader
            }
        }
        # add parsed metrics to source group
        datum["data"][self.source] = { "contents": data }
        # apply overriding properties
        datum.update(kwargs)
        return datum

    '''Prettyprint data'''
    def __repr__(self):
        prettylines = []
        for i, bl in enumerate(self):
            prettylines.append("### BLOCK "+str(i)+' ('+str(bl['nrows'])+' rows)') ## BLOCK
            prettylines += self.getTable(block=i,header=True)
            prettylines.append('')
        return "\n".join(prettylines)

    '''write to file (picard metrics format)'''
    def write(self,fi):
        prettylines = []
        prettylines += self.fileheader
        prettylines.append('')
        for i, bl in enumerate(self):
            prettylines += self[i]['descriptor']
            prettylines += self.getTable(block=i,header=True)
            prettylines.append('')
        with open(fi,'w') as fh:
            for line in prettylines:
                print >> fh, line
        return

    '''check compatibility of MetricsFile'''
    def compatible(self,other):
        # check compatibility
        try:
            # same type
            assert type(self) is type(other)
            # same colums
            assert not any([ len(set(self[i]['idxname'].values()).symmetric_difference(set(other[i]['idxname'].values()))) \
                for i,b in enumerate(self) ])
        except AssertionError:
            return False
        except:
            raise
        else:
            return True

    '''merge two metricsFile objects (checks headers for each block) !atomic'''
    def merge(self,other):
        # require blockwise compatibility (same idxnames in each block)
        try:
            assert self.compatible(other)
        except:
            raise MetricsFileError('IncompatibleFiles')
        # append blocks
        for i,b in enumerate(self):
            # append columns
            for k in b['cols'].keys():
                b['cols'][k] += other[i]['cols'][k]
                assert len(b['cols'][k]) == b['nrows'] + other[i]['nrows']
            # update descriptor (append if different else ignore)
            if any([ b['descriptor'][x]!=other[i]['descriptor'][x] for x in range(len(other[i]['descriptor'])) ]):
                b['descriptor'] += other[i]['descriptor']
            # append rownames
            b['rownames'] += other[i]['rownames']
            assert len(b['rownames']) == b['nrows'] + other[i]['nrows']
            # update nrows (add)
            b['nrows'] += other[i]['nrows']
        # merge header
        self.fileheader += other.fileheader
        return

    '''summarize into transposed tab-delimited file (more readable)'''
    def summarize(self,block=0):
        lines = []
        if len(self)>0:
            ordering = self[block]['idxname']
            fields = self[block]['cols']
            rownames = self[block]['rownames']
            # write summary to file
            for n in sorted(ordering.keys()):
                k = ordering[n]
                v = fields[k]
                if len(v) == 0:
                    continue
                for i,p in enumerate(rownames):
                    try:
                        float(v[i])
                    except ValueError:
                        lines.append('{:<30} {:<16} {:<26} {}'.format(self.source, p, k, v[i]))  # or pass
                    else:
                        lines.append('{:<30} {:<16} {:<26} {}'.format(self.source, p, k, v[i]))
        return lines

    '''get metric (last row of first matching block by default)'''
    def getMetric(self,m,row=-1,block=0):
        for block in range(len(self)):
            try:
                return self[block]['cols'][m][row].split(',')
            except:
                continue
        return '?'  # return nothing if not found


    '''check if metric within bounds (last row of first matching block by default)'''
    def metricWithinInterval(self,m,interval,row=-1):
        for block in range(len(self)):
            try:
                values = self[block]['cols'][m][row].split(',')
            except:
                continue
            else:
                result = []
                for val in values:
                    res = []
                    try:
                        assert len([ i is None or float(i) for i in interval ])==2
                    except ValueError:  # assume factor (not None or float)
                        res.append(val in interval)
                    except:
                        raise
                    else:  # assume numeric interval
                        if interval[0] is not None:
                            if float(interval[0]) <= float(val):
                                res.append(True)
                            else:
                                res.append(False)
                        if interval[1] is not None:
                            if float(val) < float(interval[1]):
                                res.append(True)
                            else:
                                res.append(False)

                    # store result
                    result.append(False if False in res else True)
                return result
        raise ValueError  # metric not found in any block


if __name__=="__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(prog="Metrics.py")
    #   configuration files
    parser.add_argument("metricsfile",metavar="METRICS", help="Picard styled metrics file")
    parser.add_argument('-b', dest="block", type=int, default=0, metavar="BLOCK", help="Metric block number")
    parser.add_argument("-n", dest="nocolnames", default=True, action='store_false', help="suppress column headers")
    parser.add_argument("-r", dest="rownames", default=False, action='store_true', help="print row names")
    parser.add_argument("-d", dest="noguess", default=False, action='store_true', help="do not guess rowname column")
    options = parser.parse_args()
    # parse metric file
    metrics = MetricsFile(options.metricsfile,options.noguess)
    m = metrics.json(somelist=[1,2,3], somedict={"animal":"dog", "count": 2}, somestr='lalala')
    # print json.dumps(m,sort_keys=True,indent=2)
    print >> sys.stdout, '\n'.join(metrics.getTable(block=options.block,header=options.nocolnames,rownames=options.rownames))
