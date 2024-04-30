#!/usr/bin/env python3

import re
import os
import json
import sys
from zipfile import ZipFile
from crimson import fastqc

from Metrics import MetricsFile
from utils import parseCustomMetrics

output_file = sys.argv[1]
input_files = os.listdir()
print(input_files)

uploadMetrics = []
for f in input_files:
    if f.endswith('fastqc.html'):
        '''FASTQC'''
        f = re.sub(r'\.html$',r'.zip',f)
        if os.path.exists(f):
            with ZipFile(f, 'r') as myzip:
                for zf in myzip.namelist():
                    if zf.endswith('fastqc_data.txt'):
                        tmpFile = myzip.extract(zf)
                        uploadMetrics.append({
                            'type': 'fastqc',
                            'source': f,
                            'data': dict(fastqc.parse(tmpFile))
                        })          
    elif f.endswith('selfSM'):
        '''verifyBamID'''
        m = MetricsFile(f)
        uploadMetrics.append(m.json(type="verifyBamID"))
            
    else:
        '''PICARD'''
        m = re.search('\.(metrics|counts)\.([^\.]+)$',f)
        print(m)
        if m:
            uploadMetrics.append({
                'type': m.group(2),
                'source': f,
                'data': dict(parseCustomMetrics(f)) ## automatically falls back to custom first block parser
                })
                
# write metrics file
with open(output_file,'w') as json_output:
    json.dump(uploadMetrics,json_output,indent=4)
