#!/usr/bin/env python2

__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT"
__version__ = "4.0.7"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"


import os
import re


from crimson import utils, picard



"""Custom metric file parser"""
def parseCustomMetrics(in_data, max_size=1024*1024*1):
    """
    Parses an input Custom Picard-style metrics file into a dictionary.
    Only parses first metrics section!

    :param in_data: Input metrics file.
    :type in_data: str or file handle
    :param max_size: Maximum allowed size of the metrics file (default: 1 MiB).
    :type max_size: int
    :returns: Parsed metrics values.
    :rtype: dict
    """
    # try standard picard parser
    try:
        return picard.parse(in_data,max_size)
    except:
        pass
    # do custom parsing
    with utils.get_handle(in_data) as fh:
        contents = fh.read(max_size)
    sections = contents.strip(os.linesep).split(os.linesep * 2)
    header = [re.compile(r"^#+\s+").sub("", x) for x in sections[0].split(os.linesep)]
    return {
        "header": { "flags": header[1], "time": header[-1] },
        "metrics": picard.parse_metrics(sections[1], os.linesep),
        "histogram": None
    }
