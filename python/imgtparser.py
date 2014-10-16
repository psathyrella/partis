#!/usr/bin/env python
import csv
import sys
import os
import re
from bs4 import BeautifulSoup

from opener import opener
import utils
# from performanceplotter import PerformancePlotter


html_doc = """
<html><head><title>The Dormouse's story</title></head>
<body>
<p class="title"><b>The Dormouse's story</b></p>

<p class="story">Once upon a time there were three little sisters; and their names were
<a href="http://example.com/elsie" class="sister" id="link1">Elsie</a>,
<a href="http://example.com/lacie" class="sister" id="link2">Lacie</a> and
<a href="http://example.com/tillie" class="sister" id="link3">Tillie</a>;
and they lived at the bottom of a well.</p>

<p class="story">...</p>
"""

class IMGTParser(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, seqfname, infname, datadir):
        self.debug = 0
        n_max_queries = 1
        queries = []

        self.germline_seqs = utils.read_germlines(datadir, remove_N_nukes=False)
        # perfplotter = PerformancePlotter(self.germline_seqs, os.getenv('www') + '/partis/imgt_performance', 'imgt')

        # get info that was passed to joinsolver
        self.seqinfo = {}
        with opener('r')(seqfname) as seqfile:
            reader = csv.DictReader(seqfile)
            iline = 0
            for line in reader:
                if len(queries) > 0 and line['unique_id'] not in queries:
                    continue
                self.seqinfo[line['unique_id']] = line
                iline += 1
                if n_max_queries > 0 and iline >= n_max_queries:
                    break

        with opener('r')(infname) as infile:
            soup = BeautifulSoup(infile)
            for unique_id in self.seqinfo:
                imgtinfo = []
                for pre in soup.find_all('pre'):  # NOTE this loops over everything an awful lot of times. shouldn't really matter for now
                    if unique_id in pre.text:
                        imgtinfo.append(pre.text)
                self.parse_query(unique_id, imgtinfo)

    # ----------------------------------------------------------------------------------------
    def parse_query(self, unique_id, query_info):
        assert len(query_info) == 4  # one for the query sequence, then one for v, d, and j
        query_seq = query_info[0].replace('>', '').replace(unique_id, '')  # strip off the unique id
        query_seq = ''.join(query_seq.split())  # strip of white space
        assert query_seq == self.seqinfo[unique_id]['seq'].lower()
        
                

iparser = IMGTParser('caches/recombinator/simu.csv', '/home/dralph/Dropbox/imgtvquest.html', datadir='./data')
