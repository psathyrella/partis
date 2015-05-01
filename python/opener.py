import bz2
import gzip
import os.path
import sys
import csv
from collections import OrderedDict
from Bio import SeqIO

# ----------------------------------------------------------------------------------------
def opener(mode):
    """
    Open a file, with optional compression based on extension
    """
    exts = {'.bz2': bz2.BZ2File,
            '.gz': gzip.open}
    def open_file(path):
        if not path or ('w' not in mode and not os.path.exists(path)):
            raise Exception('file %s not found' % path)
        if path == '-':
            if mode.startswith('r'):
                return sys.stdin
            else:
                return sys.stdout

        open_fn = exts.get(os.path.splitext(path)[1], open)
        return open_fn(path, mode)
    return open_file
