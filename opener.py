import bz2
import gzip
import os.path
import sys

#def opener(mode, *args, **kwargs):
#    """
#    Open a file, with optional compression based on extension
#    """
#    exts = {'.bz2': bz2.BZ2File,
#            '.gz': gzip.open}
#    def open_file(path):
#        if path == '-':
#            if mode.startswith('r'):
#                return sys.stdin
#            else:
#                return sys.stdout
#
#        open_fn = exts.get(os.path.splitext(path)[1], open)
#        return open_fn(path, mode, *args, **kwargs)
#    return open_file

def opener(mode):
    """
    Open a file, with optional compression based on extension
    """
    exts = {'.bz2': bz2.BZ2File,
            '.gz': gzip.open}
    def open_file(path):
        if path == '-':
            if mode.startswith('r'):
                return sys.stdin
            else:
                return sys.stdout

        open_fn = exts.get(os.path.splitext(path)[1], open)
        return open_fn(path, mode)
    return open_file
