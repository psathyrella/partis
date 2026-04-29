#!/usr/bin/env python3
"""
Main entry point for partis package
"""

import runpy
import os
import sys
from pathlib import Path

def main():
    """Main entry point for partis command line interface."""
    bin_partis = str(Path(__file__).resolve().parent.parent / 'bin' / 'partis')
    if not os.path.exists(bin_partis):
        print('Error: partis script not found at %s' % bin_partis)
        sys.exit(1)
    sys.argv[0] = bin_partis
    runpy.run_path(bin_partis, run_name='__main__')

if __name__ == '__main__':
    main()
