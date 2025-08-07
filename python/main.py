#!/usr/bin/env python3
"""
Main entry point for partis package
"""

import os
import sys

def main():
    """Main entry point for partis command line interface."""
    # Set up the partis environment 
    current_dir = os.path.dirname(os.path.realpath(__file__))
    site_packages_dir = os.path.dirname(current_dir)  # site-packages/
    
    # Let the original bin/partis script auto-detect its location
    # Don't override PARTIS_DIR - let it be set naturally
    
    # Execute the original bin/partis script
    bin_partis = os.path.join(site_packages_dir, 'bin', 'partis')
    
    if os.path.exists(bin_partis):
        import subprocess
        result = subprocess.run([sys.executable, bin_partis] + sys.argv[1:])
        sys.exit(result.returncode)
    else:
        print(f'Error: partis script not found at {bin_partis}')
        sys.exit(1)

if __name__ == '__main__':
    main()