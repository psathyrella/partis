#!/usr/bin/env python3
"""
Main entry point for partis package
"""

import os
import sys

def main():
    """Main entry point for partis command line interface."""
    # Find the partis installation directory
    # If installed via pip, this will be in site-packages/partis/
    # If in development, it should find the root directory
    current_dir = os.path.dirname(os.path.realpath(__file__))
    
    # Look for the partis root directory
    partis_dir = None
    search_paths = [
        os.path.dirname(current_dir),  # ../
        os.path.dirname(os.path.dirname(current_dir)),  # ../../ (for site-packages/partis/main.py)
    ]
    
    for search_path in search_paths:
        if os.path.exists(os.path.join(search_path, 'bin', 'partis')):
            partis_dir = search_path
            break
    
    if partis_dir is None:
        # Try to use the current directory structure assuming pip install
        partis_dir = os.path.dirname(current_dir)
        if not os.path.exists(os.path.join(partis_dir, 'data')):
            print('Error: Could not locate partis installation directory')
            sys.exit(1)
    
    # Set up environment for partis
    sys.path.insert(0, partis_dir)
    os.environ['PARTIS_DIR'] = partis_dir
    
    # Import and execute the partis main functionality
    try:
        # Try to import the original partis logic
        original_partis = os.path.join(partis_dir, 'bin', 'partis')
        if os.path.exists(original_partis):
            # Read and execute the original script, but set up the environment first
            globals_dict = {'__file__': original_partis}
            with open(original_partis, 'r') as f:
                partis_code = f.read()
            exec(partis_code, globals_dict)
        else:
            print(f'Error: partis script not found at {original_partis}')
            sys.exit(1)
    except Exception as e:
        print(f'Error running partis: {e}')
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()