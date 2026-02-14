#!/usr/bin/env python3
"""
Enhanced Satellite Constellation Visualizer
Usage examples:
  # Animate all satellites
  ./visualize_satellites.py data/WalkerDelta.csv
  
  # Animate specific satellites
  ./visualize_satellites.py data/WalkerDelta.csv --satellites sat_0 sat_3 service_1
  
  # Show transfer trajectory between satellites
  ./visualize_satellites.py data/WalkerDelta.csv --transfer sat_0 sat_3 1000 5000
  
  # List available satellites
  ./visualize_satellites.py data/WalkerDelta.csv --list-satellites
"""

import sys
import os

# Add plotting scripts to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(script_dir, 'Plotting_scripts'))

from EnhancedAnimatedResults import main

if __name__ == "__main__":
    # Modify sys.argv to handle the wrapper
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    # Pass through arguments to the main function

main()
