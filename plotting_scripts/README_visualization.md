# Enhanced Satellite Visualization

## Overview

Enhanced Python-based visualization for satellite constellation data with interactive controls and Lambert transfer trajectory capabilities.

## Features

### Core Visualization
- **3D Animation**: Real-time satellite position animation
- **Earth Reference**: Accurate Earth sphere rendering
- **Orbital Trails**: Full trajectory paths with adjustable opacity
- **Satellite Selection**: Filter to display specific satellites only
- **Interactive Controls**: Play/pause, speed control, time scrubbing

### Advanced Features
- **Transfer Trajectories**: Visualize satellite-to-satellite transfers
- **Smooth Animation**: Interpolated motion between CSV time steps
- **Color Coding**: Different colors for different satellites
- **Time Display**: Real-time simulation time shown in title

## Usage

### Basic Usage

```bash
# Activate virtual environment
source .venv/bin/activate

# Animate all satellites
./visualize_satellites.py data/WalkerDelta.csv

# Animate specific satellites only
./visualize_satellites.py data/WalkerDelta.csv --satellites sat_0 sat_3 service_1

# List available satellites
./visualize_satellites.py data/WalkerDelta.csv --list-satellites
```

### Transfer Trajectory Visualization

```bash
# Show transfer trajectory between satellites at specific times
./visualize_satellites.py data/WalkerDelta.csv --transfer sat_0 sat_3 1000 5000

# Combined with satellite selection
./visualize_satellites.py data/WalkerDelta.csv --satellites sat_0 sat_3 service_1 --transfer sat_0 sat_3 1000 5000
```

## Interactive Controls

During animation:

- **Space**: Pause/Resume animation
- **←/→ Arrow Keys**: Step forward/backward through time
- **+/- Keys**: Increase/decrease animation speed (0.1x to 5.0x)
- **R Key**: Reset to beginning
- **Q Key**: Quit animation

## Implementation Notes

### Data Format
Expects CSV with columns:
- `time_s`: Time in seconds
- `name`: Satellite identifier
- `x`, `y`, `z`: Position coordinates (ECI, meters)
- `vx`, `vy`, `vz`: Velocity components (optional, m/s)

### Lambert Transfer Trajectory
The transfer trajectory calculation provides a **visual approximation** of satellite-to-satellite transfers. This is designed for visualization purposes and uses:

- Simple elliptical arc interpolation between two positions
- Perpendicular component for orbital curvature
- Time-based position interpolation for smooth animation

**Note**: This is not a physically accurate Lambert solver - it's optimized for visual clarity and smooth animation. For precise trajectory planning, use the full Lambert solver in the C++ codebase.

### Performance
- **Interpolation**: Linear interpolation between CSV time steps for smooth motion
- **Memory**: Loads entire CSV into memory for fast access
- **Rendering**: Optimized for 60 FPS animation
- **Scalability**: Tested with 30+ satellites over 5000+ time steps

## Examples

### Walker Delta Constellation Analysis
```bash
# Visualize service satellite approach to constellation satellites
./visualize_satellites.py data/WalkerDelta.csv --satellites service_1 sat_0 sat_5 sat_10 --transfer service_1 sat_0 10000 15000

# Analyze specific orbital plane
./visualize_satellites.py data/WalkerDelta.csv --satellites sat_0 sat_1 sat_2 sat_3 sat_4

# Service satellite trajectory planning
./visualize_satellites.py data/WalkerDelta.csv --transfer service_1 sat_15 20000 25000
```

## File Structure

```
Plotting_scripts/
├── AnimatedResults.py              # Original basic animation
├── EnhancedAnimatedResults.py      # Enhanced visualization engine
└── README_visualization.md        # This file

Root/
├── visualize_satellites.py         # User-friendly wrapper script
├── .venv/                       # Python virtual environment
└── data/
    └── WalkerDelta.csv           # Satellite trajectory data
```

## Dependencies

Required Python packages (included in .venv):
- `pandas`: CSV data handling
- `numpy`: Numerical computations and interpolation
- `matplotlib`: 3D plotting and animation
- `mpl_toolkits.mplot3d`: 3D plotting support

## Troubleshooting

### Common Issues

1. **ModuleNotFoundError**: Ensure virtual environment is activated
   ```bash
   source .venv/bin/activate
   ```

2. **File not found**: Check CSV file path
   ```bash
   ./visualize_satellites.py --list-satellites data/WalkerDelta.csv
   ```

3. **No satellite data**: Verify satellite names match CSV
   ```bash
   ./visualize_satellites.py data/WalkerDelta.csv --list-satellites
   ```

### Performance Tips

- For large datasets, consider filtering to specific satellites
- Reduce animation trail opacity for better visibility
- Use smaller windows for better performance on older hardware

## Future Enhancements

Potential improvements:
- **Physically accurate Lambert solver integration**
- **Real-time delta-V calculations**
- **3D model support for satellites**
- **Export capabilities (video, images)**
- **Multiple CSV file comparison**
- **Orbital element display**

## Integration with C++ Codebase

This visualization is designed to complement the existing C++ satellite simulation:

1. **Run simulation**: Generate CSV data using existing `Nbody` executable
2. **Visualize results**: Use this enhanced visualization for analysis
3. **Plan trajectories**: Use transfer visualization for mission planning
4. **Iterate**: Adjust simulation parameters and re-run

The visualization reads the same CSV format generated by the C++ simulation, ensuring seamless integration between simulation and analysis.