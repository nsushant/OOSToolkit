#!/usr/bin/env python3
"""
Custom N-body Satellite Animation - Simplified Debug Version
Features:
- Black background
- Semi-transparent white sphere for Earth (orbits show through)
- No grid lines on sphere
- Simple dotted lines for satellite orbits
- White triangles with satellite numbers
- Arrow key controls for frame-by-frame navigation
- Zoomed in view for better visibility
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import argparse
import sys
import os

# Constants
R_EARTH = 6378137.0  # Earth radius in meters

class CustomSatelliteAnimator:
    def __init__(self, csv_path):
        self.csv_path = csv_path
        self.df = None
        self.fig = None
        self.ax = None
        self.satellite_markers = []
        self.satellite_texts = []
        self.times = None
        self.colors = {}
        self.current_frame = 0
        self.is_playing = True
        self.ani = None
        
        self.load_data()
        self.setup_plot()
        
    def load_data(self):
        """Load CSV data"""
        print(f"Loading data from: {self.csv_path}")
        try:
            self.df = pd.read_csv(self.csv_path)
            self.satellites = self.df['name'].unique()
            self.times = sorted(self.df['time_s'].unique())
            self.n_steps = len(self.times)
            
            # Assign colors to each satellite
            import matplotlib.cm as cm
            color_map = cm.get_cmap('rainbow')
            colors = color_map(np.linspace(0, 1, len(self.satellites)))
            for i, sat in enumerate(self.satellites):
                if(sat.split("_")[0] != "service"):
                    self.colors[sat] = "lime" 
                    
                else:
                    self.colors[sat] = "white" 


            print(f"Loaded {len(self.satellites)} satellites, {self.n_steps} time steps")
            
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)
    
    def setup_plot(self):
        """Setup 3D plot with custom styling"""
        # Create figure with black background
        self.fig = plt.figure(figsize=(12, 10), facecolor='black')
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor('black')
        
        # Remove all grids and axes
        self.ax.grid(False)
        self.ax.set_axis_off()
        
        # Set equal aspect ratio with zoom (smaller radius = more zoomed in)
        max_radius = self.df[['x', 'y', 'z']].abs().values.max() * 0.8  # Zoom in 40%
        for axis in 'xyz':
            getattr(self.ax, f'set_{axis}lim')([-max_radius, max_radius])
        
        # Plot simple dotted orbit lines for each satellite
        for sat_name in self.satellites:
            sat_data = self.df[self.df['name'] == sat_name].sort_values('time_s')
            color = self.colors[sat_name]
            
            # Simple dotted line - just one plot call
            self.ax.plot(sat_data['x'].values, sat_data['y'].values, sat_data['z'].values,
                        ':',  # dotted line style
                        color=color,
                        alpha=0.7,  
                        linewidth=1.5,
                        label=sat_name)
        
        # Create triangular markers for satellites
        self.satellite_markers = []
        self.satellite_texts = []
        for sat_name in self.satellites:
            
            if (sat_name.split("_")[0] != "service"):
 

                marker = self.ax.scatter([], [], [], 
                                        s=300,  # triangle size
                                        c='black',  # white color
                                        marker='^',  # triangular arrow head
                                        edgecolors='white',  # black edge for contrast
                                        linewidths=1.5,
                                        alpha=1.0)

            else:
                

                marker = self.ax.scatter([], [], [], 
                                        s=300,  # triangle size
                                        c='black',  # white color
                                        marker='^',  # triangular arrow head
                                        edgecolors='cyan',  # black edge for contrast
                                        linewidths=1.5,
                                        alpha=1.0)

            self.satellite_markers.append(marker)
            
            # Extract satellite number from name (e.g., "sat_0" -> "0")
            if sat_name.startswith("sat_"):
                sat_num = sat_name.split("_")[1]
            else:
                sat_num = sat_name.split("_")[1]
            
            # Create text label for satellite number
            text = self.ax.text(0, 0, 0, sat_num, 
                              color='white',  # black text for visibility on white triangle
                              fontsize=8, 
                              fontweight='bold',
                              ha='center', va='top',
                              alpha=1.0, 
                                 zorder=1000)
            self.satellite_texts.append(text)
        
        # Plot Earth as completely smooth sphere without any grid
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x = R_EARTH * np.outer(np.cos(u), np.sin(v))
        y = R_EARTH * np.outer(np.sin(u), np.sin(v))
        z = R_EARTH * np.outer(np.ones(np.size(u)), np.cos(v))
        
        # Plot smooth surface
        self.ax.plot_surface(x, y, z, color='white', alpha=0.3, 
                            linewidth=0, antialiased=True, shade=True)
        
        # Set title
        self.ax.set_title('N-body Satellite Animation', color='white', fontsize=16, pad=20)
        
        # Set viewing angle
        self.ax.view_init(elev=20, azim=45)
        
        # Add keyboard controls
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        # Add control instructions
        self.fig.text(0.02, 0.98, 
                    'Controls: Space=Play/Pause, ←/→=Step, ↑/↓=Jump 10 frames, Q=Quit',
                    transform=self.fig.transFigure, fontsize=10, 
                    verticalalignment='top', color='white',
                    bbox=dict(boxstyle='round', facecolor='black', alpha=0.7)
                    )
    
    def on_key_press(self, event):
        """Handle keyboard events"""
        if event.key == ' ':
            self.is_playing = not self.is_playing
            print(f"Animation {'playing' if self.is_playing else 'paused'}")
        
        elif event.key == 'right':
            self.current_frame = (self.current_frame + 1) % self.n_steps
            self.update_frame()
            print(f"Frame {self.current_frame}/{self.n_steps-1} (Time: {self.times[self.current_frame]:.0f}s)")
        
        elif event.key == 'left':
            self.current_frame = (self.current_frame - 1) % self.n_steps
            self.update_frame()
            print(f"Frame {self.current_frame}/{self.n_steps-1} (Time: {self.times[self.current_frame]:.0f}s)")
        
        elif event.key == 'up':
            self.current_frame = (self.current_frame + 10) % self.n_steps
            self.update_frame()
            print(f"Frame {self.current_frame}/{self.n_steps-1} (Time: {self.times[self.current_frame]:.0f}s)")
        
        elif event.key == 'down':
            self.current_frame = (self.current_frame - 10) % self.n_steps
            self.update_frame()
            print(f"Frame {self.current_frame}/{self.n_steps-1} (Time: {self.times[self.current_frame]:.0f}s)")
        
        elif event.key == 'q':
            plt.close('all')
            print("Animation stopped")
    
    def update_frame(self):
        """Update current frame without animation"""
        t = self.times[self.current_frame]
        current_data = self.df[self.df['time_s'] == t]
        
        # Update each satellite marker position
        for i, sat_name in enumerate(self.satellites):
            sat_current = current_data[current_data['name'] == sat_name]
            if len(sat_current) > 0:
                # Get position
                pos_x = sat_current['x'].iloc[0]
                pos_y = sat_current['y'].iloc[0]
                pos_z = sat_current['z'].iloc[0]
                
                # Get velocity for direction
                vel_x = sat_current['vx'].iloc[0]
                vel_y = sat_current['vy'].iloc[0]
                vel_z = sat_current['vz'].iloc[0]
                
                # Update triangle position
                self.satellite_markers[i]._offsets3d = ([pos_x], [pos_y], [pos_z])
                
                # Update text position to follow triangle in 3D
                self.satellite_texts[i].set_position_3d((pos_x, pos_y, pos_z))
                
                # Scale size based on speed
                vel_mag = np.sqrt(vel_x**2 + vel_y**2 + vel_z**2)
                if vel_mag > 0:
                    speed_factor = min(vel_mag / 5000, 2.0)  # Scale size with speed
                    self.satellite_markers[i].set_sizes([200 * speed_factor])
        
        # Update title with current time
        status = "PLAYING" if self.is_playing else "PAUSED"
        self.ax.set_title(f'N-body Satellite Animation - {status}\\nFrame {self.current_frame}/{self.n_steps-1} | Time = {t:.0f} s', 
                         color='white', fontsize=16, pad=20)
        
        self.fig.canvas.draw()
    
    def update(self, frame):
        """Animation update function"""
        if self.is_playing:
            self.current_frame = frame
        
        self.update_frame()
        return self.satellite_markers
    
    def animate(self):
        """Start animation"""
        print("\\n=== Animation Controls ===")
        print("Space: Play/Pause")
        print("←/→: Step forward/backward 1 frame")
        print("↑/↓: Jump forward/backward 10 frames")
        print("Q: Quit")
        print("="*30)
        print(f"Zoomed in 40% to reduce crowding")
        print(f"Starting frame {self.current_frame}/{self.n_steps-1}")
        
        # Create animation
        self.ani = FuncAnimation(self.fig, self.update, frames=self.n_steps, 
                               interval=1, blit=False, repeat=True)
        
        # Show animation
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Custom N-body Satellite Animation')
    parser.add_argument('csv_file', help='Path to satellite CSV data file')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.csv_file):
        print(f"Error: File '{args.csv_file}' not found")
        sys.exit(1)
    
    # Create and run animator
    animator = CustomSatelliteAnimator(args.csv_file)
    animator.animate()

if __name__ == "__main__":
    main()
