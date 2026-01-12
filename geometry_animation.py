#!/usr/bin/env python3
"""
Manual Geometry Sphere-Orbit Clipping Implementation
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Constants
R_EARTH = 6378137.0

def line_sphere_intersection(p1, p2, sphere_center, radius):
    """Calculate intersection points between line segment and sphere"""
    # Ensure inputs are numpy arrays
    p1 = np.array(p1)
    p2 = np.array(p2) 
    sphere_center = np.array(sphere_center)
    
    # Vector from p1 to p2
    d = p2 - p1
    f = p1 - sphere_center
    
    a = np.dot(d, d)
    b = 2 * np.dot(f, d)
    c = np.dot(f, f) - radius**2
    
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        return None, None  # No intersection
    
    sqrt_discriminant = np.sqrt(discriminant)
    
    if abs(a) < 1e-10:  # Line parallel to sphere center
        return None, None
    
    t1 = (-b - sqrt_discriminant) / (2*a)
    t2 = (-b + sqrt_discriminant) / (2*a)
    
    # Only return intersections within segment bounds [0,1]
    intersections = []
    if 0 <= t1 <= 1:
        intersections.append(p1 + t1 * d)
    if 0 <= t2 <= 1:
        intersections.append(p1 + t2 * d)
    
    if len(intersections) == 2:
        return intersections[0], intersections[1]
    elif len(intersections) == 1:
        return intersections[0], None
    else:
        return None, None
    
    t1 = (-b - sqrt_discriminant) / (2*a)
    t2 = (-b + sqrt_discriminant) / (2*a)
    
    # Only return intersections within segment bounds [0, 1]
    intersections = []
    if 0 <= t1 <= 1:
        intersections.append(p1 + t1 * d)
    if 0 <= t2 <= 1:
        intersections.append(p1 + t2 * d)
    
    if len(intersections) == 2:
        return intersections[0], intersections[1]
    elif len(intersections) == 1:
        return intersections[0], None
    else:
        return None, None

def is_point_visible(point, camera_azim, camera_elev):
    """Check if point is visible from current camera angle"""
    # Convert spherical coordinates to view direction
    theta = np.radians(camera_azim)
    phi = np.radians(90 - camera_elev)
    
    # Camera viewing direction vector
    view_x = np.cos(theta) * np.cos(phi)
    view_y = np.sin(theta) * np.cos(phi)
    view_z = np.sin(phi)
    
    # Vector from sphere center to point
    point_vector = point / np.linalg.norm(point)
    
    # Dot product determines front/back
    dot_product = view_x * point_vector[0] + view_y * point_vector[1] + view_z * point_vector[2]
    
    return dot_product > 0  # Positive = front-facing

def classify_orbit_segment(p1, p2, camera_azim, camera_elev):
    """Classify orbit segment based on visibility"""
    # Check if segment intersects sphere
    # Ensure p1 and p2 are numpy arrays for broadcasting
    p1_array = np.array(p1, dtype=float)
    p2_array = np.array(p2, dtype=float)
    intersect1, intersect2 = line_sphere_intersection(p1_array, p2_array, np.array([0,0,0]), R_EARTH)
    
    if intersect1 is None:
        # No intersection - check midpoint visibility
        midpoint = (p1 + p2) / 2
        if is_point_visible(midpoint, camera_azim, camera_elev):
            return 'visible_front'
        else:
            return 'hidden_back'
    else:
        # Intersects sphere - create visible/hidden segments
        if is_point_visible(p1, camera_azim, camera_elev):
            return 'visible_to_intersection'
        else:
            return 'hidden_to_intersection'

class GeometryClippedAnimator:
    def __init__(self, csv_path):
        self.csv_path = csv_path
        self.df = pd.DataFrame()
        self.fig = None
        self.ax = None
        self.satellite_markers = []
        self.satellite_texts = []
        self.orbit_lines = []
        self.times = []
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
                self.colors[sat] = colors[i]
            
            print(f"Loaded {len(self.satellites)} satellites, {self.n_steps} time steps")
            
        except Exception as e:
            print(f"Error loading data: {e}")
            return
    
    def update_orbit_visibility(self):
        """Update orbit visibility based on current camera angle"""
        camera_azim = getattr(self.ax, 'azim', 45)
        camera_elev = getattr(self.ax, 'elev', 20)
        
        # Clear old orbit lines
        for line in self.orbit_lines:
            line.remove()
        self.orbit_lines = []
        
        # Plot each orbit with proper clipping
        for sat_name in self.satellites:
            sat_data = self.df[self.df['name'] == sat_name].sort_values('time_s')
            color = self.colors[sat_name]
            
            # Get orbit points
            x_vals = sat_data['x'].values
            y_vals = sat_data['y'].values
            z_vals = sat_data['z'].values
            orbit_points = np.array([x_vals, y_vals, z_vals])
            
            # Process each orbit segment
            for i in range(len(orbit_points) - 1):
                p1 = orbit_points[i]
                p2 = orbit_points[i + 1]
                
                # Classify segment visibility
                segment_type = classify_orbit_segment(p1, p2, camera_azim, camera_elev)
                
                if segment_type in ['visible_front', 'visible_to_intersection']:
                    # Draw visible segment
                    self.orbit_lines.append(
                        self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                                  ':', color=color, alpha=0.8, linewidth=2, zorder=1)
                    )
                else:
                    # Don't draw hidden segments
                    pass
    
    def setup_plot(self):
        """Setup 3D plot with geometry clipping"""
        # Create figure with black background
        self.fig = plt.figure(figsize=(12, 10), facecolor='black')
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor('black')
        
        # Remove grids and axes
        self.ax.grid(False)
        self.ax.set_axis_off()
        
        # Set viewing limits
        max_radius = self.df[['x', 'y', 'z']].abs().values.max() * 0.8
        for axis in 'xyz':
            getattr(self.ax, f'set_{axis}lim')([-max_radius, max_radius])
        
        # Plot Earth sphere
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x = R_EARTH * np.outer(np.cos(u), np.sin(v))
        y = R_EARTH * np.outer(np.sin(u), np.sin(v))
        z = R_EARTH * np.outer(np.ones(np.size(u)), np.cos(v))
        
        self.ax.plot_surface(x, y, z, color='white', alpha=0.4, 
                            linewidth=0, antialiased=True, shade=True, zorder=0)
        
        # Initialize orbit lines (will be updated dynamically)
        self.orbit_lines = []
        
        # Create satellite markers
        self.satellite_markers = []
        self.satellite_texts = []
        for sat_name in self.satellites:
            if sat_name.split("_")[0] != "service":
                marker = self.ax.scatter([], [], [],
                                        s=200, c='black', marker='^',
                                        edgecolors='white', linewidths=1, alpha=1.0, zorder=2)
            else:
                marker = self.ax.scatter([], [], [],
                                        s=200, c='black', marker='^',
                                        edgecolors='cyan', linewidths=1, alpha=1.0, zorder=2)
            
            self.satellite_markers.append(marker)
            
            # Create text label
            sat_num = sat_name.split("_")[1] if "_" in sat_name else sat_name[:2]
            text = self.ax.text(0, 0, 0, sat_num,
                              color='white', fontsize=10, fontweight='bold',
                              ha='center', va='center', alpha=1.0, zorder=3)
            self.satellite_texts.append(text)
        
        # Initial orbit visibility setup
        self.update_orbit_visibility()
        
        # Set viewing angle
        self.ax.view_init(elev=20, azim=45)
        
        # Add keyboard controls
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        # Add instructions
        self.fig.text(0.02, 0.98,
                    'Controls: Space=Play/Pause, ←/→=Step, ↑/↓=Jump 10 frames, Q=Quit\\nGeometry: Real sphere clipping',
                    transform=self.fig.transFigure, fontsize=9,
                    verticalalignment='top', color='white',
                    bbox=dict(boxstyle='round', facecolor='black', alpha=0.7))
        
        self.ax.set_title('Geometry-Clipped Satellite Animation', color='white', fontsize=16, pad=20)
    
    def on_key_press(self, event):
        """Handle keyboard events"""
        if event.key == ' ':
            self.is_playing = not self.is_playing
        elif event.key == 'right':
            self.current_frame = (self.current_frame + 1) % self.n_steps
        elif event.key == 'left':
            self.current_frame = (self.current_frame - 1) % self.n_steps
        elif event.key == 'up':
            self.current_frame = (self.current_frame + 10) % self.n_steps
        elif event.key == 'down':
            self.current_frame = (self.current_frame - 10) % self.n_steps
        elif event.key == 'q':
            plt.close('all')
            return
        
        self.update_frame()
    
    def update_frame(self):
        """Update frame with geometry clipping"""
        t = self.times[self.current_frame]
        current_data = self.df[self.df['time_s'] == t].iloc
        
        # Update orbit visibility for current viewing angle
        self.update_orbit_visibility()
        
        # Update satellites
        for i, sat_name in enumerate(self.satellites):
            sat_current = current_data[current_data['name'] == sat_name]
            if len(sat_current) > 0:
                pos_x = sat_current['x'].iloc[0]
                pos_y = sat_current['y'].iloc[0]  
                pos_z = sat_current['z'].iloc[0]
                
                # Update triangle and text
                self.satellite_markers[i]._offsets3d = ([pos_x], [pos_y], [pos_z])
                self.satellite_texts[i].set_position_3d((pos_x, pos_y, pos_z))
        
        # Update title
        status = "PLAYING" if self.is_playing else "PAUSED"
        self.ax.set_title(f'Geometry-Clipped - {status}\\nFrame {self.current_frame}/{self.n_steps-1} | Time = {t:.0f} s',
                         color='white', fontsize=16, pad=20)
        
        self.fig.canvas.draw()
    
    def update(self, frame):
        """Animation update function"""
        if self.is_playing:
            self.current_frame = frame
        self.update_frame()
        return []
    
    def animate(self):
        """Start animation"""
        print("Starting geometry-clipped animation...")
        print("Features: Real sphere-orbit intersection, perspective-based visibility")
        
        self.ani = FuncAnimation(self.fig, self.update, frames=self.n_steps,
                               interval=50, blit=False, repeat=True)
        plt.show()

def main():
    import argparse
    import sys
    import os
    
    parser = argparse.ArgumentParser(description='Geometry-Clipped Satellite Animation')
    parser.add_argument('csv_file', help='Path to satellite CSV data file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.csv_file):
        print(f"Error: File '{args.csv_file}' not found")
        sys.exit(1)
    
    animator = GeometryClippedAnimator(args.csv_file)
    animator.animate()

if __name__ == "__main__":
    main()