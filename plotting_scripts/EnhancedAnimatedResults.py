import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import argparse
import sys
import os

# Constants
MU_EARTH = 3.986004418e14  # m^3/s^2
R_EARTH = 6378137.0  # m

class SatelliteVisualizer:
    def __init__(self, csv_path, selected_satellites=None):
        self.csv_path = csv_path
        self.selected_satellites = selected_satellites
        self.df = None
        self.fig = None
        self.ax = None
        self.scatter = None
        self.times = None
        self.current_frame = 0
        self.is_paused = False
        self.animation_speed = 1.0
        self.lambert_trajectory = None
        
        self.load_data()
        self.setup_plot()
        
    def load_data(self):
        """Load and filter CSV data"""
        print(f"Loading data from: {self.csv_path}")
        try:
            self.df = pd.read_csv(self.csv_path)
            
            # Filter by selected satellites if specified
            if self.selected_satellites:
                print(f"Filtering to satellites: {self.selected_satellites}")
                self.df = self.df[self.df['name'].isin(self.selected_satellites)]
                if self.df.empty:
                    print(f"Warning: No data found for satellites {self.selected_satellites}")
                    all_sats = list(pd.read_csv(self.csv_path)['name'].unique())
                    print(f"Available satellites: {all_sats}")
                    sys.exit(1)
            
            self.satellites = self.df['name'].unique()
            self.times = sorted(self.df['time_s'].unique())
            self.n_steps = len(self.times)
            
            print(f"Loaded {len(self.satellites)} satellites, {self.n_steps} time steps")
            
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)
    
    def setup_plot(self):
        """Setup the 3D plot"""
        self.fig = plt.figure(figsize=(12, 10))
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # Plot Earth as reference sphere
        u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
        x = R_EARTH*np.cos(u)*np.sin(v)
        y = R_EARTH*np.sin(u)*np.sin(v)
        z = R_EARTH*np.cos(v)
        self.ax.plot_surface(x, y, z, color='lightblue', alpha=0.3, linewidth=0)
        
        # Set equal aspect ratio
        max_radius = self.df[['x', 'y', 'z']].abs().values.max() * 1.1
        for axis in 'xyz':
            getattr(self.ax, f'set_{axis}lim')([-max_radius, max_radius])
        
        # Labels
        self.ax.set_xlabel('X (m)')
        self.ax.set_ylabel('Y (m)')
        self.ax.set_zlabel('Z (m)')
        self.ax.set_title('Satellite Constellation Animation')
        
        # Plot full trajectories for context
        colors = plt.cm.rainbow(np.linspace(0, 1, len(self.satellites)))
        for i, name in enumerate(self.satellites):
            sat_data = self.df[self.df['name'] == name]
            style = '--' if name == "service_1" else '-'
            self.ax.plot(sat_data['x'], sat_data['y'], sat_data['z'], 
                       style, alpha=0.3, linewidth=1, label=name, color=colors[i])
        
        # Create scatter object for current positions
        self.scatter = self.ax.scatter([], [], [], s=50, c='red', marker='o')
        
        # Add legend
        self.ax.legend(loc='upper right', fontsize=8)
        
        # Setup keyboard controls
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        # Instructions
        self.fig.text(0.02, 0.98, 
                    'Controls: Space=Pause/Play, ←/→=Step, +/-=Speed, R=Reset, Q=Quit',
                    transform=self.fig.transFigure, fontsize=10, 
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    def get_satellite_position(self, sat_name, time):
        """Get satellite position at specific time with interpolation"""
        sat_data = self.df[self.df['name'] == sat_name].sort_values('time_s')
        
        if sat_data.empty:
            return None
        
        # Find closest time points
        times = sat_data['time_s'].values
        x = sat_data['x'].values
        y = sat_data['y'].values
        z = sat_data['z'].values
        
        if time <= times[0]:
            return np.array([x[0], y[0], z[0]])
        if time >= times[-1]:
            return np.array([x[-1], y[-1], z[-1]])
        
        # Simple linear interpolation
        idx = np.searchsorted(times, time)
        if idx >= len(times):
            idx = len(times) - 1
        if idx == 0:
            idx = 1
            
        t1, t2 = times[idx-1], times[idx]
        p1 = np.array([x[idx-1], y[idx-1], z[idx-1]])
        p2 = np.array([x[idx], y[idx], z[idx]])
        
        # Linear interpolation
        alpha = (time - t1) / (t2 - t1)
        return p1 * (1 - alpha) + p2 * alpha
    
    def plot_transfer_trajectory(self, sat1, sat2, t1, t2):
        """Plot simple transfer trajectory between two satellites"""
        print(f"Computing transfer trajectory: {sat1} -> {sat2} (t={t1:.0f} to {t2:.0f})")
        
        # Get positions
        r1 = self.get_satellite_position(sat1, t1)
        r2 = self.get_satellite_position(sat2, t2)
        
        if r1 is None or r2 is None:
            print("Could not get satellite positions for transfer trajectory")
            return
        
        # Generate simple elliptical transfer arc
        n_points = 30
        trajectory = []
        
        for i in range(n_points):
            # Simple interpolation for visualization (not physically accurate but good for display)
            alpha = i / (n_points - 1)
            
            # Start with linear interpolation
            pos = r1 * (1 - alpha) + r2 * alpha
            
            # Add some curvature to make it look like an orbit
            # Project to Earth center and add perpendicular component
            to_center = -pos
            to_center_norm = to_center / (np.linalg.norm(to_center) + 1e-10)
            
            # Calculate perpendicular direction
            cross = np.cross(r1, r2)
            if np.linalg.norm(cross) > 1e-10:
                perp = np.cross(cross, pos)
                perp_norm = perp / (np.linalg.norm(perp) + 1e-10)
                
                # Add arc height
                arc_height = np.linalg.norm(pos) * 0.1 * np.sin(alpha * np.pi)
                pos = pos + perp_norm * arc_height
            
            trajectory.append(pos)
        
        trajectory = np.array(trajectory)
        
        # Plot the transfer trajectory
        self.ax.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2], 
                    'r--', linewidth=3, alpha=0.8, label=f'{sat1}→{sat2} Transfer')
        
        # Mark start and end points
        self.ax.scatter(*r1, s=100, c='green', marker='>', label=f'{sat1} at t={t1:.0f}s')
        self.ax.scatter(*r2, s=100, c='orange', marker='s', label=f'{sat2} at t={t2:.0f}s')
        
        print(f"Transfer trajectory plotted with {n_points} points")
        
        return trajectory
    
    def update(self, frame):
        """Animation update function"""
        if not self.is_paused:
            # Apply animation speed
            actual_frame = int(frame * self.animation_speed) % self.n_steps
            self.current_frame = actual_frame
        else:
            actual_frame = self.current_frame
        
        t = self.times[actual_frame]
        current_data = self.df[self.df['time_s'] == t]
        
        if not current_data.empty:
            # Update satellite positions
            self.scatter._offsets3d = (current_data['x'], current_data['y'], current_data['z'])
            
            # Update title with current time and speed
            status = "PAUSED" if self.is_paused else "PLAYING"
            self.ax.set_title(f'Satellite Constellation - {status}\nTime = {t:.0f} s | Speed = {self.animation_speed:.1f}x')
        
        return self.scatter,
    
    def on_key_press(self, event):
        """Handle keyboard events"""
        if event.key == ' ':
            self.is_paused = not self.is_paused
            print(f"Animation {'paused' if self.is_paused else 'resumed'}")
        
        elif event.key == 'right':
            self.current_frame = (self.current_frame + 1) % self.n_steps
            print(f"Time step: {self.times[self.current_frame]:.0f} s")
        
        elif event.key == 'left':
            self.current_frame = (self.current_frame - 1) % self.n_steps
            print(f"Time step: {self.times[self.current_frame]:.0f} s")
        
        elif event.key == '+' or event.key == '=':
            self.animation_speed = min(5.0, self.animation_speed + 0.5)
            print(f"Animation speed: {self.animation_speed:.1f}x")
        
        elif event.key == '-' or event.key == '_':
            self.animation_speed = max(0.1, self.animation_speed - 0.5)
            print(f"Animation speed: {self.animation_speed:.1f}x")
        
        elif event.key == 'r':
            self.current_frame = 0
            self.animation_speed = 1.0
            print("Reset to beginning")
        
        elif event.key == 'q':
            plt.close('all')
            print("Animation stopped")
    
    def animate(self):
        """Start the animation"""
        print("\n=== Animation Controls ===")
        print("Space: Pause/Play")
        print("←/→: Step forward/backward")
        print("+/-: Speed up/slow down")
        print("R: Reset to beginning")
        print("Q: Quit")
        print("="*30)
        
        # Create animation
        self.ani = FuncAnimation(self.fig, self.update, frames=self.n_steps, 
                               interval=50, blit=False, repeat=True)
        
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Enhanced Satellite Constellation Visualizer')
    parser.add_argument('csv_file', help='Path to satellite CSV data file')
    parser.add_argument('--satellites', nargs='+', help='Specific satellites to display')
    parser.add_argument('--transfer', nargs=4, metavar=('SAT1', 'SAT2', 'T1', 'T2'),
                       help='Show transfer trajectory between SAT1 and SAT2 from time T1 to T2')
    parser.add_argument('--list-satellites', action='store_true', 
                       help='List all available satellites in CSV file')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.csv_file):
        print(f"Error: File '{args.csv_file}' not found")
        sys.exit(1)
    
    # List satellites if requested
    if args.list_satellites:
        df = pd.read_csv(args.csv_file)
        satellites = df['name'].unique()
        print(f"Available satellites in {args.csv_file}:")
        for sat in satellites:
            count = len(df[df['name'] == sat])
            print(f"  {sat}: {count} data points")
        sys.exit(0)
    
    # Create visualizer
    viz = SatelliteVisualizer(args.csv_file, args.satellites)
    
    # Add transfer trajectory if specified
    if args.transfer:
        sat1, sat2, t1_str, t2_str = args.transfer
        try:
            t1 = float(t1_str)
            t2 = float(t2_str)
            viz.plot_transfer_trajectory(sat1, sat2, t1, t2)
        except ValueError:
            print("Error: Times must be numeric values")
            sys.exit(1)
    
    # Start animation
    viz.animate()

if __name__ == "__main__":
    main()