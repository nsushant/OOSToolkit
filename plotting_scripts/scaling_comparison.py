#!/usr/bin/env python3
"""
Scaling Comparison Plot Generator
Generates two-panel plot comparing local search vs exact method:
1. Panel 1: Execution time vs instance size 
2. Panel 2: Accuracy of local search vs exact method
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import ScalarFormatter
import matplotlib.style as mplstyle

# Use matplotlib default style
mplstyle.use('default')

def load_and_process_data():
    """Load and process scaling comparison data"""
    # Folder of the current script
    this_dir = os.path.dirname(os.path.abspath(__file__))
    # Project root
    root = os.path.dirname(this_dir)
    # Path to data file
    data_path = os.path.join(root, "data", "scaling_comparison_results.csv")
    
    print(f"Loading data from: {data_path}")
    
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data file not found: {data_path}")
    
    df = pd.read_csv(data_path)
    print(f"Loaded data:\n{df}")
    
    # Filter successful runs
    df_success = df[df['success'] == True].copy()
    
    # Separate exact and local search data
    df_exact = df_success[df_success['method'] == 'exact']
    df_local = df_success[df_success['method'] == 'local_search']
    
    # Calculate accuracy metrics
    accuracy_data = []
    for visits in df_exact['visits'].unique():
        exact_row = df_exact[df_exact['visits'] == visits]
        local_row = df_local[df_local['visits'] == visits]
        
        if not exact_row.empty and not local_row.empty:
            exact_dv = exact_row.iloc[0]['deltaV']
            local_dv = local_row.iloc[0]['deltaV']
            quality_gap = local_row.iloc[0]['quality_gap']
            
            accuracy_data.append({
                'visits': visits,
                'quality_gap': quality_gap,
                'exact_deltaV': exact_dv,
                'local_deltaV': local_dv
            })
    
    df_accuracy = pd.DataFrame(accuracy_data)
    
    return df_exact, df_local, df_accuracy

def create_scaling_plot(df_exact, df_local, df_accuracy):
    """Create two-panel comparison plot"""
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: Execution Time vs Instance Size
    ax1.plot(df_exact['visits'], df_exact['time_ms'], 'o-', 
             color='red', linewidth=2, markersize=8, label='Exact Method')
    ax1.plot(df_local['visits'], df_local['time_ms'], 's-', 
             color='blue', linewidth=2, markersize=8, label='Local Search')
    
    ax1.set_xlabel('Number of Satellite Visits', fontsize=12)
    ax1.set_ylabel('Execution Time (ms)', fontsize=12)
    ax1.set_title('Execution Time Scaling Comparison', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Use log scale for better visualization of exponential growth
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim(2, 13)
    
    # Add growth rate annotations
    ax1.text(0.05, 0.95, 'Log-Log Scale\nRed: Exponential growth\nBlue: Polynomial growth', 
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=9)
    
    # Panel 2: Accuracy of Local Search
    ax2.plot(df_accuracy['visits'], df_accuracy['quality_gap'], 'o-', 
             color='green', linewidth=2, markersize=8, label='Local Search Accuracy Gap')
    ax2.fill_between(df_accuracy['visits'], 0, df_accuracy['quality_gap'], 
                    alpha=0.3, color='green')
    
    ax2.set_xlabel('Number of Satellite Visits', fontsize=12)
    ax2.set_ylabel('Quality Gap (%)', fontsize=12)
    ax2.set_title('Local Search Accuracy vs Exact Method', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(2, 13)
    
    # Add horizontal lines for accuracy thresholds
    ax2.axhline(y=1, color='orange', linestyle='--', alpha=0.7, label='1% threshold')
    ax2.axhline(y=5, color='red', linestyle='--', alpha=0.7, label='5% threshold')
    ax2.axhline(y=10, color='purple', linestyle='--', alpha=0.7, label='10% threshold')
    
    # Add accuracy region shading
    ax2.axhspan(0, 1, alpha=0.1, color='green', label='Excellent (<1%)')
    ax2.axhspan(1, 5, alpha=0.1, color='yellow', label='Good (1-5%)')
    ax2.axhspan(5, 10, alpha=0.1, color='orange', label='Acceptable (5-10%)')
    
    # Add statistics text
    if len(df_accuracy) > 0:
        avg_gap = df_accuracy['quality_gap'].mean()
        max_gap = df_accuracy['quality_gap'].max()
        stats_text = f'Avg Gap: {avg_gap:.2f}%\nMax Gap: {max_gap:.2f}%'
        ax2.text(0.05, 0.95, stats_text, 
                 transform=ax2.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5),
                 fontsize=9)
    
    plt.suptitle('Local Search vs Exact Method: Performance Comparison', 
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    return fig

def save_results(fig, output_formats=['pdf']):
    """Save the plot in specified formats"""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(this_dir)
    
    for fmt in output_formats:
        output_path = os.path.join(root, f"scaling_comparison.{fmt}")
        fig.savefig(output_path, format=fmt, dpi=300, bbox_inches='tight')
        print(f"Plot saved as: {output_path}")

def main():
    """Main function to generate scaling comparison plots"""
    try:
        # Load and process data
        df_exact, df_local, df_accuracy = load_and_process_data()
        
        # Create the plot
        fig = create_scaling_plot(df_exact, df_local, df_accuracy)
        
        # Save results
        save_results(fig, ['pdf', 'png'])
        
        # Display summary statistics
        print("\n=== Scaling Comparison Summary ===")
        print(f"Exact method tested on {len(df_exact)} instances")
        print(f"Local search tested on {len(df_local)} instances")
        
        if len(df_accuracy) > 0:
            print(f"Average quality gap: {df_accuracy['quality_gap'].mean():.2f}%")
            print(f"Maximum quality gap: {df_accuracy['quality_gap'].max():.2f}%")
            print(f"Minimum quality gap: {df_accuracy['quality_gap'].min():.2f}%")
            
            # Time complexity comparison
            if len(df_exact) >= 3 and len(df_local) >= 3:
                # Rough growth rate estimation
                exact_times = df_exact['time_ms'].values
                local_times = df_local['time_ms'].values
                visits = df_exact['visits'].values
                
                print(f"\n=== Time Complexity Analysis ===")
                print(f"Exact method: {exact_times[-1]/exact_times[0]:.1f}x slower from {visits[0]} to {visits[-1]} visits")
                print(f"Local search: {local_times[-1]/local_times[0]:.1f}x slower from {visits[0]} to {visits[-1]} visits")
        
        plt.show()
        
    except Exception as e:
        print(f"Error generating plots: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())