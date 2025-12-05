import sys
import ast
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button

# --- Parsing Logic (Same as before) ---
def extract_balanced_array(text, start_index):
    if start_index >= len(text) or text[start_index] != '[':
        return None, start_index + 1
    open_count = 0
    for i in range(start_index, len(text)):
        if text[i] == '[': open_count += 1
        elif text[i] == ']': open_count -= 1
        if open_count == 0: return text[start_index : i+1], i + 1
    return None, len(text)

def parse_array_string(array_str):
    safe_str = re.sub(r'\bnan\b', '0.0', array_str, flags=re.IGNORECASE)
    safe_str = re.sub(r'\binf\b', '0.0', safe_str, flags=re.IGNORECASE)
    try:
        return np.array(ast.literal_eval(safe_str))
    except:
        return None

def get_data_from_stdin():
    try:
        content = sys.stdin.read()
    except:
        return None
    
    if not content: return None

    pattern = re.compile(r'\[\s*\[\s*\[')
    candidates = []
    
    for match in pattern.finditer(content):
        block_str, _ = extract_balanced_array(content, match.start())
        if not block_str: continue
        arr = parse_array_string(block_str)
        if arr is not None and arr.ndim == 3:
            candidates.append(arr)
            
    if not candidates: return None
    # Return the candidate with the most non-zero elements
    return max(candidates, key=lambda x: np.count_nonzero(np.abs(x) > 1e-6))

# --- Interactive Visualization Logic ---
def main():
    print("Reading and parsing data...")
    data = get_data_from_stdin()
    
    if data is None:
        print("No valid 3D data found. Pipe file content: cat out.txt | python visualize_interactive.py")
        return

    print(f"Data Loaded. Shape: {data.shape}")

    # Initial Setup
    z_dim, y_dim, x_dim = data.shape
    
    # Create the figure and 3D axis
    fig = plt.figure(figsize=(12, 8))
    # Leave room at the bottom for sliders
    ax = fig.add_axes([0.05, 0.25, 0.9, 0.7], projection='3d')
    
    # Global variable to store the scatter plot
    scatter_plot = [None] 

    def update_plot(val=None):
        # 1. Get Slider Values
        x_min, x_max = int(s_xmin.val), int(s_xmax.val)
        y_min, y_max = int(s_ymin.val), int(s_ymax.val)
        z_min, z_max = int(s_zmin.val), int(s_zmax.val)
        threshold = s_thresh.val

        # 2. Clear previous plot
        if scatter_plot[0]:
            scatter_plot[0].remove()
        
        # 3. Filter Data
        # Create a mask for the coordinate cuts
        # We start with all False (masked) or True? 
        # Easier to build indices.
        
        # Get all non-zero indices first to save compute
        z_idx, y_idx, x_idx = np.where(np.abs(data) > threshold)
        
        # Apply coordinate cuts
        mask = (
            (x_idx >= x_min) & (x_idx <= x_max) &
            (y_idx >= y_min) & (y_idx <= y_max) &
            (z_idx >= z_min) & (z_idx <= z_max)
        )
        
        xf, yf, zf = x_idx[mask], y_idx[mask], z_idx[mask]
        values = data[zf, yf, xf]

        # 4. Plot
        if len(values) > 0:
            scatter_plot[0] = ax.scatter(xf, yf, zf, c=values, cmap='viridis', 
                                         marker='s', s=40, edgecolors='k', alpha=0.8)
        else:
            # Plot dummy invisible point to prevent error
            scatter_plot[0] = ax.scatter([0], [0], [0], alpha=0)

        ax.set_xlim(0, x_dim)
        ax.set_ylim(0, y_dim)
        ax.set_zlim(0, z_dim)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        fig.canvas.draw_idle()

    # --- GUI Widget Setup ---
    # Define axes for sliders [left, bottom, width, height]
    axcolor = 'lightgoldenrodyellow'
    
    # X Axis Sliders
    ax_xmin = plt.axes([0.1, 0.15, 0.3, 0.03], facecolor=axcolor)
    ax_xmax = plt.axes([0.1, 0.11, 0.3, 0.03], facecolor=axcolor)
    
    # Y Axis Sliders
    ax_ymin = plt.axes([0.55, 0.15, 0.3, 0.03], facecolor=axcolor)
    ax_ymax = plt.axes([0.55, 0.11, 0.3, 0.03], facecolor=axcolor)

    # Z Axis Sliders
    ax_zmin = plt.axes([0.1, 0.05, 0.3, 0.03], facecolor=axcolor)
    ax_zmax = plt.axes([0.1, 0.01, 0.3, 0.03], facecolor=axcolor)
    
    # Threshold Slider
    ax_thresh = plt.axes([0.55, 0.03, 0.3, 0.03], facecolor=axcolor)

    # Instantiate Sliders
    s_xmin = Slider(ax_xmin, 'X Min', 0, x_dim-1, valinit=0, valstep=1)
    s_xmax = Slider(ax_xmax, 'X Max', 0, x_dim-1, valinit=x_dim-1, valstep=1)
    
    s_ymin = Slider(ax_ymin, 'Y Min', 0, y_dim-1, valinit=0, valstep=1)
    s_ymax = Slider(ax_ymax, 'Y Max', 0, y_dim-1, valinit=y_dim-1, valstep=1)
    
    s_zmin = Slider(ax_zmin, 'Z Min', 0, z_dim-1, valinit=0, valstep=1)
    s_zmax = Slider(ax_zmax, 'Z Max', 0, z_dim-1, valinit=z_dim-1, valstep=1)

    max_val = np.max(data)
    s_thresh = Slider(ax_thresh, 'Min Value', 0.0, max_val, valinit=1e-6)

    # Connect update function
    s_xmin.on_changed(update_plot)
    s_xmax.on_changed(update_plot)
    s_ymin.on_changed(update_plot)
    s_ymax.on_changed(update_plot)
    s_zmin.on_changed(update_plot)
    s_zmax.on_changed(update_plot)
    s_thresh.on_changed(update_plot)

    # Initial plot
    update_plot()
    plt.show()

if __name__ == "__main__":
    main()
