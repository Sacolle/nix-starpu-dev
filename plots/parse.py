import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# ------------------------------------------------------------
# Parse the dimensions line (./main VTI sizeX sizeY sizeZ ...)
# ------------------------------------------------------------
def parse_dimensions(text):
    """
    Extracts sizeX, sizeY, sizeZ, absorbWidth from the line that starts with ./main
    """
    line = None
    for l in text.splitlines():
        if l.strip().startswith("./main"):
            line = l.strip()
            break

    if line is None:
        print("ERROR: Could not find './main ...' line in input.")
        return None

    parts = line.split()
    if len(parts) < 6:
        print("ERROR: './main' line does not contain enough parameters.")
        return None

    try:
        sizeX = int(parts[2])
        sizeY = int(parts[3])
        sizeZ = int(parts[4])
        absorbWidth = int(parts[5])
    except:
        print("ERROR: Failed to parse sizeX,sizeY,sizeZ,absorbWidth.")
        return None

    Nx = sizeX + 2 * absorbWidth + 8
    Ny = sizeY + 2 * absorbWidth + 8
    Nz = sizeZ + 2 * absorbWidth + 8

    return Nx, Ny, Nz


# ------------------------------------------------------------
# Extract a 1D stream of numbers (no brackets)
# ------------------------------------------------------------
def extract_number_stream(text):
    """
    Extract all floating-point numbers from the text after the './main' line.
    """
    lines = text.splitlines()
    data_started = False
    numbers = []

    for line in lines:
        if data_started:
            # Extract numbers from the line
            nums = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            numbers.extend(nums)
        if line.strip().startswith("./main"):
            data_started = True  # Start reading numbers after this line

    if not numbers:
        return None

    return np.array(list(map(float, numbers)), dtype=float)


# ------------------------------------------------------------
# Load + reshape
# ------------------------------------------------------------
def load_data():
    content = sys.stdin.read()

    dims = parse_dimensions(content)
    if dims is None:
        return None
    Nx, Ny, Nz = dims

    arr = extract_number_stream(content)
    if arr is None:
        print("ERROR: no numeric data found.")
        return None

    total = Nx * Ny * Nz
    if arr.size != total:
        print(f"ERROR: number count = {arr.size}, expected {total}.")
        return None

    data = arr.reshape((Nz, Ny, Nx))
    return data


# ------------------------------------------------------------
# Visualization — unchanged
# ------------------------------------------------------------
def visualize(data):
    z_dim, y_dim, x_dim = data.shape

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_axes([0.05, 0.25, 0.9, 0.7], projection='3d')

    scatter_plot = [None]

    def update_plot(val=None):
        x_min, x_max = int(s_xmin.val), int(s_xmax.val)
        y_min, y_max = int(s_ymin.val), int(s_ymax.val)
        z_min, z_max = int(s_zmin.val), int(s_zmax.val)
        threshold = s_thresh.val

        if scatter_plot[0]:
            scatter_plot[0].remove()

        z_idx, y_idx, x_idx = np.where(np.abs(data) > threshold)

        mask = (
            (x_idx >= x_min) & (x_idx <= x_max) &
            (y_idx >= y_min) & (y_idx <= y_max) &
            (z_idx >= z_min) & (z_idx <= z_max)
        )

        xf, yf, zf = x_idx[mask], y_idx[mask], z_idx[mask]
        values = data[zf, yf, xf]

        if len(values) > 0:
            scatter_plot[0] = ax.scatter(
                xf, yf, zf, c=values, cmap='viridis',
                marker='s', s=40, edgecolors='k', alpha=0.8
            )
        else:
            scatter_plot[0] = ax.scatter([0], [0], [0], alpha=0)

        ax.set_xlim(0, x_dim)
        ax.set_ylim(0, y_dim)
        ax.set_zlim(0, z_dim)
        fig.canvas.draw_idle()

    # Sliders UI
    axcolor = 'lightgoldenrodyellow'
    ax_xmin = plt.axes([0.1, 0.15, 0.3, 0.03], facecolor=axcolor)
    ax_xmax = plt.axes([0.1, 0.11, 0.3, 0.03], facecolor=axcolor)
    ax_ymin = plt.axes([0.55, 0.15, 0.3, 0.03], facecolor=axcolor)
    ax_ymax = plt.axes([0.55, 0.11, 0.3, 0.03], facecolor=axcolor)
    ax_zmin = plt.axes([0.1, 0.05, 0.3, 0.03], facecolor=axcolor)
    ax_zmax = plt.axes([0.1, 0.01, 0.3, 0.03], facecolor=axcolor)
    ax_thresh = plt.axes([0.55, 0.03, 0.3, 0.03], facecolor=axcolor)

    s_xmin = Slider(ax_xmin, 'X Min', 0, x_dim-1, valinit=0, valstep=1)
    s_xmax = Slider(ax_xmax, 'X Max', 0, x_dim-1, valinit=x_dim-1, valstep=1)
    s_ymin = Slider(ax_ymin, 'Y Min', 0, y_dim-1, valinit=0, valstep=1)
    s_ymax = Slider(ax_ymax, 'Y Max', 0, y_dim-1, valinit=y_dim-1, valstep=1)
    s_zmin = Slider(ax_zmin, 'Z Min', 0, z_dim-1, valinit=0, valstep=1)
    s_zmax = Slider(ax_zmax, 'Z Max', 0, z_dim-1, valinit=z_dim-1, valstep=1)

    s_thresh = Slider(ax_thresh, 'Min Value', 0.0, np.max(np.abs(data)), valinit=1e-6)

    for slider in (s_xmin, s_xmax, s_ymin, s_ymax, s_zmin, s_zmax, s_thresh):
        slider.on_changed(update_plot)

    update_plot()
    plt.show()


def main():
    print("Reading and parsing 1D array...")
    data = load_data()
    if data is None:
        return
    print(f"Loaded data shape = {data.shape}")
    visualize(data)


if __name__ == "__main__":
    main()
