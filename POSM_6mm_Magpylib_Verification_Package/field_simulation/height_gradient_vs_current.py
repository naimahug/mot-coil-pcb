
from pathlib import Path
import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt
from scipy.optimize import brentq

folder = Path.cwd() / "DENSE_06mm_8Cu_7WINDING_LAYERS"
vertices = np.loadtxt(folder / "vertices_closed_m.txt")
coil = magpy.current.Polyline(current=1.0, vertices=vertices)

currents = np.linspace(0.1, 1.0, 19)
heights_mm = np.linspace(0.5, 20.0, 300)
dz = 1e-5

def bz_gauss(height_mm, current_A):
    point = np.array([[0, 0, height_mm / 1000]])
    return coil.getB(point)[0, 2] * 1e4 * current_A

def gradient_g_cm(height_mm, current_A):
    z = height_mm / 1000
    points = np.array([
        [0, 0, z - dz],
        [0, 0, z + dz]
    ])
    bz = coil.getB(points)[:, 2] * 1e4 * current_A
    return (bz[1] - bz[0]) / (2 * dz * 100)

# Plot 1: zero-point height vs current
zero_heights = []

for current in currents:
    fields = [bz_gauss(h, current) for h in heights_mm]
    zero = np.nan

    for i in range(len(heights_mm) - 1):
        if fields[i] * fields[i + 1] < 0:
            zero = brentq(
                lambda h: bz_gauss(h, current),
                heights_mm[i],
                heights_mm[i + 1]
            )
            break

    zero_heights.append(zero)

plt.figure(figsize=(8, 5))
plt.plot(currents, zero_heights, "o-")
plt.xlabel("Current (A)")
plt.ylabel("Zero-point height (mm)")
plt.title("Magnetic-Field Zero-Point Height vs. Current")
plt.grid(True)
plt.tight_layout()
plt.savefig("height_vs_current.png", dpi=200)
plt.show()

# Plot 2: axial gradient vs current at 6 mm
gradients = [
    gradient_g_cm(6.0, current) for current in currents
]

plt.figure(figsize=(8, 5))
plt.plot(currents, gradients, "o-")
plt.xlabel("Current (A)")
plt.ylabel("Axial gradient dBz/dz (G/cm)")
plt.title("Axial Gradient vs. Current at 6 mm")
plt.grid(True)
plt.tight_layout()
plt.savefig("gradient_vs_current.png", dpi=200)
plt.show()

print("Saved height_vs_current.png")
print("Saved gradient_vs_current.png")
print("\nCurrent (A) | Zero height (mm) | Gradient (G/cm)")
for current, height, grad in zip(currents, zero_heights, gradients):
    print(f"{current:10.2f} | {height:16.4f} | {grad:14.4f}")


