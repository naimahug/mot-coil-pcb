from pathlib import Path
import json
import numpy as np
import magpylib as magpy
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp
import plotly.graph_objects as go

root = Path(__file__).resolve().parent
folder = root / "DENSE_06mm_8Cu_7WINDING_LAYERS"

# Load the supplied 3D routed current path and metadata.
vertices_m = np.loadtxt(folder / "vertices_closed_m.txt")
vertices_mm = vertices_m * 1000.0
meta = json.loads((folder / "metadata.json").read_text())
current = float(meta["current_20G_A"])

print(f"Loaded {len(vertices_m):,} routed vertices.")
print(f"Current: {current:.6f} A")
print("Building the 3D magnetic-field grid; this may take a few minutes...")

# Downsample only the path used for the visualization field calculation.
# The full-resolution path remains untouched for the independent verification.
# Keep roughly 3,000-5,000 vertices for a reasonable visual calculation time.
stride = max(1, len(vertices_m) // 4000)
vfield = vertices_m[::stride]
if not np.allclose(vfield[-1], vertices_m[-1]):
    vfield = np.vstack([vfield, vertices_m[-1]])

coil = magpy.current.Polyline(current=current, vertices=vfield)

# Set a 3D volume around the real routed geometry.
# Coordinates for the field grid are in millimetres.
lo = vertices_mm.min(axis=0)
hi = vertices_mm.max(axis=0)
span = hi - lo
pad = max(float(np.max(span)) * 0.18, 2.0)
bounds_lo = lo - pad
bounds_hi = hi + pad

# Grid resolution is deliberately moderate to keep this workable on a laptop.
nx, ny, nz = 23, 23, 27
gx = np.linspace(bounds_lo[0], bounds_hi[0], nx)
gy = np.linspace(bounds_lo[1], bounds_hi[1], ny)
gz = np.linspace(bounds_lo[2], bounds_hi[2], nz)

# Magpylib uses metres and returns tesla. Evaluate small batches to limit memory.
grid_points_mm = np.array(
    np.meshgrid(gx, gy, gz, indexing="ij")
).reshape(3, -1).T
grid_points_m = grid_points_mm / 1000.0
field_T = np.empty((len(grid_points_m), 3), dtype=float)

batch = 32
for start in range(0, len(grid_points_m), batch):
    stop = min(start + batch, len(grid_points_m))
    field_T[start:stop] = coil.getB(grid_points_m[start:stop])
    if start == 0 or stop == len(grid_points_m) or start // batch % 30 == 0:
        print(f"  Field grid: {stop:,}/{len(grid_points_m):,} points")

# Convert to gauss and reshape to the 3D grid.
field_G = (field_T * 1e4).reshape(nx, ny, nz, 3)
interps = [
    RegularGridInterpolator(
        (gx, gy, gz), field_G[..., j],
        bounds_error=False, fill_value=np.nan
    )
    for j in range(3)
]

def B_at_mm(point):
    p = np.asarray(point, dtype=float).reshape(1, 3)
    return np.array([f(p)[0] for f in interps])

# Trace field lines through the interpolated 3D vector field.
# Seed them around the coil/aperture region rather than only in one plane.
def trace_line(seed, sign=1, max_length=35.0):
    def rhs(t, r):
        b = B_at_mm(r)
        mag = np.linalg.norm(b)
        if not np.isfinite(mag) or mag < 1e-5:
            return np.zeros(3)
        return sign * b / mag

    def outside(t, r):
        return min(
            r[0] - bounds_lo[0], bounds_hi[0] - r[0],
            r[1] - bounds_lo[1], bounds_hi[1] - r[1],
            r[2] - bounds_lo[2], bounds_hi[2] - r[2]
        )
    outside.terminal = True
    outside.direction = -1

    try:
        sol = solve_ivp(
            rhs, (0, max_length), np.asarray(seed, dtype=float),
            max_step=0.22, rtol=2e-3, atol=1e-4,
            events=outside
        )
        if sol.y.shape[1] >= 8:
            return sol.y.T
    except (ValueError, FloatingPointError):
        pass
    return np.empty((0, 3))

# Start field lines around the aperture and around the winding edges.
rxy = max(float(np.max(np.abs(vertices_mm[:, :2]))), 1.0)
zmid = float(np.median(vertices_mm[:, 2]))
seeds = []
for radius in np.linspace(0.6 * rxy, 1.05 * rxy, 4):
    for angle in np.linspace(0, 2*np.pi, 12, endpoint=False):
        for zoff in (-1.5, 0.0, 1.5):
            seeds.append([
                radius*np.cos(angle),
                radius*np.sin(angle),
                zmid + zoff
            ])

# Also seed near the working region on small circles.
target_z = float(meta["target_height_mm"])
for radius in (0.5, 1.0, 1.5):
    for angle in np.linspace(0, 2*np.pi, 8, endpoint=False):
        seeds.append([
            radius*np.cos(angle), radius*np.sin(angle), target_z
        ])

field_lines = []
for seed in seeds:
    for sign in (1, -1):
        line = trace_line(seed, sign)
        if len(line) > 10:
            field_lines.append(line)

print(f"Traced {len(field_lines)} 3D field-line branches.")

fig = go.Figure()


# Actual routed current path in copper/gold.
fig.add_trace(go.Scatter3d(
    x=vertices_mm[:, 0], y=vertices_mm[:, 1], z=vertices_mm[:, 2],
    mode="lines",
    line=dict(color="#B87333", width=5),
    name="Routed current path",
    hovertemplate="x=%{x:.3f} mm<br>y=%{y:.3f} mm"
                  "<br>z=%{z:.3f} mm<extra></extra>"
))

# Arrows showing the actual current direction through the routed path.
step = max(1, len(vertices_mm) // 28)
idx = np.arange(0, len(vertices_mm)-1, step)
d = vertices_mm[idx+1] - vertices_mm[idx]
norm = np.linalg.norm(d, axis=1)
good = norm > 0
idx, d, norm = idx[good], d[good], norm[good]
d = d / norm[:, None] * 0.65

fig.add_trace(go.Cone(
    x=vertices_mm[idx, 0], y=vertices_mm[idx, 1],
    z=vertices_mm[idx, 2],
    u=d[:, 0], v=d[:, 1], w=d[:, 2],
    anchor="tail", sizemode="absolute", sizeref=0.65,
    colorscale=[[0, "#B87333"], [1, "#B87333"]],
    showscale=False, name="Current direction"
))

# Calculate |B| along every traced field line.
# The interpolators return field components in gauss.
line_magnitudes = []
for line in field_lines:
    components = np.column_stack([f(line) for f in interps])
    line_magnitudes.append(np.linalg.norm(components, axis=1))

finite_values = np.concatenate([
    values[np.isfinite(values)]
    for values in line_magnitudes
    if np.any(np.isfinite(values))
])
bmin = float(np.min(finite_values))
bmax = float(np.max(finite_values))
if bmax <= bmin:
    bmax = bmin + 1e-9

# Rainbow field lines: color represents the calculated |B| in gauss.
for i, (line, bmag) in enumerate(zip(field_lines, line_magnitudes)):
    fig.add_trace(go.Scatter3d(
        x=line[:, 0], y=line[:, 1], z=line[:, 2],
        mode="lines",
        line=dict(
            color=bmag,
            colorscale="Turbo",
            cmin=bmin,
            cmax=bmax,
            width=5
        ),
        name="3D magnetic field lines",
        showlegend=(i == 0),
        hovertemplate="x=%{x:.2f} mm<br>y=%{y:.2f} mm"
                      "<br>z=%{z:.2f} mm<extra></extra>"
    ))

    # One arrow on each line indicates the direction of the magnetic field.
    k = max(2, len(line) // 2)
    tangent = line[min(k+1, len(line)-1)] - line[max(k-1, 0)]
    mag = np.linalg.norm(tangent)
    if mag > 0:
        tangent = tangent / mag * 0.5
        fig.add_trace(go.Cone(
            x=[line[k, 0]], y=[line[k, 1]], z=[line[k, 2]],
            u=[tangent[0]], v=[tangent[1]], w=[tangent[2]],
            anchor="tail", sizemode="absolute", sizeref=0.5,
            colorscale=[[0, "#168BFF"], [1, "#168BFF"]],
            showscale=False, showlegend=False, hoverinfo="skip"
        ))

# Colorbar showing the actual field-magnitude range.
# This invisible marker trace supplies the scale without adding a visible point.
fig.add_trace(go.Scatter3d(
    x=[bounds_hi[0]], y=[bounds_hi[1]], z=[bounds_hi[2]],
    mode="markers",
    marker=dict(
        size=0.1,
        color=[bmin],
        colorscale="Turbo",
        cmin=bmin,
        cmax=bmax,
        showscale=True,
        colorbar=dict(title=dict(text="|B| (G)"), x=1.02),
        opacity=0
    ),
    showlegend=False,
    hoverinfo="skip"
))

# Intended working point.
fig.add_trace(go.Scatter3d(
    x=[0], y=[0], z=[target_z],
    mode="markers+text",
    marker=dict(size=6, color="magenta"),
    text=[f"Target z={target_z:g} mm"],
    textposition="top center",
    name="Intended working point"
))

fig.update_layout(
    title=(
        "POSM 6 mm MOT Coil — 3D Routed Current and Magnetic Field"
        f"<br><sup>I = {current:.4f} A | "
        "3D field lines traced from the routed-coil model</sup>"
    ),
    scene=dict(
        xaxis_title="x (mm)",
        yaxis_title="y (mm)",
        zaxis_title="z (mm)",
        aspectmode="data",
        camera=dict(eye=dict(x=1.6, y=1.6, z=1.2))
    ),
    margin=dict(l=0, r=0, t=90, b=0)
)

out = root / "POSM_6mm_TRUE_3D_field_lines.html"
fig.write_html(str(out), include_plotlyjs=True)
print(f"\nCreated: {out}")
print("Open this HTML file in a browser to rotate and zoom.")
print("Note: field lines use a 3D interpolated grid and a downsampled path for visualization.")
