"""Independent Magpylib 5 cross-check for the supplied 6 mm POSM MOT coil.

Run from any directory:
    python verify_magpylib_6mm.py

Requirements:
    python -m pip install 'magpylib>=5,<6' numpy

Uses the bundled routed conductor vertices (metres) and metadata (SI units).
The saved validation.json is the reference from the package's analytic segment solver.
"""
from pathlib import Path
import json
import sys
import numpy as np
import magpylib as magpy

HERE = Path(__file__).resolve().parent
DESIGN = HERE / "DENSE_06mm_8Cu_7WINDING_LAYERS"
VALIDATION = HERE / "validation.json"
META_FILE = DESIGN / "metadata.json"
VERTEX_FILE = DESIGN / "vertices_closed_m.txt"


def main():
    version = getattr(magpy, "__version__", "unknown")
    major = version.split(".")[0]
    print(f"Magpylib version: {version}")
    if major != "5":
        raise RuntimeError(
            f"This script requires Magpylib 5.x (found {version}). Run: "
            "python -m pip install --upgrade 'magpylib>=5,<6' numpy"
        )

    for path in (META_FILE, VERTEX_FILE, VALIDATION):
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path.name}. Keep this script beside validation.json and "
                "the DENSE_06mm_8Cu_7WINDING_LAYERS folder from the supplied ZIP."
            )

    meta = json.loads(META_FILE.read_text())
    validations = json.loads(VALIDATION.read_text())
    ref = next((item for item in validations if item.get("name") == DESIGN.name), None)
    if ref is None:
        raise ValueError(f"No {DESIGN.name} entry found in validation.json")

    vertices = np.loadtxt(VERTEX_FILE, comments="#", dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 2:
        raise ValueError(f"Unexpected vertex array shape: {vertices.shape}")

    current = float(meta["current_20G_A"])
    target_mm = float(meta["target_height_mm"])
    q = np.array([0.0, 0.0, target_mm / 1000.0])  # metres, board origin at F.Cu midplane
    d = 1e-6  # 1 micrometre offset, in metres
    offsets = np.vstack((q, q + np.eye(3) * d, q - np.eye(3) * d))

    coil = magpy.current.Polyline(current=current, vertices=vertices)
    b_tesla = np.asarray(coil.getB(offsets), dtype=float)
    if b_tesla.shape != (7, 3):
        raise RuntimeError(f"Unexpected Magpylib field shape: {b_tesla.shape}")

    # B at q is in tesla; multiply by 1e4 for gauss.
    b_target_g = b_tesla[0] * 1e4
    # Rows are B components; columns are x,y,z derivative directions.
    # dB/dx is T/m, convert to G/cm using 1 T/m = 100 G/cm.
    grad_g_cm = ((b_tesla[1:4] - b_tesla[4:7]) / (2.0 * d)).T * 100.0

    ref_b_g = np.asarray(ref["B_at_target_G"], dtype=float)
    ref_grad = np.asarray(ref["gradient_tensor_G_per_cm"], dtype=float)
    b_diff_g = b_target_g - ref_b_g
    grad_diff = grad_g_cm - ref_grad

    np.set_printoptions(precision=9, suppress=True)
    print("\n--- Input geometry/model ---")
    print(f"Design folder: {DESIGN.name}")
    print(f"Routed vertices loaded: {len(vertices):,}")
    print(f"Current: {current:.12f} A")
    print(f"Target point: (x,y,z) = (0, 0, {target_mm:.6f}) mm")
    print("Model: closed path, including the package's representative local OUT-to-IN return")

    print("\n--- Magpylib result ---")
    print("B at target [G] [Bx, By, Bz]:")
    print(b_target_g)
    print("Gradient tensor [G/cm] (rows Bx,By,Bz; columns x,y,z):")
    print(grad_g_cm)
    print(f"Axial gradient dBz/dz: {grad_g_cm[2, 2]:.9f} G/cm")

    print("\n--- Difference from saved analytic-segment reference ---")
    print("Reference B at target [G]:")
    print(ref_b_g)
    print("Magpylib minus reference B [G]:")
    print(b_diff_g)
    print(f"Maximum absolute field difference: {np.max(np.abs(b_diff_g)):.9g} G")
    print("Reference gradient tensor [G/cm]:")
    print(ref_grad)
    print("Magpylib minus reference gradient tensor [G/cm]:")
    print(grad_diff)
    print(f"Maximum absolute gradient-tensor difference: {np.max(np.abs(grad_diff)):.9g} G/cm")

    # The reference reports a field discrepancy <= 2.50e-16 T, or 2.50e-12 G.
    # Allow a slightly looser threshold for floating-point/platform variation.
    field_tol_g = 1e-8
    grad_tol_g_cm = 1e-3
    field_ok = np.max(np.abs(b_diff_g)) <= field_tol_g
    grad_ok = np.max(np.abs(grad_diff)) <= grad_tol_g_cm
    print("\n--- Pass/fail ---")
    print(f"Field comparison tolerance: {field_tol_g:g} G -> {'PASS' if field_ok else 'CHECK'}")
    print(f"Gradient comparison tolerance: {grad_tol_g_cm:g} G/cm -> {'PASS' if grad_ok else 'CHECK'}")
    print("Saved reference axial gradient:", f"{ref_grad[2, 2]:.9f} G/cm")
    if not (field_ok and grad_ok):
        print("CHECK: compare the Magpylib version, vertex file, current, and target coordinates.")
        sys.exit(2)
    print("PASS: Magpylib independently reproduces the saved analytic-segment field calculation within the stated tolerances.")
    print("Scope: this validates solver agreement for the supplied routed centerline model, not manufactured-PCB performance or thermal behavior.")


if __name__ == "__main__":
    main()
