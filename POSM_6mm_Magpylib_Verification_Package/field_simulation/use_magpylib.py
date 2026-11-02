"""pip install 'magpylib>=5,<6' numpy
Evaluate one exported 1 A coil at a chosen drive current; SI units throughout.
"""
from pathlib import Path
import json
import numpy as np
import magpylib as magpy
if int(magpy.__version__.split('.')[0]) != 5:
    raise RuntimeError('This example targets Magpylib 5 SI units.')
folder=Path(__file__).resolve().parent/'DENSE_06mm_8Cu_7WINDING_LAYERS'
meta=json.loads((folder/'metadata.json').read_text())
vertices=np.loadtxt(folder/'vertices_closed_m.txt')
current=meta['current_20G_A']  # Change to any current; fields scale linearly.
coil=magpy.current.Polyline(current=current,vertices=vertices)
q=np.array([0.,0.,meta['target_height_mm']/1000])
d=1e-6
points=np.vstack([q,q+np.eye(3)*d,q-np.eye(3)*d])
B=coil.getB(points)  # tesla
J=((B[1:4]-B[4:7])/(2*d)).T  # tesla/metre
print('B at target [G]:',B[0]*1e4)
print('Gradient tensor [G/cm]:\n',J*100)
# Never pass a huge grid in one call: memory scales as observers x segments.
# Batch observers in small groups if calculating a map with Magpylib.
