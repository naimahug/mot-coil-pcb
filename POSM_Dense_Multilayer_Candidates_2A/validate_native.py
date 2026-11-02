"""Run with KiCad's Python/pcbnew environment: python validate_native.py BOARD.
Attaches the included project so its 0.30mm clearance rules and local footprint
library are used. This does not establish thermal or manufacturer approval.
"""
from pathlib import Path
import json,sys
import pcbnew
p=Path(sys.argv[1]).resolve();project=p.with_suffix('.kicad_pro')
sm=pcbnew.GetSettingsManager();assert sm.LoadProject(str(project),True)
b=pcbnew.LoadBoard(str(p));b.SetProject(sm.GetProject(str(project)))
assert b.GetCopperLayerCount()==8
print(p.stem,'loaded;',len(b.GetTracks()),'track/arc/via items;',b.GetCopperLayerCount(),'copper layers',flush=True)
out=p.parent/'native_drc_final.txt';assert pcbnew.WriteDRCReport(b,str(out),pcbnew.EDA_UNITS_MILLIMETRES,False)
report=out.read_text();assert 'Severity: error' not in report,report[:3000]
assert 'Found 0 unconnected pads' in report
print('PASS: zero native DRC errors and zero unconnected pads',flush=True)
