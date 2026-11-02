# Native and independent validation

Four revised dense multilayer boards were loaded in KiCad 7.0.11 with their included project files. Each project has 0.30 mm minimum copper clearance and track width, 0.50 mm copper-to-edge clearance, 0.60/0.30 mm via pad/drill and the included POSM terminal footprint library.

All four final native reports contain zero DRC errors and zero unconnected pads. Individual native_drc_final.txt reports are included. These checks apply to the proposed design rules, not manufacturer production approval or thermal qualification.

Independent checks verify a continuous series path across every active layer, two external endpoints, all adjacent-turn gaps at least 0.30 mm, no unintended via-to-winding contact on other layers, backside path clearance, backside-to-unrelated-via clearance, at least 0.50 mm clearance to the entire board perimeter and central aperture, and native arc approximation error. Through-via accidental shorts are checked geometrically even though every winding conductor has the same physical net.

The routed current model includes the proposed layer depths, all active pairs, series bridges, interlayer links, used via barrels and a representative local terminal closure. Each candidate reaches 20 G/cm below 2 A in that model. Exact current, gradient, loss, voltage, alignment correction and stackup are recorded in its parameter JSON and stackup CSV.

Native F.Cu exports and generated front/back previews are included for visual review. More than one active copper layer overlaid in the PCB editor will look denser than the single-layer previews.

The baseline reference is unchanged from its prior terminated release. It is not a new two-ampere thermal qualification.
