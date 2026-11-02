POSM MOT coil PCB — review draft, not fabrication-approved

Geometry: 60 mm circular outer diameter; 35 mm circular through aperture.
Copper: 8-layer KiCad stack; seven spiral winding layers F.Cu through In6.Cu; B.Cu carries only terminal leadouts.
Trace: 0.30 mm width, 5.15 spiral revolutions per winding layer, ~0.583 mm pitch, ~0.283 mm nominal same-layer spacing.
Spiral radii: 18.95 to 21.95 mm (centerlines); interlayer through vias at radii 18.1 and 22.6 mm.
Edge connector band: nominal radii 23–30 mm (7 mm); two 2.2 mm plated through-hole wire pads at (-15.0, -22.0) and (15.0, -22.0) mm.
Interlayer joins: 6 through-hole vias, with separate entry and exit vias; each consecutive winding reverses radial traversal (in-to-out / out-to-in), while current rotates CCW on all seven layers to make their magnetic fields ADD rather than cancel.

IMPORTANT:
- This is ONE PCB design to manufacture twice. For anti-Helmholtz, mount boards coaxially and wire their current circulation in opposing directions as seen along the same fixed axis. Mirroring/flipping the boards changes the apparent direction; check experimentally before energizing.
- 7 active winding layers requires an 8-copper-layer manufactured board. The eighth B.Cu is not unused copper: it carries the two leadouts.
- Board-file text parentheses and geometry checks were run. KiCad pcbnew/kicad-cli and KiCad DRC were NOT available, so KiCad opening, manufacturing rules, full DRC, shorts, connectivity, thermal and magnetic performance are NOT verified.
- All interlayer vias are standard through vias, which physically pass through all copper layers; the geometry separates the via rings from other-layer spiral rings, but fabrication DRC must still be performed.
- Both pads are assigned one coil net to permit the physical conductor connecting them; verify continuity and lack of unintended connections in KiCad.
- The 35 mm aperture + 7 mm edge band leave only 5.5 mm radial space for winding, so this is about 5 turns per layer, NOT 11 turns.
- The exact turn count is 5.15 geometric revolutions, not exactly five; this staggers via angles and keeps the current continuous.
- Do NOT order boards or run significant current without checking KiCad DRC, copper weight, current/temperature, voltage drop, coil gradient, board separation and driver suitability.
