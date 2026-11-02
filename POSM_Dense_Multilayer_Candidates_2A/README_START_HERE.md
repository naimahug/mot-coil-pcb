# Revised POSM dense-winding MOT PCB candidates

## What changed

This package supersedes the broad, few-turn, high-current variants. All four new designs use **0.30 mm winding width and at least 0.30 mm clear turn spacing**, 70 um copper, at least ten turns per annulus per winding layer, a 2 mm central aperture and the accepted plated wire terminals with strain relief. The true radial pitch is 0.6008 mm; its 0.0008 mm guard maintains normal gap after curve approximation.

Each candidate is ONE eight-copper-layer PCB. Six or seven layers carry series-connected opposing coil pairs, and B.Cu carries the interlayer series links and terminal routing. This is not an assembly of separate PCBs. Per-layer and total turn counts are distinct: for example, 11 + 10 turns repeated on seven active winding layers means 77 inner and 70 outer turns in the series current path.

The layer count and total thickness are explicit design inputs. These files require an actual multilayer fabrication quote; they must not be ordered as two-layer boards. The 70 um copper target applies to the inner layers as well as the outside copper. Confirm the fabricator's inner-copper capability and 0.30/0.30 mm line/space for this stackup before release. No manufacturer's production acceptance is asserted here.

The original terminated 13 + 13-turn board is included unchanged as REFERENCE_BASELINE. It remains a single winding-plane reference with the previously accepted connections. Its original 0.20 mm width and approximately 0.20 mm clear spacing have not been silently widened. See its reference-current calculation; 2 A is not a newly established thermal rating for that board.

## Gradient and current range

The requested 10-20 G/cm range is the AXIAL gradient. In an approximately axisymmetric quadrupole, the two radial gradients have roughly half the magnitude and opposite sign. The conical-mirror manuscript, arXiv:2610.05983v1, Fig. 4, uses an axial gradient of 17 G/cm for several measurements. That supports this as a useful tuning range, but the experimental optimum still depends on light intensity, detuning, polarization and optical/magnetic alignment.

The hardware screening requirement is 20 G/cm below 2 A, with a nominal search target below 1.8 A to leave current margin. The final routed candidates need only about 0.55-1.01 A for 20 G/cm. A 2 A driver limit is a ceiling, not an instruction to operate every board at 2 A: several would substantially exceed 20 G/cm and their listed dissipation at that current.

For any fixed board and low-field nonmagnetic environment, the field gradient scales linearly with current. The designed zero does not move when the whole series current is changed. Thus set the operating current continuously between the listed I10 and I20 values for a nominal 10-20 G/cm range.

| Target height above F.Cu | Turns per annulus on each active layer | Active winding layers / copper layers | I at 10 G/cm | I at 20 G/cm | Loss at 10 G/cm | Loss at 20 G/cm | Voltage at 20 G/cm | Main body side |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 mm | 11 + 10 | 7 / 8 | 0.276 A | 0.552 A | 0.62 W | 2.48 W | 4.49 V | 45 mm |
| 8 mm | 14 + 11 | 6 / 8 | 0.454 A | 0.907 A | 1.87 W | 7.47 W | 8.23 V | 49 mm |
| 10 mm | 19 + 16 | 7 / 8 | 0.400 A | 0.801 A | 3.48 W | 13.91 W | 17.36 V | 75 mm |
| 12 mm | 22 + 17 | 7 / 8 | 0.505 A | 1.010 A | 6.76 W | 27.03 W | 26.76 V | 85 mm |

The integral 8 x 10 mm terminal tab adds 8 mm to the overall board width. For example, a 45 mm main body has an overall size of 53 x 45 mm. The 10 and 12 mm candidates remain large boards; mechanical fit to the conical assembly or vacuum flange has not been established. The optical aperture itself is 2 mm on every board.

**Power at the upper end is important.** Doubling the gradient from 10 to 20 G/cm doubles current and approximately quadruples resistive loss. The 12 mm design is approximately 6.8 W at 10 G/cm but 27 W at 20 G/cm. It meets the current limit; it is not a low-power thermally qualified design at the upper gradient. Use a shorter working distance if possible, or separately investigate thicker copper, more layers, larger copper volume or active cooling if the full upper range is required there. Do not treat 2 A as a guarantee of low heating.

Losses use copper resistivity 1.724e-8 ohm*m at 20 C. They include all winding traces, local coil links, B.Cu interlayer links and estimated used via barrels with 25 um wall plating. External leads, contacts and the driver are excluded. Hot copper resistance rises approximately by 1 + 0.00393*(T - 20 C). No safe steady-state temperature, continuous current rating or cooling arrangement has been established. Driver voltage must exceed the listed IR voltage with regulation headroom: the tallest candidate needs about 27 V at 20 G/cm before driver and lead overhead.

## Electrical connections

Use only the two external terminals J1.1 (IN) and J1.2 (OUT). Both are 1 mm finished plated holes with 2.2 mm copper pads. There are two 2 mm nonplated strain-relief holes in the integral tab. Set drill and wire details with the actual fabricator and selected conductor bundle.

The first winding layer carries current from its inner via through the inner spiral, across its radial intercoil bridge, through the opposite-handed outer spiral and out to its outer via. B.Cu connects that outer via to the next layer's inner via. This repeats through all active layers; the final outer via goes to OUT. All active pairs have equal current and reinforcing quadrupole orientation. No extra external jumpers between layers are needed.

Through vias are 0.60 mm copper diameter and 0.30 mm finished drill. The inner vias lie on a 1.85 mm-radius circle around the aperture. The winding patterns are rotated between layers, so each pair has its own inner and outer via. A through via deliberately clears all unrelated winding layers. Backside links rotate by one spoke as they travel inward; simply drawing radial lines through other layers' inner vias would short the series chain. Independent checks cover these exact accidental-contact risks, even though all copper has one conductor net.

The winding copper and both terminals use MOT_DISTRIBUTED_WINDING as a physical-conductor net. That permits the intended continuous coil. If integrating a schematic, represent the distributed winding/net tie deliberately; do not overwrite all net assignments blindly from a schematic symbol with unrelated end nets.

Keep IN/OUT wires paired and provide mechanical strain relief. A reference local current closure between terminal centers is included in the field model. The actual remote cable geometry, ambient fields and nearby hardware must be checked in the assembled experiment.

## Stackup and layout inspection

All new boards have eight copper layers and 2.0 mm nominal total thickness. Native Board Setup contains 70 um copper on every layer and an equally spaced dielectric proposal matching the model. Dielectric material is TO_BE_CONFIRMED; the stated dielectric permittivity/loss are placeholders, not a chosen substrate or vacuum approval. Each folder includes proposed_stackup.csv with the modeled copper midplane depths and the purpose of each layer. The fabricator must propose a realizable stackup, and changed layer depths must be put back into the field model before finalizing the radii.

For seven active winding layers, F.Cu and In1.Cu through In6.Cu carry coil pairs; B.Cu handles routing. For six active winding layers, F.Cu and In1.Cu through In5.Cu carry pairs, In6.Cu is empty except through-via pads, and B.Cu handles routing. Do not add a copper pour on the spare layer.

Look at one winding layer at a time in KiCad to see the 0.30 mm clear gap correctly. Overlaying all rotated layers visually fills the gaps even though they are insulated layers. Each folder includes a front/back preview. The comparison image shows the F.Cu pattern only; other active layers have the same turn counts rotated around the center.

The native board uses connected copper TRACK ARCS for the spirals, grouped by layer and coil. They are not DXF graphics or isolated circles. Five-degree arc fits reproduce the intended Archimedean centerline to substantially below one micrometer; exact fit errors and native object counts are recorded in parameters_and_checks.json. Pads, vias and interconnects are native electrical objects.

## Field model and alignment

Heights are measured from the F.Cu midplane to the intended atoms. They are not measured from the mirror's outer face, flange, viewport or nominal apex. The mirror body is about 19.5 mm tall in the supplied paper, but that does not determine the atom-to-PCB distance. Check a real mechanical cross section before selecting a height.

The model includes every active layer at its actual proposed depth, all coil bridges, both terminal escapes, every used via barrel length and all B.Cu interlayer links, closed by a short reference conductor between the external terminals. The outer starting radius is retuned after this complete routed geometry is included. It does not assume all copper planes coincide or simply multiply a single-layer field by the layer count.

Rotating the active layers also reduces transverse field asymmetry. The remaining modeled compensation and uncompensated zero are:

| Height | Small model bias Bx, By, Bz at 20 G/cm | Uncompensated null x, y, z |
|---:|---|---|
| 6 mm | +0.01891 G, +0.00048 G, -0.00026 G | +0.01895 mm, +0.00049 mm, +6.00012 mm |
| 8 mm | +0.02508 G, +0.00124 G, -0.00033 G | +0.02514 mm, +0.00124 mm, +8.00014 mm |
| 10 mm | +0.01544 G, +0.00086 G, -0.00041 G | +0.01545 mm, +0.00087 mm, +10.00020 mm |
| 12 mm | +0.01632 G, +0.00099 G, -0.00049 G | +0.01634 mm, +0.00099 mm, +12.00024 mm |

These small uniform bias corrections would center the model zero exactly on the aperture axis. They scale with current and are separate from ambient-field compensation. No bias-compensation coils are included in this package. The model does not include a full finite-width current-density solution, the actual external wire paths, conductive-mirror switch-off transients or a complete optical MOT model. Calculated heights and gradients are targets to verify experimentally, not guaranteed trapping positions.

More turns increase inductance and can affect field ramp/switch-off. The 10-20 G/cm capture range does not replace the need to turn the field down for optical molasses/PGC. No inductance, flyback clamp, driver dynamics or switch-off time has been qualified here.

## Optimization scope

The search holds width at 0.30 mm, required turn gap at 0.30 mm, copper thickness at 70 um, minimum intercoil edge gap at 1 mm and at least ten turns in each annulus on each active winding layer. It enumerates integer turn counts, selected inner starting radii and layer counts, solving the outer start radius continuously to put the magnetic zero at the target. It screens for I20 <=1.8 A, then selects low estimated copper loss over tested board envelopes. Full routed-field tuning and resistance calculation follow. This is a bounded candidate search, not proof of a global optimum or a fabrication-cost optimum.

The geometry therefore trades multilayer fabrication cost and, for taller targets, board footprint against current and copper loss. It deliberately avoids the previous one-layer wide-trace/few-turn solutions. A smaller permitted board, fewer copper layers, cheaper copper stackup, specified thermal budget or a different actual atom distance can change the appropriate choice.

## Files and verification

- Four DENSE_* folders: native KiCad board, detailed parameters/checks, proposed stackup, front/back preview and native validation outputs.
- REFERENCE_BASELINE: original terminated single-plane design, included as a reference rather than a new 2 A thermal qualification.
- comparison.csv and dense_results.json: calculated operating ranges and physical dimensions.
- dense_comparison.png: visual comparison of one winding layer of each board.
- search_dense.py and build_dense.py: parameter search and complete native geometry source.
- POSM.pretty and per-folder fp-lib-table: the supplied wire-terminal footprint.

Independent checks cover native-file syntax, layer-aware series-path connectivity with only IN and OUT as endpoints, no accidental contact between through vias and unrelated winding layers, winding turn clearance, B.Cu link-to-link and link-to-unrelated-via clearance, central aperture clearance, arc-fit accuracy and a closed-loop field solution. Native KiCad load/export/DRC results are recorded separately. All designs still need final mechanical, substrate, thermal and experiment-level review before manufacturing.

Python dependencies: numpy scipy shapely matplotlib sexpdata. Regenerate with:

    python search_dense.py
    python build_dense.py
    python write_notes.py

Regeneration overwrites output boards; save Naima's edits separately. Changing copper layer depths, turn counts or radii requires rerunning magnetic checks. Changing copper thickness also requires a new resistance and stackup check.

## Primary references

1. Chen et al., Planar-Integrated Magneto-Optical Trap, Physical Review Applied 17, 034031 (2022): https://arxiv.org/abs/2107.07367 . Original planar-annulus geometry and 0.20 mm line/space reference.
2. Nagase et al., Sub-Doppler cooling in a compact conical hollow mirror integrated with anti-Helmholtz coils: https://arxiv.org/abs/2610.05983 . Section II and Figs. 2-4 provide the 2 mm aperture and the example 17 G/cm operating gradient; they do not establish these PCB alternatives' fit.
3. Gao et al., Low-SWaP magneto-optical trap enabled by planar photonic and magnetic components, Science Advances (2026): https://pmc.ncbi.nlm.nih.gov/articles/PMC13488907/ . A multilayer planar-coil precedent; its geometry and reported operating figures are not copied into these candidates.
