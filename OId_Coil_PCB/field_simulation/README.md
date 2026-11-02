# POSM routed-coil magnetic-field handoff

Generated 2026-10-17 from the four latest dense KiCad candidates, not the superseded wide-track candidates. These are one eight-copper-layer PCB per design. This simulation package does not establish a manufacturable vacuum-compatible stackup.

## Start here

1. Open `field_comparison.png` for the four x-z field sections.
2. Each design folder has `field_plots_20G.png`: |B| contours with projected field streamlines, all three on-axis components, transverse component cuts, and axial gradient. Red crosses mark intended working points, not fitted exact nulls.
3. For Magpylib 5, run `python use_magpylib.py` after `pip install 'magpylib>=5,<6' numpy`.
4. For the older vuthalab `biot-savart==0.3`, use `biot_savart_1A_closed_cm.txt`, `biot_savart_10G_closed_cm.txt`, or `biot_savart_20G_closed_cm.txt`. These match its file parser, but see the numerical caution below.
5. Regenerate everything with `pip install numpy scipy matplotlib sexpdata numba 'magpylib>=5,<6'`, then `python export_and_plot.py`. Its default input is the bundled `source_boards` directory. All paths are relative to the script location. An optional `--source /path/to/mot_dense_family` selects another matching board family.

## Files, units and coordinates

- `biot_savart_*_cm.txt`: header-free comma-separated x,y,z,I; lengths in centimetres, current in amperes. Current on each row flows to the next vertex. Last current value is unused.
- `vertices_*_m.txt`: whitespace-separated x,y,z in metres, with a comment header. Use `numpy.loadtxt` directly with Magpylib 5. Assign current in Python.
- `axis_field_20G.csv`: axis position in mm; B in gauss; axial gradient in G/cm.
- `field_maps_20G.npz`: x_mm, z_mm, B_G[iz,ix,component], current_A. Slice is y=0.
- `metadata.json`: final board geometry, proposed copper depths and drive currents.
- `validation.json`: gradient tensors, numerical refinement checks, independent Magpylib cross-checks.
- `source_boards/`: exact KiCad copper geometry and original final parameter JSON used for extraction; reference copies for reproducibility, not a complete fabrication project.

Origin is aperture centre at the F.Cu copper midplane. +z points toward the intended atoms; +x points toward the terminal tab. +y is the physical board's upper direction, the inverse of native KiCad screen y. Heights are from that copper midplane, not the mirror apex or mirror backside. Layer depths are the proposed stackup; they must change if the fabricator changes the stackup.

Each conductor is followed continuously from IN to OUT. Native copper arcs are sampled with a maximum 2.5-degree step; straight traces are retained. All active winding layers, necks, B-C bridges, outer escapes, backside links and current-carrying via barrel portions are included. Dangling portions of through-via barrels carry no DC series current and are omitted. Finite pad spreading and finite track cross sections are not resolved.

`open` contains the physical PCB current path only. `closed` adds a straight local OUT-to-IN external return at the backside copper midplane. That final line is a representative external wire, **not an extra PCB trace**. Replace it with the actual feed/return cable geometry when mounting is known. Keep the external pair close together to reduce its unwanted field. Do not fabricate or route the artificial closing segment onto the PCB.

At the supplied positive current the axial gradient is about -20 G/cm and transverse gradients about +10 G/cm. This is the same quadrupole strength as the opposite sign convention. Reverse current to reverse every component. At half current, fields and gradients halve, and resistive power quarters. These files do not impose any magnetic compensation bias; small transverse residual fields shift the exact null slightly from the intended axis.

## Validation and interpretation

The field solver uses the analytic finite-straight-segment Biot-Savart expression. At the target and 1 micrometre offsets used for gradient differences, Magpylib 5 agrees to below 1e-12 T. Refining native arc sampling from 2.5 degrees to 0.5 degrees changes the largest tensor scale by less than 0.01% for every board. This is a discretization check, not a bound on fabrication, current, stackup or environmental uncertainty. Map grids use 0.08 mm steps in x and z; centre derivatives use 1 micrometre offsets.

Coils are free-space, infinitesimal centreline filaments. No mirror, chamber, magnetic hardware, ambient fields, eddy currents, thermal change, finite-width current redistribution or compensation coils are included. A centerline model is appropriate for initial fields millimetres away from 0.30 mm traces, but is not a near-copper field/current-density model. Projected streamlines in the y=0 slice are illustrative; the actual routed assembly has a small By component and is not perfectly axisymmetric. Treat the plots as predictions for an ideal proposed coil assembly, not measurements or validation of a complete MOT.

## Why prefer Magpylib over the old 0.3 package here?

The older package's documented input convention matches these TXT files: https://pypi.org/project/biot-savart/0.3/ . Its source has an import of `scipy.integrate.trapz`, which is absent in the modern SciPy used here. A temporary alias allowed inspection and a numerical comparison; this alias is not shipped as a production fix.

Its `slice_coil` and `calculate_field` numerical pipeline is sensitive to segment subdivision for these nonuniform, densely sampled paths. In a three-point check on the 6 mm board, its 0.01 cm slicing produced a maximum component discrepancy of about 0.362 G versus the analytic segment result; bypassing slicing gave about 0.0028 G. The exact comparison is in `legacy_package_check.json`. This is a local benchmark for this input and setting, not a blanket claim about every use of that package. Do not assume a finer `coil_resolution` automatically improves it. Use Magpylib or the included analytic solver for the reported quantitative results.

Magpylib 5 uses SI units: metres, amperes, tesla. Older major versions used different conventions; use the supplied version-5 example. https://magpylib.readthedocs.io/en/stable/_autogen/magpylib.current.html

## Mirror simulation path

For room-temperature aluminum, magnetic permeability is very close to that of free space. Expect a small direct DC magnetization effect, but do not infer a numerical distortion bound for the complete assembly without its alloy, coatings, mounts and geometry. Ferromagnetic screws or nickel-containing parts need separate treatment.

Aluminum conducts electricity, so changing coil current induces eddy currents. These oppose the change and can leave a decaying residual magnetic field after coil current is switched off. PCB-under-mirror placement does not eliminate this. The actual current ramp, mirror alloy conductivity, wall shape and nearby conductors determine the response.

Start with the included full routed free-space model. Then use an axisymmetric ring approximation in FEMM to estimate the mirror's DC and single-frequency effects; FEMM is not a full 3D solver. Use a transient finite-element model such as COMSOL AC/DC for the switch-off waveform and residual gradient. Specify the actual driver current waveform, not merely a nominal switching frequency. Validate mesh refinement, outer air boundary and decay versus time before predicting a cooling-sequence delay. No mirror decay constant has been calculated in this handoff.

Primary documentation: https://www.femm.info/wiki/FAQ and https://www.comsol.com/model/induction-currents-from-circular-coils-124 . COMSOL's conductor example includes aluminum: https://www.comsol.com/model/download/890951/models.acdc.eddy_currents.pdf .

## Substrate and fabrication findings

As checked on 2026-10-17, JLCPCB's published flex capability lists 1, 2 and 4 copper layers and outer copper up to 35 micrometres. The four-layer option lists inner/outer copper up to 1 oz. That does not match these proposed eight-layer, 70-micrometre-every-layer candidates. Its published high-frequency product page lists 2 layers and 1 oz copper. Do not substitute either service's material into this eight-layer design without a new confirmed stackup and recalculation.

- Flex capability: https://jlcpcb.com/capabilities/flex-pcb-capabilities
- High-frequency product: https://jlcpcb.com/pcb-fabrication/high-frequency-pcb

I did not verify a standard JLC eight-layer alumina offering. Solid alumina/AlN substrates and multilayer ceramic circuits require a supplier-specific metallization, via and layer process. Ceramic-filled organic RF laminate is not the same construction as a dense alumina substrate. An RF material designation alone says nothing sufficient about UHV suitability.

For in-vacuum use, select exact laminate, bondply/adhesive, coverlay or mask, finish, solder and cleaning process with the vacuum team and fabricator. Request material-specific outgassing evidence and bake compatibility, and qualify the finished board for the target pressure and temperature. NASA's ASTM E595 database provides total mass loss and collected condensable-material screening; that alone does not certify a complete MOT board for UHV or establish its gas load. A first-order chamber estimate is P_extra = Q_board / S_effective, using measured gas load at operating/bake history and effective pumping speed at the board.

NASA database: https://etd.gsfc.nasa.gov/capabilities/outgassing-database/
ASTM screening scope: https://www.astm.org/e0595-15r21.html

If mechanical placement permits locating the PCB outside vacuum, that avoids introducing its polymer surfaces into the chamber, at the cost of additional coil-to-atom distance and possibly different gradients. Confirm the mechanical section before optimizing the next geometry. Current candidates remain proposed electrical/magnetic geometries, not approved vacuum or thermal assemblies.
