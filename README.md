# Raman-Spectra-Phonopy

Gamma-point phonons, symmetry-resolved Raman mode assignment, and Raman intensities for crystals,
using VASP + phonopy (finite displacements) and a frozen-phonon dielectric-function route for intensities.
Worked example: YB4 (P4/mbm, 20 atoms), a metallic boride.

**Status.** Frequencies, symmetry labels, Y/B participation and mode animations are implemented and were run
for YB4. The intensity scripts (`prepare_raman.py`, `raman_intensities.py`) are written and unit-checked on
synthetic input only; they have not yet been run on real VASP output and are not validated. Laser wavelength,
Drude damping and optics convergence (k-mesh, NBANDS, CSHIFT) are placeholders you must set.

## What it does
1. Tight relaxation of the structure (max force < 1e-4 eV/A).
2. `phonopy-init -d --dim "1 1 1"`: displaced cells; one static VASP run each.
3. `phonopy-init -f`: FORCE_SETS; `analyze_gamma_modes.py`: Gamma frequencies and eigenvectors, irreps
   (Mulliken labels), Raman/IR/silent flags, Y vs B participation, extxyz animations for OVITO.
4. `plot_raman_sticks.py`: stick spectrum of Raman-active modes (optionally over a measured `.wdf` spectrum).
5. `prepare_raman.py` + `raman_intensities.py`: Raman tensor by finite differences of eps(omega) along each mode
   (VASP LOPTICS plus a Drude term), powder-averaged intensities.

Full physics, every VASP tag, and the command sequence: [docs/WORKFLOW_GUIDE.md](docs/WORKFLOW_GUIDE.md).
Quick-start checklist and how to explain the method to a collaborator: [docs/HOWTO_AND_TALKING_POINTS.md](docs/HOWTO_AND_TALKING_POINTS.md).

## Layout
- `scripts/` analysis and preparation scripts, SLURM array template
- `templates/` INCAR (relax, static), KPOINTS, relax job script. Supply your own POTCAR.
- `examples/YB4/` tight-relaxed structure, and the mode table and plot from the first-pass phonon run (a slightly looser structure, max force 3.7e-4 eV/A; the tight-structure rerun is pending)
- `docs/` workflow guide

## Requirements
Python with phonopy >= 4 (uses `phonopy-init`), numpy, matplotlib; VASP 6. Set `VASP_BIN` for the SLURM
templates and add your site's account/QOS/module lines. Classification of Raman/IR activity in
`analyze_gamma_modes.py` assumes point group D4h and stops if phonopy returns other labels; edit the label
sets for other crystals.

## Limits
Harmonic, non-resonant first-order Raman. No electron-phonon shifts, resonance or Fano shapes. Metals need
the free-carrier term; its OUTCAR parser has not been tested on a real file.

No license is set yet; choose one before others reuse the code.
