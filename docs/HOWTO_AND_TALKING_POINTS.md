# How to use this, and how to explain it

Part A is a quick-start for running it. Part B is what to say when someone asks "how did you get this spectrum?".
Physics and every VASP tag are in [WORKFLOW_GUIDE.md](WORKFLOW_GUIDE.md).

---------------------------------------------------------------------------------------------------

## Part A. Using the workflow (checklist)

Work in a directory with `00_relax/ 01_displacements/ 02_analysis/ 03_raman/` and this repo's `scripts/` next to them.
Use an environment with phonopy >= 4 (`phonopy-init`), numpy, matplotlib. You supply the POTCAR.

1. **Relax tightly.** Start from `templates/INCAR.relax`. Pass if the OUTCAR says "reached required accuracy" and
   the last TOTAL-FORCE block has max force < 1e-4 eV/A. If VASP stops with `ZBRENT`, copy CONTCAR to POSCAR and
   restart with IBRION=1, ISIF=2, EDIFF=1E-9.
2. **Displacements.** `cp 00_relax/CONTCAR 01_displacements/POSCAR`; put `templates/INCAR.static`, KPOINTS and
   POTCAR in `01_displacements`; run `bash scripts/step1_3_displacements.sh`. Pass if it prints the number of
   displaced cells (10 for YB4).
3. **Forces.** Edit the account/QOS/module lines in `01_displacements/disp_array.sbatch`, set `VASP_BIN`,
   `sbatch disp_array.sbatch`. Pass if every `disp-NNN/OUTCAR` has "General timing".
4. **Collect and analyse.** `bash scripts/step1_3_displacements.sh collect`, then in `02_analysis`:
   `python ../scripts/analyze_gamma_modes.py`. Pass if: three modes at 0 cm^-1, no negative frequencies, the printed
   space group is what you expect, and the irrep counts add up to 3N modes. Outputs: `gamma_modes.csv/txt`,
   `irreps.yaml`, `animations/*.extxyz` (load in OVITO).
5. **Plot.** `python ../scripts/plot_raman_sticks.py --spot1 a.wdf --spot2 b.wdf` overlays measured Renishaw spectra;
   without the flags it draws the computed sticks only. Equal stick heights are not intensities.
6. **Intensities (experimental code, not validated).** In `03_raman`: `python ../scripts/prepare_raman.py --select pilot`,
   submit `eps_array.sbatch`, then `python ../scripts/raman_intensities.py --laser-nm <true value> --drude-gamma <value>`.
   First read the printed `nonlinearity` (should be well below 1) and confirm the Drude term was found (no warning).

If you change the crystal: update the label sets in `analyze_gamma_modes.py` (D4h assumed) and the electron/band
count in `prepare_raman.py` (`--nbands`).

---------------------------------------------------------------------------------------------------

## Part B. Explaining it to a collaborator

### The 30-second version
"We relaxed the structure to essentially zero force with DFT, displaced each symmetry-distinct atom by 0.01 A to
get the force constants, and diagonalised them at the zone centre. That gives every vibration's frequency and
atomic motion. Group theory then tells us which of those vibrations can scatter light, which is the set of
Raman-active peaks, and we compare their positions with the measured spectrum."

### The 3-minute version (say it in this order)
1. **Why DFT phonons?** Raman peaks are vibrations of the crystal. DFT gives the energy surface; its curvature gives
   the vibrational frequencies, with no fitting to the experiment.
2. **Why the zone centre only?** Light carries almost no momentum, so first-order Raman only sees phonons at q = 0.
   There, a one-cell calculation is exact, so no supercell is needed.
3. **How are the force constants obtained?** Displace atoms by 0.01 A (phonopy finds the few independent
   displacements by symmetry), run one static DFT calculation per displacement, and fit forces against
   displacements. Same cutoff, k-mesh and smearing as the relaxation, because otherwise the curvature
   is not consistent with the structure.
4. **What do we get?** For YB4 (20 atoms): 60 modes, three acoustic at zero, none imaginary (which would mean an
   unstable structure). Each mode has a symmetry label.
5. **Which peaks can appear?** The symmetry label decides. In the tetragonal point group D4h, A1g, B1g, B2g and Eg
   modes are Raman active; A2u and Eu are infrared active; the rest are silent. That gave 19 distinct Raman-active
   frequencies (5 A1g, 4 B1g, 4 B2g, 6 Eg).
6. **Where is the motion?** Projecting each eigenvector onto the atoms shows which modes are mostly Y (low frequency,
   heavy atoms) and which are boron-network modes (high frequency).
7. **How do we compare with the experiment?** By peak position and peak count, with the symmetry label attached to each
   computed peak. The computed sticks have equal height: intensities are a separate, harder calculation.

### Likely questions and honest answers
- **"Where are the intensities?"** The phonon calculation gives positions and symmetry selection rules, not strengths.
  Strengths need the Raman tensor (how the dielectric response changes along each mode). We have a frozen-phonon
  pipeline written for it, but it is not yet validated, and this material is a metal, where it is harder. Until it is
  validated, we show positions only.
- **"Why is a metal harder?"** Free carriers screen the field and absorb light, so insulator methods (dielectric
  constants from DFPT, LO-TO splitting) do not apply, and electron-phonon coupling can shift and broaden peaks. A harmonic
  calculation cannot capture those shifts.
- **"Is the structure converged?"** The first phonon pass used a structure with max force 3.7e-4 eV/A and showed no
  imaginary modes. A second pass on a structure at 1.8e-5 eV/A is the one to quote; check that frequencies agree to
  within a few cm^-1.
- **"Why PBE?"** It is the standard functional used for the other results in this project, but boride frequencies can
  depend on it; a PBEsol or r2SCAN check is a sensible robustness test. Do not claim it has been done unless it has.
- **"Could a computed peak be missing from the experiment?"** Yes: equal-height sticks include modes that may be
  very weak, and some active modes may be too weak to see. A missing peak is not evidence against the calculation
  until the intensities are in.
- **"Could the measured spectrum contain peaks that are not in your list?"** Yes: second-order scattering, surface
  layers, impurities or other phases, and non-zero-q defect-activated modes are all outside a first-order Gamma
  calculation.
- **"How would we know the method is right?"** Run the same pipeline on a metal with published Raman data, and compare
  polarized spectra irrep by irrep. Both are planned and neither is done yet.

### What not to claim
- Do not call the black sticks a predicted spectrum; they are predicted peak positions of Raman-active modes.
- Do not quote intensity ratios until the intensity route has been validated and its placeholders (laser line,
  Drude damping, band count, k-mesh) are settled.
- Do not say the assignment of a measured peak to a symmetry label is established from positions alone.
