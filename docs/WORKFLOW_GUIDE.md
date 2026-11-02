# From a relaxed crystal to a Raman spectrum: what we did and why

Worked example: YB4, 20 atoms in the cell (4 Y, 16 B), space group P4/mbm (#127), point group D4h.
Every step lists the command, the VASP tags (where VASP is involved), the physics, the check that tells
you it worked, and what went wrong for us. Status labels: DONE (ran and checked), BUILT (written and
tested only on synthetic input), PLACEHOLDER (a value I chose that you must replace).

Directory map (all under `YB_exploration/phonon_modes/`):
- `YB4/` first pass: structure from a looser relax. Frequencies and plots live here.
- `YB4_tight/` second pass on the tightly relaxed structure, plus the Raman-intensity machinery.

---------------------------------------------------------------------------------------------------

## 0. The physics in one page

1. **Born-Oppenheimer + harmonic approximation.** Atoms move on the DFT energy surface E(R). Near the minimum,
   E ~ E0 + 1/2 sum Phi_{ia,jb} u_ia u_jb, where Phi = d2E/du_ia du_jb are the force constants (the second
   derivative of the energy with respect to atomic displacements u). First derivatives vanish only at a true
   minimum, which is why Step 1 (tight relaxation) matters.
2. **Normal modes.** Newton's equations give m_i u_i'' = - sum Phi u. For a mode with frequency w the
   displacements are u_i ~ e_i / sqrt(m_i), where e solves the eigenproblem D e = w^2 e with the
   mass-weighted dynamical matrix D_{ia,jb} = Phi_{ia,jb} / sqrt(m_i m_j).
3. **Why only the Gamma point.** Light carries almost no momentum, so first-order Raman scattering only
   creates or absorbs one phonon with wavevector q ~ 0. So we need D(q = 0) only. At q = 0, a 1x1x1 cell is
   exact (a periodic image moves identically to the atom itself), so no supercell is needed.
4. **Finite displacements.** The force constants are measured, not derived: displace one atom by a small
   amount, compute the forces on all atoms with DFT (Hellmann-Feynman), and since F = -Phi u, forces are
   linear in the displacement. Symmetry means only a few displacements are independent.
5. **Group theory (irreps).** Each eigenmode transforms as an irreducible representation of the crystal's
   point group (labels like A1g, Eg). Raman activity is decided by symmetry alone: a mode is Raman active if
   the polarizability tensor component it induces is nonzero. For D4h: A1g, B1g, B2g, Eg are Raman active
   (g = even under inversion); A2u, Eu are infrared active; A1u, A2g, B1u, B2u are silent. Acoustic modes
   (w = 0) belong to A2u + Eu.
6. **Raman intensity (what we do not get from the phonons alone).** Intensity needs the Raman tensor
   alpha_ij = d(chi_ij)/dQ, the change of the electronic polarizability (dielectric response) when atoms move
   along mode Q. Placzek's non-resonant picture gives I ~ |e_s . alpha . e_i|^2 times frequency and thermal
   factors. The phonon calculation gives w and e (so where the peaks are and which are allowed); the extra
   dielectric-function calculations give how strong they are.
7. **Metals.** YB4 is metallic. Free carriers screen and absorb light, so (a) dielectric methods built for
   insulators (LEPSILON, DFPT dielectric constants) are not valid, (b) the dielectric function has an
   interband part (LOPTICS computes it) and an intraband free-carrier (Drude) part that we must add by hand,
   (c) electron-phonon coupling can shift and broaden peaks, which the harmonic calculation cannot capture.

---------------------------------------------------------------------------------------------------

## 1. Tight relaxation of the structure (DONE)

**Why.** Forces on the equilibrium structure must be ~0, otherwise the extracted force constants carry a spurious
linear term and soft modes can appear. Criterion: max force < 1e-4 eV/A.

**First attempt** (`YB4/00_relax`): ISIF=3 (relax ions + cell shape + volume), IBRION=2 (conjugate gradient),
EDIFF=1E-8, EDIFFG=-1E-4. It crashed with `ZBRENT: fatal error in bracketing` at max force 3.7e-4 eV/A. The
energy was already converged to ~1e-7 eV, so the line search hit the noise floor.

**Fix** (`YB4/00_relax_2`, copied into `YB4_tight/00_relax`): restart from the CONTCAR with IBRION=1
(quasi-Newton), ISIF=2 (ions only, because the cell stress was only 0.08 kB), EDIFF=1E-9. Result: converged in
4 ionic steps, max force 1.8e-5 eV/A.

Key INCAR tags used in every VASP step (keep them identical throughout):
| tag | value | meaning |
|---|---|---|
| ENCUT | 910 | plane-wave cutoff in eV; 1.3 x the largest ENMAX (700 eV for the B_h POTCAR) |
| PREC | Accurate | avoids wrap-around errors in the FFT grid |
| ALGO | Normal | robust electronic minimiser (Davidson) |
| ISMEAR / SIGMA | 1 / 0.05 | Methfessel-Paxton smearing, 0.05 eV; required for a metal |
| LREAL | .FALSE. | projectors in reciprocal space, more accurate for forces |
| EDIFF | 1E-8 (static), 1E-9 (relax) | electronic convergence; forces need it tight |
| EDIFFG | -1E-4 | negative = stop when all forces < 1e-4 eV/A |
| IBRION / ISIF / NSW | 2 or 1 / 3 or 2 / 200 | ionic algorithm, what is relaxed, max steps |
| LWAVE / LCHARG | .FALSE. | do not write huge WAVECAR/CHGCAR |

KPOINTS: 9x9x9 Gamma-centred mesh. POTCAR: Y_sv + B_h (PAW_PBE). Functional: PBE.
Check: `grep "reached required accuracy" OUTCAR` and read the last TOTAL-FORCE block.

---------------------------------------------------------------------------------------------------

## 2. Displacements (DONE)

```
conda activate ELFnet            # phonopy 4.4.0 lives here
cp 00_relax/CONTCAR 01_displacements/POSCAR
cd 01_displacements
phonopy-init -d --dim="1 1 1" -c POSCAR
```
- `--dim "1 1 1"`: supercell = unit cell (exact at Gamma, see 0.3).
- Output: `POSCAR-001..010` (displaced cells), `phonopy_disp.yaml` (the bookkeeping file, needed later),
  `SPOSCAR`. Phonopy used amplitude 0.01 A (we checked the yaml) and found 10 symmetry-independent cells for
  20 atoms.
- Gotcha: in phonopy 4, `phonopy -d` was replaced by `phonopy-init -d`; the same for `-f`.

## 3. Forces in every displaced cell (DONE)

One static VASP run per `POSCAR-NNN`, in its own directory `disp-NNN` with the same KPOINTS and POTCAR and:

| tag | value | why |
|---|---|---|
| IBRION / NSW | -1 / 0 | single-point, no relaxation (atoms must stay displaced) |
| EDIFF | 1E-8 | forces converged to far below the 0.01 A displacement signal |
| ADDGRID | .TRUE. | extra support grid, cleaner forces |
| ISMEAR, SIGMA, ENCUT, PREC, LREAL | as in Step 1 | identical to the relaxation, or force constants shift |

Run as a SLURM array (`disp_array.sbatch`, `--array=1-10%4`: ten tasks, four at once, 16 cores each, about
25 min each). The displaced cells lose symmetry, which is why they use 365 k-points against 75 for the relaxed cell.
Check: each `disp-NNN/OUTCAR` contains "General timing"; no `vasprun.xml` truncation.

## 4. Force constants -> Gamma frequencies and eigenvectors (DONE)

```
phonopy-init -f disp-*/vasprun.xml      # writes FORCE_SETS; fails loudly if a run is unfinished
```
`scripts/analyze_gamma_modes.py` then uses the phonopy Python API:
`phonopy.load("phonopy_disp.yaml", force_sets_filename="FORCE_SETS", is_nac=False, produce_fc=True)` builds
Phi by fitting forces to displacements under the crystal symmetry. `run_qpoints([[0,0,0]], with_eigenvectors=True)`
diagonalises D(0). Frequencies are in THz; multiply by 33.35641 for cm^-1.
`is_nac=False`: no non-analytic (LO-TO) correction, because a metal has no long-range macroscopic field.
Result (first pass): 60 modes = 3 x 20 atoms; 3 acoustic at 0; no imaginary frequencies.
Check: three zeros, zero negatives. A negative frequency means an unstable structure or force noise.

## 5. Symmetry labels (DONE)

`ph.set_irreps([0,0,0], degeneracy_tolerance=1e-4)`; labels come from `ph.irreps` (private attribute
`_ir_labels`, because phonopy 4.4 exposes no public accessor for them). Modes whose frequencies agree within the
tolerance are grouped (doublets such as Eg and Eu).
**Our error, caught by the script:** I had assumed space group Pbam (D2h) and wrote D2h Raman rules. phonopy's
labels (Eg, A2u, ...) are tetragonal, and `Symmetry` reported P4/mbm. The script stops if labels fall outside the
expected point group; always check `space group` in the printout before trusting activity flags.
Result: Raman active A1g x5, B1g x4, B2g x4, Eg x6 (12 modes counting both doublet partners) = 19 distinct peaks.
The three acoustic modes are unlabelled by phonopy; we flag them from w ~ 0.

## 6. Atom participation (DONE)

For mode m, the participation of atom i is p_i = sum_xyz |e_{i,xyz}|^2 with e the unit-normalised, mass-weighted
eigenvector (so sum_i p_i = 1). Y% = sum over Y atoms. The dominant Wyckoff site comes from spglib's
`wyckoffs` and `equivalent_atoms` (labels like 4g(Y), 8j(B)). Finding: below ~185 cm^-1 modes are 90-100% Y;
above ~290 cm^-1 they are >97% B.

## 7. Animations (DONE, not yet viewed in OVITO)

Frame k: r_i + A cos(2 pi k / N) u_i, with u_i = e_i / sqrt(m_i), rescaled so the largest atomic amplitude is
A = 0.4 A and N = 20 frames. Eigenvectors are complex in general; we remove a global phase per column and keep
the real part (valid at Gamma, where eigenvectors can be chosen real; the residual imaginary part was 1e-16).
Written as extended XYZ with a per-atom displacement column.

## 8. Comparison with experiment (DONE as an overlay)

`scripts/plot_raman_sticks.py` draws the 19 Raman-active positions as black sticks of equal height over the
measured spectra (Spot 1 blue, Spot 2 orange). Equal heights are not intensities. The experimental `.wdf`
files are Renishaw binary; `scripts/wdf_min.py` is a minimal reader (header at byte 512, tagged blocks, `XLST` =
x-axis, `DATA` = spectra, float32), validated by reproducing Spot 1's text export exactly. We plotted the `_Mod`
files (processed, background removed; check this with your collaborators).

---------------------------------------------------------------------------------------------------

## 9. Raman intensities (BUILT, not yet run on VASP output)

**Frozen-phonon Raman tensor.** For a mode with unit mass-weighted eigenvector e, define the mode coordinate Q
through atomic displacements dr_i = Q e_i / sqrt(m_i) (units A sqrt(amu)). Compute the dielectric function
eps_ij(omega) at structures displaced by +Q and -Q and at the reference; then

    alpha_ij = [eps_ij(+Q) - eps_ij(-Q)] / (2 Q),    evaluated at the laser energy.

The ± pair gives a central difference, cancelling even-order terms. The printout `nonlinearity` =
|eps(+Q)+eps(-Q)-2 eps(0)| / |eps(+Q)-eps(-Q)| must be well below 1.

**VASP tags for eps(omega)** (`prepare_raman.py`, `INCAR_OPTICS`):
| tag | value | meaning |
|---|---|---|
| LOPTICS | .TRUE. | compute the frequency-dependent dielectric function from interband transitions; written to vasprun.xml `<dielectricfunction>` (imag from sum over transitions, real by Kramers-Kronig) |
| NEDOS | 4001 | number of energy grid points for eps(omega) |
| CSHIFT | 0.05 | complex shift (eV) used in the Kramers-Kronig transform, i.e. the broadening of eps |
| NBANDS | 160 | many empty bands, because optics needs transitions to high states (46 occupied bands for 92 electrons). PLACEHOLDER: converge it |
Other tags as in Step 3.

**Free carriers (metals).** LOPTICS covers interband transitions only. The Drude term
eps_D(omega) = - wp^2 / (omega^2 + i gamma omega), where wp is the plasma frequency (VASP prints wp^2 in the
OUTCAR for metals) and gamma a damping, is added in `raman_intensities.py`. The OUTCAR parser was written without
a sample: it has never read a real file (the script warns if it finds nothing). gamma = 0.1 eV is a
PLACEHOLDER.

**Intensity.** For a polycrystal, with a = tr(alpha)/3 (isotropic part) and g^2 the anisotropy
g^2 = 1/2[|a_xx-a_yy|^2 + |a_yy-a_zz|^2 + |a_zz-a_xx|^2 + 6(|a_xy|^2 + |a_yz|^2 + |a_zx|^2)]:

    I_par = 45|a|^2 + 4 g^2,   I_perp = 3 g^2,   I_total = 45|a|^2 + 7 g^2
    I = (omega_L - omega)^4 / omega * (n(omega,T) + 1) * I_total,     n = 1 / (exp(hc omega / k T) - 1)

(Stokes scattering at temperature T; omega in cm^-1; relative intensities only.) Eg partners are summed.
For a single crystal, only the tensor components allowed by each irrep (A1g: xx, yy, zz; B1g: xx-yy; B2g: xy;
Eg: xz, yz) survive; this is the group theory of Step 5 acting on alpha.

**Assumptions that bound the claim** (state them in any paper): non-resonant, first-order, harmonic,
interband-plus-Drude dielectric response; PBE; no electron-phonon shifts or Fano asymmetry; laser wavelength
532 nm is a PLACEHOLDER for the true value; depth sensitivity of a metal surface not modelled.

---------------------------------------------------------------------------------------------------

## 10. Replicating without help: the full command sequence

```
# once
source ~/miniconda3/etc/profile.d/conda.sh; conda activate ELFnet

# 1  tight relax: INCAR above, ISIF=3 first, then IBRION=1/ISIF=2/EDIFF=1E-9 if ZBRENT appears
sbatch runjob                                  # check max force < 1e-4 eV/A

# 2  displacements
cd 01_displacements && cp ../00_relax/CONTCAR POSCAR
phonopy-init -d --dim="1 1 1" -c POSCAR
# make disp-NNN/ with POSCAR-NNN as POSCAR, plus INCAR.static, KPOINTS, POTCAR
sbatch disp_array.sbatch                       # one VASP run per cell

# 3  force sets and analysis
phonopy-init -f disp-*/vasprun.xml
cd ../02_analysis && cp ../01_displacements/{FORCE_SETS,phonopy_disp.yaml} .
python ../scripts/analyze_gamma_modes.py       # table + animations
python ../scripts/plot_raman_sticks.py         # overlay with experiment

# 4  intensities
cd ../03_raman
python ../scripts/prepare_raman.py --select pilot   # then --select all
sbatch eps_array.sbatch
python ../scripts/raman_intensities.py --laser-nm <real value> --temp 300 --drude-gamma <fit>
```
(`scripts/step1_3_displacements.sh` automates the displacement and collect steps in 2 and 3.)

## 11. What still has to be true before this is publishable

1. Frequencies: tight-structure rerun agrees with the first pass within a few cm^-1; cross-check with DFPT (IBRION=8).
2. Intensities: converge k-mesh, NBANDS, CSHIFT and displacement size on one mode; confirm the OUTCAR plasma
   parser; use the true laser line; fit or justify the Drude damping.
3. Validate the whole pipeline on a metal with published Raman intensities before applying it to YB4.
4. Polarized spectra from the experiment would let us compare irrep by irrep instead of peak positions only.

## 12. Mistakes we made on the way (so you recognise them)
- VASP exited at once: no POSCAR in the run directory (we had copied only a CONTCAR).
- `ZBRENT: fatal error in bracketing`: conjugate-gradient line search at the noise floor; restart with IBRION=1.
- `phonopy -d` rejected: phonopy 4 renamed setup commands to `phonopy-init`.
- Assumed the wrong space group (Pbam vs P4/mbm); the label check in the script caught it.
- Treated the 3.7e-4 eV/A structure as good enough for a first pass; the tight rerun exists to prove it.
