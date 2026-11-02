#!/usr/bin/env python3
"""
Build the frozen-phonon cells for Raman intensities: for each selected Gamma mode, structures
displaced by +Q and -Q along the mode, plus the undisplaced reference, each with an optics
(LOPTICS) INCAR. Run in ELFnet from 03_raman after the tight-structure phonons are analysed:

    python ../scripts/prepare_raman.py --select pilot      # the 3 strong A1g modes (>700 cm^-1)
    python ../scripts/prepare_raman.py --select all        # every Raman-active mode (incl. both Eg partners)

Mode coordinate: displacement dr_i = Q * e_i / sqrt(m_i), e = unit-norm mass-weighted eigenvector, so
Q has units A*sqrt(amu) and the Raman tensor is alpha = d(eps)/dQ. --delta sets the largest atomic
displacement (A); Q follows from it. Degenerate modes (Eg) are re-orthonormalised as real vectors, which
leaves the doublet's summed intensity unchanged.
"""
import argparse, csv, json, shutil
from pathlib import Path
import numpy as np
import phonopy

HERE = Path(__file__).resolve().parent
INCAR_OPTICS = """SYSTEM = YB4_raman_eps
ISTART = 0
ICHARG = 2
ENCUT = 910
PREC = Accurate
ALGO = Normal
ISMEAR = 1
SIGMA = 0.05
LREAL = .FALSE.
ADDGRID = .TRUE.
IBRION = -1
NSW = 0
EDIFF = 1E-8
LOPTICS = .TRUE.
NEDOS = 4001
CSHIFT = 0.05
NBANDS = {nbands}
LWAVE = .FALSE.
LCHARG = .FALSE.
"""
SBATCH = """#!/bin/bash
#SBATCH --job-name=YB4t_raman_eps
#SBATCH --ntasks=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --array=1-{n}%4
#SBATCH --output=logs/eps_%A_%a.log
#SBATCH --error=logs/eps_%A_%a.err
export OMP_NUM_THREADS=1
module purge
# load your site VASP modules here
cd "$SLURM_SUBMIT_DIR/$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" jobs.txt)"
srun --mpi=pmix --cpu-bind=none "${VASP_BIN:-vasp_std}"
"""


def write_poscar(path, cell, frac, symbols):
    order = []
    for s in symbols:
        if s not in order:
            order.append(s)
    cnt = [sum(1 for s in symbols if s == o) for o in order]
    idx = [i for o in order for i, s in enumerate(symbols) if s == o]
    L = ["Raman frozen-phonon cell", "1.0"]
    L += [" ".join(f"{x:.16f}" for x in v) for v in cell]
    L += [" ".join(order), " ".join(map(str, cnt)), "Direct"]
    L += [" ".join(f"{x:.16f}" for x in frac[i]) for i in idx]
    Path(path).write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", default="../02_analysis")
    ap.add_argument("--select", default="pilot", help="pilot | all | comma list of mode numbers")
    ap.add_argument("--delta", type=float, default=0.01, help="largest atomic displacement (A)")
    ap.add_argument("--nbands", type=int, default=160)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    A, out = Path(args.analysis), Path(args.outdir)

    rows = list(csv.DictReader(open(A / "gamma_modes.csv")))
    raman = [r for r in rows if r["activity"] == "Raman"]
    if args.select == "pilot":
        sel = [r for r in raman if r["irrep"] == "A1g" and float(r["freq_cm1"]) > 700]
    elif args.select == "all":
        sel = raman
    else:
        want = {int(x) for x in args.select.split(",")}
        sel = [r for r in rows if int(r["mode"]) in want]
    print(f"{len(sel)} modes selected:", [(r["mode"], r["freq_cm1"], r["irrep"]) for r in sel])

    ph = phonopy.load(str(A / "phonopy_disp.yaml"), force_sets_filename=str(A / "FORCE_SETS"),
                      is_nac=False, produce_fc=True, log_level=0)
    uc = ph.unitcell
    cell, masses, symbols = uc.cell, np.array(uc.masses), list(uc.symbols)
    frac0, nat = uc.scaled_positions, len(uc)
    ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
    qp = ph.get_qpoints_dict()
    freqs_cm = qp["frequencies"][0] * 33.35641
    ev = qp["eigenvectors"][0]

    # real, orthonormal eigenvectors (phase fixed per column, then QR inside degenerate groups)
    E = np.zeros_like(ev, dtype=float)
    worst_imag = 0.0
    for m in range(ev.shape[1]):
        v = ev[:, m]
        v = v * np.exp(-1j * np.angle(v[np.argmax(np.abs(v))]))
        worst_imag = max(worst_imag, float(np.abs(v.imag).max()))
        E[:, m] = v.real
    groups, used = [], set()
    for m in range(len(freqs_cm)):
        if m in used:
            continue
        g = [k for k in range(len(freqs_cm)) if abs(freqs_cm[k] - freqs_cm[m]) < 0.05]
        used.update(g); groups.append(g)
        if len(g) > 1:
            q, _ = np.linalg.qr(E[:, g]); E[:, g] = q
    print(f"max imaginary residual after phase fix: {worst_imag:.2e}")

    (out / "logs").mkdir(exist_ok=True)
    jobs, meta = [], []

    def make(name, frac):
        d = out / name; d.mkdir(exist_ok=True)
        write_poscar(d / "POSCAR", cell, frac, symbols)
        (d / "INCAR").write_text(INCAR_OPTICS.format(nbands=args.nbands))
        for f in ("KPOINTS", "POTCAR"):
            shutil.copy(A.parent / "01_displacements" / f, d / f)
        jobs.append(name)

    make("ref", frac0)
    for r in sel:
        m = int(r["mode"]) - 1
        u = E[:, m].reshape(nat, 3) / np.sqrt(masses)[:, None]       # A per unit Q
        s = args.delta / np.linalg.norm(u, axis=1).max()             # Q (A sqrt(amu))
        for sign, tag in ((+1, "p"), (-1, "m")):
            dr = sign * s * u
            make(f"m{r['mode']:0>3}_{tag}", frac0 + dr @ np.linalg.inv(cell))
        meta.append({"mode": int(r["mode"]), "freq_cm1": float(r["freq_cm1"]), "irrep": r["irrep"],
                     "Q_A_sqrt_amu": float(s), "max_disp_A": args.delta})
    (out / "raman_meta.json").write_text(json.dumps(
        {"modes": meta, "volume_A3": float(abs(np.linalg.det(cell)))}, indent=2))
    (out / "jobs.txt").write_text("\n".join(jobs) + "\n")
    (out / "eps_array.sbatch").write_text(SBATCH.format(n=len(jobs)))
    print(f"{len(jobs)} VASP runs listed in jobs.txt; submit with: sbatch eps_array.sbatch")


if __name__ == "__main__":
    main()
