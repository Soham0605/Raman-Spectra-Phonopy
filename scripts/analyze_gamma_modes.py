#!/usr/bin/env python3
"""
Steps 4-7 of the YB4 Gamma-point phonon workflow, using the phonopy Python API.

Run in the ELFnet conda env, from 02_analysis (needs phonopy_disp.yaml, FORCE_SETS):

    conda activate ELFnet
    python ../scripts/analyze_gamma_modes.py            # defaults
    python ../scripts/analyze_gamma_modes.py --amp 0.4 --frames 20

Writes into --outdir (default: this directory):
    gamma_modes.csv / gamma_modes.txt   frequency, irrep, activity, Y%, B%, dominant Wyckoff site
    irreps.yaml                         phonopy's own irrep output (written in the current directory)
    animations/mode_XXX_<freq>cm_<irrep>.extxyz   one per Raman-active mode (Step 7)

Activity classification assumes point group D4h (YB4 is P4/mbm, #127; the label set is checked):
    A1g, B1g, B2g, Eg -> Raman;  A2u, Eu -> IR;  A1u, A2g, B1u, B2u -> silent.
    The three acoustic modes (omega ~ 0, A2u + Eu) are flagged "acoustic".
If phonopy returns labels outside D4h the script stops rather than guess.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
import phonopy

THZ_TO_CM = 33.35641  # 1 THz = 33.35641 cm^-1
# D4h (P4/mbm, YB4). Raman: A1g, B1g, B2g, Eg. IR: A2u, Eu. Silent: A1u, A2g, B1u, B2u.
RAMAN = {"A1g", "B1g", "B2g", "Eg"}
IR = {"A2u", "Eu"}
SILENT = {"A1u", "A2g", "B1u", "B2u"}


def classify(label, freq_cm, acoustic_cut=5.0):
    if label is None:
        # phonopy cannot label the 3 acoustic modes (A2u + Eu); everything else unlabelled is a problem
        return "acoustic" if abs(freq_cm) < acoustic_cut else "unlabelled"
    base = label.replace("'", "").replace("′", "")
    if base in RAMAN:
        return "Raman"
    if base in IR:
        return "acoustic" if abs(freq_cm) < acoustic_cut else "IR"
    if base in SILENT:
        return "silent"
    return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disp-yaml", default="phonopy_disp.yaml")
    ap.add_argument("--force-sets", default="FORCE_SETS")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--amp", type=float, default=0.4, help="max atomic displacement in animation (A)")
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--tol", type=float, default=1e-4, help="degeneracy tolerance (THz)")
    args = ap.parse_args()

    out = Path(args.outdir)
    (out / "animations").mkdir(parents=True, exist_ok=True)

    # Steps 4: Gamma frequencies and eigenvectors (1x1x1 supercell, no NAC for a metal)
    ph = phonopy.load(args.disp_yaml, force_sets_filename=args.force_sets,
                      is_nac=False, produce_fc=True, log_level=0)
    uc = ph.unitcell
    n_atoms = len(uc)
    ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
    qp = ph.get_qpoints_dict()
    freqs_thz = qp["frequencies"][0]
    eigvecs = qp["eigenvectors"][0]          # (3N, 3N), columns are modes, mass-weighted, norm 1
    freqs_cm = freqs_thz * THZ_TO_CM
    print(f"{len(freqs_cm)} modes; lowest 6 (cm^-1): {np.round(freqs_cm[:6], 3)}")
    n_imag = int((freqs_cm < -1.0).sum())
    print(f"imaginary modes (< -1 cm^-1): {n_imag}")

    # Step 5: irreps
    ph.set_irreps([0, 0, 0], degeneracy_tolerance=args.tol)
    ir = ph.irreps
    ph.write_yaml_irreps(show_irreps=False)  # phonopy writes irreps.yaml in the current directory
    labels_by_set = getattr(ir, "_ir_labels", None)
    mode_label = [None] * len(freqs_cm)
    if labels_by_set is not None:
        for lab, idx in zip(labels_by_set, ir.band_indices):
            for i in idx:
                mode_label[i] = lab
    pg = getattr(ir, "_pointgroup_symbol", None) or getattr(ir, "pointgroup_symbol", None)
    if labels_by_set is None or all(l is None for l in mode_label):
        raise SystemExit("phonopy returned no irrep labels (point group not recognised); "
                         "see gamma_irreps.yaml characters and label by hand.")
    unknown = {l for l in mode_label if l and l.replace("'", "") not in RAMAN | IR | SILENT}
    if unknown:
        raise SystemExit(f"Labels {unknown} are not D4h; extend classify() before trusting the flags.")

    # Step 6: participation. Eigenvector columns are already sqrt(m)-weighted and normalised to 1.
    symbols = uc.symbols
    masses = np.array(uc.masses)
    from phonopy.structure.symmetry import Symmetry
    sym = Symmetry(uc, symprec=1e-3)
    ds = sym.dataset
    wyck = ds.wyckoffs if hasattr(ds, "wyckoffs") else ds["wyckoffs"]
    equiv = ds.equivalent_atoms if hasattr(ds, "equivalent_atoms") else ds["equivalent_atoms"]
    sg = ds.international if hasattr(ds, "international") else ds["international"]
    print(f"space group {sg}")
    # site label like "4g" = multiplicity + Wyckoff letter, tagged with the species
    mult = {e: int((equiv == e).sum()) for e in set(equiv)}
    site_of = [f"{mult[equiv[i]]}{wyck[i]}({symbols[i]})" for i in range(n_atoms)]
    sites = sorted(set(site_of))

    rows = []
    for m in range(len(freqs_cm)):
        e = eigvecs[:, m].reshape(n_atoms, 3)
        w = (np.abs(e) ** 2).sum(axis=1)         # per-atom participation, sums to 1
        w /= w.sum()
        y = float(sum(w[i] for i in range(n_atoms) if symbols[i] == "Y"))
        per_site = {s: float(sum(w[i] for i in range(n_atoms) if site_of[i] == s)) for s in sites}
        dom = max(per_site, key=per_site.get)
        rows.append({
            "mode": m + 1, "freq_cm1": round(float(freqs_cm[m]), 2), "freq_THz": round(float(freqs_thz[m]), 4),
            "irrep": mode_label[m], "activity": classify(mode_label[m], freqs_cm[m]),
            "Y_pct": round(100 * y, 1), "B_pct": round(100 * (1 - y), 1),
            "dominant_site": f"{dom} {100 * per_site[dom]:.0f}%",
        })

    with open(out / "gamma_modes.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    hdr = f"{'#':>3} {'cm-1':>9} {'irrep':>6} {'activity':>9} {'Y%':>6} {'B%':>6}  dominant site"
    lines = [hdr] + [f"{r['mode']:>3} {r['freq_cm1']:>9.2f} {str(r['irrep']):>6} {r['activity']:>9} "
                     f"{r['Y_pct']:>6.1f} {r['B_pct']:>6.1f}  {r['dominant_site']}" for r in rows]
    counts = {}
    for r in rows:
        counts[r["irrep"]] = counts.get(r["irrep"], 0) + 1
    lines += ["", "irrep counts (acoustic included): " + ", ".join(f"{k}:{v}" for k, v in sorted(counts.items(), key=lambda kv: str(kv[0]))),
              f"Raman-active modes (counting degenerate partners): {sum(r['activity'] == 'Raman' for r in rows)}",
              f"Raman-active distinct frequencies (peaks): {len({(r['irrep'], r['freq_cm1']) for r in rows if r['activity'] == 'Raman'})}",
              f"imaginary modes: {n_imag}"]
    (out / "gamma_modes.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    # Step 7: extxyz for each Raman-active mode. Displacement u_i ~ e_i / sqrt(m_i), scaled so the
    # largest atomic amplitude equals --amp; frame k = r + amp * cos(2 pi k/N) * u_hat.
    lat = uc.cell
    pos0 = uc.positions
    lat_str = " ".join(f"{x:.8f}" for x in lat.flatten())
    for r in rows:
        if r["activity"] != "Raman":
            continue
        m = r["mode"] - 1
        e = eigvecs[:, m].reshape(n_atoms, 3)
        u = e / np.sqrt(masses)[:, None]
        # eigenvectors are complex; at Gamma use the real part after fixing the global phase
        k = np.argmax(np.abs(u).ravel())
        u = (u * np.exp(-1j * np.angle(u.ravel()[k]))).real
        u *= 1.0 / np.linalg.norm(u, axis=1).max()
        name = f"mode_{r['mode']:03d}_{r['freq_cm1']:.0f}cm_{r['irrep']}.extxyz"
        with open(out / "animations" / name, "w") as f:
            for fr in range(args.frames):
                d = args.amp * np.cos(2 * np.pi * fr / args.frames) * u
                f.write(f"{n_atoms}\n")
                f.write(f'Lattice="{lat_str}" Properties=species:S:1:pos:R:3:disp:R:3 pbc="T T T" '
                        f'frame={fr} freq_cm1={r["freq_cm1"]} irrep={r["irrep"]}\n')
                for i in range(n_atoms):
                    p = pos0[i] + d[i]
                    f.write(f"{symbols[i]} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} "
                            f"{d[i][0]:.6f} {d[i][1]:.6f} {d[i][2]:.6f}\n")
    print(f"wrote animations to {out / 'animations'}")


if __name__ == "__main__":
    main()
