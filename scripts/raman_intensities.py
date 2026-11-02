#!/usr/bin/env python3
"""
Nonresonant first-order Raman intensities from frozen-phonon dielectric functions.

For each mode, alpha_ij = [eps_ij(+Q) - eps_ij(-Q)] / (2Q) evaluated at the laser energy, where eps is
VASP's frequency-dependent dielectric function (LOPTICS, interband) and Q the mode coordinate in
A*sqrt(amu) (see prepare_raman.py). For a metal the free-carrier term is added when the plasma frequency
can be read from the OUTCAR: eps_D = -wp^2/(w^2 + i*gamma*w).  Intensity of a mode, powder average:
    a = tr(alpha)/3,  g^2 = 1/2[|axx-ayy|^2+|ayy-azz|^2+|azz-axx|^2 + 6(|axy|^2+|ayz|^2+|azx|^2)]
    I_par = 45|a|^2 + 4 g^2,  I_perp = 3 g^2,  I_total = 45|a|^2 + 7 g^2
    I = (wL - w)^4 / w * (n(w,T) + 1) * I_total       (arbitrary units, relative across modes)
Degenerate partners (same frequency) are summed into one peak.

    python ../scripts/raman_intensities.py --laser-nm 532 --temp 300 --drude-gamma 0.1
"""
import argparse, json, re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

KB_CM = 0.6950348  # cm^-1 per K
EV_PER_NM = 1239.84198  # eV*nm


def read_eps(vasprun, e_laser):
    """complex 3x3 eps at e_laser (eV) from vasprun.xml <dielectricfunction>."""
    df = ET.parse(vasprun).getroot().find(".//dielectricfunction")
    if df is None:
        raise RuntimeError(f"no <dielectricfunction> in {vasprun} (was LOPTICS on?)")
    def block(tag):
        rows = [[float(x) for x in r.text.split()] for r in df.find(tag).iter("r")]
        return np.array(rows)
    im, re_ = block("imag"), block("real")
    e = im[:, 0]
    comp = lambda arr: [np.interp(e_laser, e, arr[:, k]) for k in range(1, 7)]
    r6, i6 = comp(re_), comp(im)
    eps6 = np.array(r6) + 1j * np.array(i6)           # xx yy zz xy yz zx
    xx, yy, zz, xy, yz, zx = eps6
    return np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])


def read_plasma(outcar):
    """3x3 plasma frequency squared (eV^2) from OUTCAR, or None if not found."""
    p = Path(outcar)
    if not p.exists():
        return None
    t = p.read_text(errors="replace")
    m = re.search(r"plasma frequency squared[^\n]*\n((?:\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*\n){3})", t, re.I)
    if not m:
        return None
    return np.array([[float(x) for x in l.split()] for l in m.group(1).strip().split("\n")])


def eps_total(d, e_laser, gamma, use_drude):
    eps = read_eps(d / "vasprun.xml", e_laser)
    drude = None
    if use_drude:
        wp2 = read_plasma(d / "OUTCAR")
        if wp2 is not None:
            drude = -wp2 / (e_laser**2 + 1j * gamma * e_laser)
            eps = eps + drude
    return eps, drude


def invariants(a):
    tr = np.trace(a) / 3
    g2 = 0.5 * (abs(a[0, 0] - a[1, 1])**2 + abs(a[1, 1] - a[2, 2])**2 + abs(a[2, 2] - a[0, 0])**2
                + 6 * (abs(a[0, 1])**2 + abs(a[1, 2])**2 + abs(a[2, 0])**2))
    return abs(tr)**2, g2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".")
    ap.add_argument("--laser-nm", type=float, default=532.0)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--drude-gamma", type=float, default=0.1, help="Drude damping (eV); assumption")
    ap.add_argument("--no-drude", action="store_true")
    args = ap.parse_args()
    D = Path(args.dir)
    meta = json.loads((D / "raman_meta.json").read_text())
    e_l = EV_PER_NM / args.laser_nm
    wl = 1e7 / args.laser_nm
    print(f"laser {args.laser_nm} nm = {e_l:.3f} eV; T = {args.temp} K; Drude gamma {args.drude_gamma} eV")

    eps0, _ = eps_total(D / "ref", e_l, args.drude_gamma, not args.no_drude)
    rows = []
    for m in meta["modes"]:
        n = f"m{m['mode']:03d}"
        ep, dp = eps_total(D / f"{n}_p", e_l, args.drude_gamma, not args.no_drude)
        em, _ = eps_total(D / f"{n}_m", e_l, args.drude_gamma, not args.no_drude)
        Q = m["Q_A_sqrt_amu"]
        alpha = (ep - em) / (2 * Q)
        a2, g2 = invariants(alpha)
        nonlin = np.abs(ep + em - 2 * eps0).max() / max(np.abs(ep - em).max(), 1e-30)
        w = m["freq_cm1"]
        bose = 1.0 / np.expm1(w / (KB_CM * args.temp)) + 1.0
        pref = (wl - w)**4 / w * bose
        rows.append({**m, "a2": a2, "g2": g2, "I_total": pref * (45 * a2 + 7 * g2),
                     "I_par": pref * (45 * a2 + 4 * g2), "I_perp": pref * 3 * g2,
                     "nonlinearity": float(nonlin), "drude_used": dp is not None})
    # sum degenerate partners into one peak
    peaks = {}
    for r in rows:
        k = (r["irrep"], round(r["freq_cm1"], 1))
        peaks.setdefault(k, {"irrep": r["irrep"], "freq_cm1": r["freq_cm1"], "I_total": 0.0, "I_par": 0.0,
                             "I_perp": 0.0, "nonlinearity": 0.0, "drude_used": r["drude_used"]})
        for f in ("I_total", "I_par", "I_perp"):
            peaks[k][f] += r[f]
        peaks[k]["nonlinearity"] = max(peaks[k]["nonlinearity"], r["nonlinearity"])
    out = sorted(peaks.values(), key=lambda p: p["freq_cm1"])
    mx = max(p["I_total"] for p in out)
    lines = ["freq_cm1,irrep,I_rel,I_par_rel,I_perp_rel,nonlinearity,drude_used"]
    for p in out:
        lines.append(f"{p['freq_cm1']:.2f},{p['irrep']},{p['I_total']/mx:.4f},{p['I_par']/mx:.4f},"
                     f"{p['I_perp']/mx:.4f},{p['nonlinearity']:.3f},{p['drude_used']}")
    (D / "raman_intensities.csv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("nonlinearity = |eps(+Q)+eps(-Q)-2 eps(0)| / |eps(+Q)-eps(-Q)|; should be << 1 for a valid finite difference")
    if not any(p["drude_used"] for p in out):
        print("WARNING: Drude (free-carrier) term NOT included: plasma frequency not found in OUTCAR.")


if __name__ == "__main__":
    main()
