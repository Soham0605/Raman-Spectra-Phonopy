#!/bin/bash
# Run AFTER 00_relax finishes and its max force is < 1e-4 eV/A.
# Steps 1-2: generate displaced cells and a SLURM array that runs one static VASP job each.
# Step 3 (after the array finishes): bash step1_3_displacements.sh collect
set -e
source ~/miniconda3/etc/profile.d/conda.sh; conda activate ELFnet
D=${D:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}   # workflow root; override with D=...
cd $D/01_displacements
if [ "$1" = collect ]; then
  phonopy-init -f disp-*/vasprun.xml        # fails loudly on unfinished runs
  cp FORCE_SETS ../02_analysis/; cp phonopy_disp.yaml POSCAR ../02_analysis/
  exit 0
fi
cp ../${RELAX_DIR:-00_relax}/CONTCAR POSCAR   # set RELAX_DIR=00_relax_2 to use the tighter restart
phonopy-init -d --dim="1 1 1" -c POSCAR     # writes POSCAR-001..., phonopy_disp.yaml
N=$(ls POSCAR-[0-9]* | wc -l); echo "$N displaced cells"
for f in POSCAR-[0-9]*; do
  id=${f#POSCAR-}; mkdir -p disp-$id/logs
  cp $f disp-$id/POSCAR; cp INCAR.static disp-$id/INCAR; cp KPOINTS POTCAR disp-$id/
done
sed "s/__N__/$N/" ../scripts/disp_array.sbatch > disp_array.sbatch
echo "submit with: cd $D/01_displacements && sbatch disp_array.sbatch"
