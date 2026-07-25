# Install

The base tool (Scanner + Calculator, CLI and web) needs only the core dependencies:

```bash
git clone https://github.com/VishnuTejasT/LOCKR.git
cd LOCKR
conda create -n igem python=3.10
conda activate igem
pip install -e .
```

## Grafting feature (optional)

The Scanner's "Graft into lucCage" step needs PyRosetta, a ~500MB structural
modeling package that isn't installed by default. Everything else in the
tool works without it, if it's missing, the graft button is simply disabled
with a note pointing back here.

```bash
conda activate igem
pip install pyrosetta \
  --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
# East coast mirror if West is slow:
# pip install pyrosetta \
#   --find-links https://graylab.jhu.edu/download/PyRosetta4/archive/release-quarterly/release

# Python 3.10 is required (already used by the igem conda env)
# Verify:
python3 -c "from pyrosetta import init; init(); print('OK')"
```

No registration or license token is needed for the public quarterly mirror.
Restart `lockr serve` after installing so the API picks it up.
