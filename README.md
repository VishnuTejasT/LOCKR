# LOCKR Thermodynamics Design Tool


This software is a CPU-only thermodynamic engine for LOCKR biosensors. It essentially works by scanning a binding region in the intended latch, and identifies liabilities and problems that can restrict cage-key thermodynamic interactioons. Before it gets to charge liabilities, it first checks whether the binder is even shaped right, since the lucCage latch is a helix, so a binder has to be helical to thread into it. It provides a fix, and outputs a final K_CK (Cage-Key) score. Then, the user can use the K_CK score from the scanner in the Calculator section to calculate fold-change and other variables to determine whether their LOCKR design is a good fit. Both sections utlize the same physics, so they share the same thermodynamic engine.

## Prerequisites

To use this software, please install the following:

```bash
1. Conda/Miniconda     # Visit this site on information to download Conda: https://www.anaconda.com/docs/getting-started/installation
2. Python 3.10         # Visit this site on information to download Python 3.10: https://www.python.org/downloads/
```

## Install

Copy and paste these commands into your Terminal/Command-Prompt

```bash
git clone https://github.com/VishnuTejasT/LOCKR.git
cd LOCKR
conda create -n igem python=3.10
conda activate igem
pip install -e .
```

## Grafting feature (optional)

The Scanner's "Graft into lucCage" step needs PyRosetta, and it is not installed by default because it's a ~500mb download.. Everything else in the software works without it, and if it is not downloaded, then the graft button is disabled.

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

## Using the tool

Activate the environment first:

```bash
conda activate igem
lockr --help          # top-level help
```

### Option 1, web app (recommended)

```bash
lockr serve
lockr serve --port 8001    # if 8000 is already taken by something else
```

It should provide a local link for you to paste into your browser to get started! The web app has two tabs: Scanner (structure pre-check, charge scanning, binder suggestions, and grafting) and Calculator (fold-change predictions).

Before running the charge scan, the Scanner checks whether the binder is helical enough to be grafted into LucCage's latch, using per-residue helix propensity, capping, and amphipathicity to flag sequences that look cyclic instead of linear. It's a propensity estimate from sequence alone, not a structural prediction, so a low score is merely a warning.

Grafting needs PyRosetta, which isn't installed by default since it's a ~500MB download. Everything else works fine without it, the graft button just stays disabled until it's there. See `INSTALL.md` for the one-line install if you want that feature.

### Option 2, command line

Scan a binder for cage-key liabilities. This also runs the same helical structure pre-check as the web app first, so the output leads with a Structure line (helix confidence, band, and any shape issues) before the charge liability numbers:

```bash
lockr scan SEQUENCE
lockr scan SEQUENCE --ph 7.4                     # pH of the intended solution
lockr scan SEQUENCE --window 1:17                # only scan this residue window
lockr scan SEQUENCE --preserve 1,2,11,12,15      # never mutate these binding residues
lockr scan SEQUENCE --policy conservative        # D->N, E->Q (keeps shape and H-bonding)
lockr scan SEQUENCE --policy neutralizing        # D->A, E->A
lockr scan SEQUENCE --no-suggest                 # skip the variant suggestion
lockr scan SEQUENCE --json                       # machine-readable output
lockr scan --file sequences.fasta                # FASTA, raw, or mixed input
lockr scan --file sequences.fasta --json
```

Predict fold-change:

```bash
lockr fc --k-ck 10 --k-open 0.001 --pull 10 --luckey 500
lockr fc --k-ck 10 --k-open 0.001 --pull 10 --luckey 500 --k-target 50 --target 5
lockr fc --k-ck 10 --k-open 0.001 --pull 10 --luckey 500 --json
```

Leave `--k-target`/`--target` off to assume the target is fully saturating. All values are in nM except from `--k-open` and `--pull`, which are dimensionless.

Graft a binder into the lucCage latch (needs PyRosetta, see above):

```bash
lockr graft SEQUENCE                             # scans every latch position and keeps the best
lockr graft SEQUENCE --position 327              # graft at one specific position instead of scanning all possible residues
lockr graft SEQUENCE --latch-start 325 --latch-end 359  # override the latch window on the template
lockr graft SEQUENCE --out grafted.pdb           # saves the grafted structure to a specific path
lockr graft SEQUENCE --json
lockr graft --status                             # check whether PyRosetta and the template are even available
```

### More help

```bash
lockr scan --help
lockr fc --help
lockr serve --help
lockr graft --help
```
