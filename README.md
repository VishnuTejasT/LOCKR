# LOCKR Thermodynamics Design Tool


This software is a CPU-only thermodynamic engine for LOCKR biosensors. It essentially works by scanning a binding region in the intended latch, and identifies liabilities and problems that can restrict cage-key thermodynamic interactioons. It provides a fix, and outputs a final K_CK (Cage-Key) score. Then, the user can use the K_CK score from the scanner in the Calculator section to calculate fold-change and other variables to determine whether their LOCKR design is a good fit. Both sections utlize the same physics, so they share the same thermodynamic engine.

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

## Using the tool

Activate the environment first:

```bash
conda activate igem
lockr --help          # top-level help
```

### Option 1, web app (recommended)

```bash
lockr serve
```

It should provide a local link for you to paste into your browser to get started!

### Option 2, command line

Scan a binder for cage-key liabilities:

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

### More help

```bash
lockr scan --help
lockr fc --help
lockr serve --help
```
