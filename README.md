# LOCKR Thermodynamics Design Tool

A CPU-only thermodynamic toolkit for designing and troubleshooting LOCKR
biosensors. It has three parts that share one physics engine:

1. **Scanner** — paste a binder sequence; it flags the charged residues that
   weaken cage–key binding, proposes a cleaned-up variant, and reports an
   estimated **K_CK** (cage–key affinity).
2. **Calculator** — take that K_CK (or your own numbers) and predict the
   sensor's **fold-change**: how much brighter it gets when the target appears.
   It tells you, in plain English, whether that number is good and what single
   change is most likely to improve it.
3. **Assembly Check** — verify that a full LOCKR sequence is put together
   correctly (binder, spacer, protected region, and latch all in the right
   place).

Both the Scanner and Calculator use the same engine, so a K_CK from the Scanner
drops straight into the Calculator.


## Install

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

### Option 1 — web app (recommended)

```bash
lockr serve
```

It should provide a local link for you to paste into your browser to get started!

### Option 2 — command line

Scan a binder for cage–key liabilities:

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

Leave `--k-target`/`--target` off to assume the target is fully saturating. All values are in nM except from `--k-open` and `--pull`.

### More help

```bash
lockr scan --help
lockr fc --help
lockr serve --help
```
