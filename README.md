# LOCKR Thermodynamics Design Tool

This software is a CPU-only thermodynamic engine for LOCKR biosensors. It essentially works by scanning a binding region in the intended latch, and identifies liabilities and problems that can restrict cage-key thermodynamic interactioons. Before it gets to charge liabilities, it first checks whether the binder is even shaped right, since the lucCage latch is a helix, so a binder has to be helical to thread into it. It provides a fix, and outputs a final K_CK (Cage-Key) score. Then, the user can use the K_CK score from the scanner in the Calculator section to calculate fold-change and other variables to determine whether their LOCKR design is a good fit. Both sections utlize the same physics, so they share the same thermodynamic engine.

## Install

The core tool needs Python 3.10 or newer. Conda is a convenient way to keep its dependencies separate.

```bash
git clone https://github.com/VishnuTejasT/LOCKR.git
cd LOCKR
conda create -n lockr python=3.10
conda activate lockr
pip install -e .
```

To install the optional test dependencies:

```bash
pip install -e ".[test]"
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

## Start the web app

Activate the environment and run:

```bash
conda activate lockr
lockr serve
```

Open the local address printed in the terminal (usually `http://127.0.0.1:8000`). If port 8000 is in use, choose another:

```bash
lockr serve --port 8001
```

The interactive API reference is at `http://127.0.0.1:8000/docs` while the server is running. Stop the server with Ctrl+C.

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

## Scanner

Open the **Scanner** tabs and select a single-sequence option or a batch mode. The sequences can only be inputted with the standard 20 amino acids.

The original and charge-optimzied binders can be grafted into the LOCKR system via the PyRosetta module, which is ian optional dependancy.

### Scan Procedure

Paste the binder sequence and select the solution's pH (defaulted to 7.4). The scan can be run and the following results are included:

— **Structural Check:** This checks whether the sequence is a helix or not by estimating helix propensity, hydrophobic moments, and potential stabilizing salt bridges. Lower score -> less confidence in grafting compatibility

— **Liability Score and Band:** Summarizes acidic residues counted in a selected window of the binder. A chart shows contirbutions per resdiues and which residues specifically drive liability. Acidic residues in a preserved positions are not counted for the score.

— **Net charge:** Calculated with the slected pH value

— **Estimated K_CK** A calibrated estimate that was derived from the libaility model. It's calibrated using the ECLIPSE data set.

— **Suggested Variant:** Changes the flagged amino acid residues while leaving the preserved psoitions alone. There are 2 methods for substitution: Conservative(D->N and E->Q) and Neutralizing(D->A and D->E)



## Calculator & Dosing

### Fold-change calculator

The **Calculator & Dosing** tab starts with example ECLIPSE values. You can choose a preset or enter your own:


The calculator reports fold-change, a best-case value for the selected pull, key-to-cage dominance ratio, limiting regime, recommendations, and warnings. It also plots how changing LucKey over time affects fold-change in the form of a graph.

To model a finite target concentration, input both K_target (nM) and target concentration (nM). K_target can be entered directly, or estimated from a docking free-energy score. Without a reference anchor, that conversion is a rough directional estimate. The optional clinical threshold adds a comparison against a user-specified fold-change cutoff.


### Reporter-key dosing

Expand **Reporter-key dosing** below the calculator. Enter K_CK, total cage concentration, reporter-key stock concentration, and assay volume. Choose either the desired fraction of cage bound or a desired brightness/fold-change jump. The tool calculates the key concentration and stock volume to add, and plots the titration curve. Check the warnings and the volume against your pipetting range before preparing an assay.

**Reporter Key Dosing** can be expanded as a section in the Clauclator tab and you enter specific parameters like K_CK, Cage cocnentration, key stock concentration, and assay volume. These values, alongside a desired fold-change, will work together to calculate for a key concentration with a stock volume to add, while also plotting it in a titration curve.

## Design Map

Presets cover cage tightness versus pull, key concentration and affinity, temperature response, and minimum detectable target. Set the preset and input the desired values to plot the design map. The design map is utilzied to compare the effects of 2 parameters on fold change like k_open and pull strength.


The **Mark suggested designs** control runs inverse design for a target fold-change and overlays combinations that meet it when the map axes are cage looseness and key concentration. The API also exposes inverse solving for other supported parameter combinations.