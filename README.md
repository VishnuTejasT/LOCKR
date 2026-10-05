# LOCKR design tool

LOCKR is a local tool for exploring thermodynamic designs of LOCKR biosensors. It combines a peptide sequence scanner with fold-change, dosing, parameter-sweep, and design-map calculations. The calculations run on your machine; the web interface is served locally by the Python API.

The sequence scanner estimates charge-related cage–key (K_CK) liabilities and checks whether a binder sequence looks compatible with a helical latch. These are screening estimates, not structural predictions or substitutes for measurements. The thermodynamic tools use the values you provide, or the ECLIPSE starting values shown in the interface.

## What you can do

- **Scanner:** inspect one binder or a batch, review acidic-residue contributions, get a charge-optimized suggestion, and optionally thread a binder into the lucCage latch.
- **Calculator & Dosing:** estimate fold-change, see which parameters limit the signal, estimate reporter-key volume to add, and explore latch tuning.
- **Design Map:** plot how two parameters affect a selected metric, apply feasibility rules, and mark parameter combinations that meet a target fold-change.
- **API:** use the same calculations from scripts through the local FastAPI service. Assembly verification, inverse design, mutation maps, and supporting calculations are available through the API even when they do not have a separate web form.

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

### Optional: enable structure grafting

Grafting uses PyRosetta, which is a large optional dependency. The scanner, calculator, dosing tools, and design map work without it. Install PyRosetta in the same environment if you need grafting:

```bash
conda activate lockr
pip install pyrosetta \
  --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
```

An alternate mirror is available if the first is slow:

```bash
pip install pyrosetta \
  --find-links https://graylab.jhu.edu/download/PyRosetta4/archive/release-quarterly/release
```

Check the install with `python -c "from pyrosetta import init; init(); print('OK')"`, then restart the server. See [INSTALL.md](INSTALL.md) for the same setup notes.

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

## Scanner

Open the **Scanner** tab and choose single-sequence or batch mode. Sequences must use the 20 standard amino-acid letters. The scanner normalizes lowercase input to uppercase; it rejects gaps, stop symbols, and non-standard residues.

### Single sequence

Paste a binder sequence, choose the solution pH (default 7.4), then run **Scan**. Results include:

- **Structure check:** a sequence-based helix-confidence estimate, helix propensity, hydrophobic moment, potential stabilizing salt bridges, and any flagged issues. It is a propensity check, not a prediction of the folded structure. A low score is a warning about graft compatibility.
- **Liability score and band:** summarizes acidic residues counted in the selected sensitive window. Per-position contributions show which residues drive the score. Acidic residues outside the window can still appear in the report but do not add to this score.
- **Net charge:** calculated at the selected pH.
- **Estimated K_CK:** a calibrated estimate derived from the liability model. Calibration uses a small ECLIPSE data set; treat values far outside the usual 1–100 nM range as directional warnings, not precise affinities.
- **Suggested variant:** changes flagged acidic residues while leaving preserved positions alone. Conservative substitutions use D→N and E→Q; neutralizing substitutions use D→A and E→A. Review the sequence and structural warnings before using a suggestion.

The sensitive window and preserved positions are numbered from 1 at the start of the pasted binder. Use the window for the part of the binder whose charge should count toward the cage–key estimate. Preserve target-contact residues if the suggested substitutions must not change them.

### Batch scan

Switch to **Batch**, then paste FASTA, one sequence per line, or a mixture of the two. You can also upload a `.fa`, `.fasta`, or `.txt` file. FASTA headers become result IDs. Select a row to inspect its full scan, or export the summary table as CSV. The best candidate strip ranks sequences by liability score; it is a ranking aid, not a claim that the top sequence will bind or function experimentally.

### Graft into lucCage

After scanning, the graft panel can thread either the original sequence or the suggested variant into the bundled lucCage latch template. By default it scans latch positions 325–359 and reports the best-scoring position and Rosetta score. Enable the specific-position option to try one position. The output includes a grafted sequence and a PDB structure file to download.

Grafting is unavailable until PyRosetta is installed. The Rosetta score is a structural-fit heuristic; inspect the resulting structure and validate the design experimentally.

## Calculator & Dosing

### Fold-change calculator

The **Calculator & Dosing** tab starts with example ECLIPSE values. You can choose a preset or enter your own:

| Input | Meaning | Units |
|---|---|---|
| K_CK | Dissociation constant for reporter key binding to open cage; lower means tighter binding | nM |
| K_open | Basal opening equilibrium in the absence of target | Dimensionless |
| Pull | Allosteric coupling between target binding and cage opening | Dimensionless |
| lucKey | Reporter-key concentration in the assay | nM |

The calculator reports fold-change, a best-case value for the selected pull, key-to-cage dominance ratio, limiting regime, recommendations, and warnings. It also plots how changing lucKey or K_open affects the result. These plots vary one quantity at a time while holding the others at their current values.

For a saturating-target assumption, leave the target-affinity and target-concentration fields blank. To model a finite target concentration, supply both K_target (nM) and target concentration (nM). You may enter K_target directly, or estimate it from a docking free-energy score using a reference score and measured affinity. Without a reference anchor, that conversion is a rough directional estimate. The optional clinical threshold adds a comparison against a user-specified fold-change cutoff.

The **How do I get K_open?** helper offers three routes: use the base lucCage estimate (0.001), estimate a change from latch mutations and a destabilization assumption, or choose a value from the published-variant table. These are estimates; use an experimental measurement when available.

The **Analyze K_open Tuning** panel compares hypothetical latch destabilization scenarios against your current settings. You can optionally paste a latch sequence to identify candidate I/V positions and mark positions to preserve. Choose a conservative, moderate, or optimistic destabilization per mutation. This analysis does not design or validate mutations; it estimates how changing K_open could affect the model.

### Reporter-key dosing

Expand **Reporter-key dosing** below the calculator. Enter K_CK, total cage concentration, reporter-key stock concentration, and assay volume. Choose either the desired fraction of cage bound or a desired brightness/fold-change jump. The tool calculates the key concentration and stock volume to add, and plots the titration curve. Check the warnings and the volume against your pipetting range before preparing an assay.

## Design Map

Open **Design Map** to compare two parameters across a grid. Start with a preset question or choose **Custom**. Presets cover cage tightness versus pull, key concentration and affinity, temperature response, and minimum detectable target. Set the remaining design values under **Your design**, choose the metric, and draw the map.

In custom mode, choose x and y parameters, their ranges, number of steps, and linear or logarithmic spacing. Log spacing is useful when values span orders of magnitude. The map can show contours and can apply feasibility constraints such as a minimum fold-change. Hover or pin a point to read its values. Export the plot as PNG or the underlying grid as CSV.

The temperature preset needs heat-response data to make a temperature-dependent prediction. Under **Advanced**, supply the relevant enthalpy (dH), heat-capacity change (dCp), and, where available, melting temperature (Tm) values for the latch/cage, key, or target interaction. If the needed data are absent, the interface warns that the calculation remains isothermal; changing the temperature slider alone does not create a temperature model. Temperature values are in kelvin in the API and are also shown in Celsius in the interface.

The **Mark suggested designs** control runs inverse design for a target fold-change and overlays combinations that meet it when the map axes are cage looseness and key concentration. The API also exposes inverse solving for other supported parameter combinations.

## Command line

The CLI covers serving the app, scanning sequences, fold-change calculation, and grafting. Run `lockr --help` or a subcommand's `--help` for the full option list.

### Scan a sequence

```bash
lockr scan LISDAELEAIFAEELDC
lockr scan LISDAELEAIFAEELDC --ph 7.4
lockr scan LISDAELEAIFAEELDC --window 1:17 --preserve 1,2,11,12,15
lockr scan LISDAELEAIFAEELDC --policy neutralizing
lockr scan LISDAELEAIFAEELDC --no-suggest
lockr scan --file sequences.fasta
lockr scan --file sequences.fasta --json
```

`--window START:END` and `--preserve` use 1-based binder positions. The default policy is `conservative`. Input files may contain FASTA records or plain sequences, one per line. `--json` emits a `results` array for scripts instead of the readable summary.

### Calculate fold-change

```bash
lockr fc --k-ck 10 --k-open 0.001 --pull 10 --luckey 500
lockr fc --k-ck 10 --k-open 0.001 --pull 10 --luckey 500 --k-target 50 --target 5
lockr fc --k-ck 10 --k-open 0.001 --pull 10 --luckey 500 --json
```

K_CK, lucKey, K_target, and target concentration are in nM. K_open and pull are dimensionless. `--k-target` and `--target` must be supplied together; omit both for saturating target. The CLI accepts the same thermodynamic core calculation as the web calculator, but the temperature controls are available through the API and Design Map.

### Graft a binder

```bash
lockr graft LISDAELEAIFAEELDC
lockr graft LISDAELEAIFAEELDC --position 327
lockr graft LISDAELEAIFAEELDC --latch-start 325 --latch-end 359
lockr graft LISDAELEAIFAEELDC --out grafted.pdb
lockr graft --status
```

Binder length is limited to 35 residues. Without `--position`, the command scans the latch window for the best-scoring position. `--status` reports whether PyRosetta and the bundled template are available. Grafting requires the optional dependency.

## API and additional analysis

The local service exposes JSON endpoints for the web app and for scripts. Start it with `lockr serve`, then inspect the exact request and response models at `/docs` or `/openapi.json`.

Beyond the web workflows, the API includes:

- **Parameter sweeps:** vary one supported parameter and return fold-change values across the requested range.
- **Assembly verification:** check a full sequence against a latch window, graft coordinates, spacer/linker positions, protected motifs, and an expected total length. It can also screen candidate variants against a protected region. Positions are 1-based; translate binder-local mutation positions to full-construct positions before submitting them.
- **Sequence checks:** locate the configured SmBiT reporter motif and report its position and sequence warnings.
- **Inverse design:** find parameter combinations that meet a target fold-change or range, with feasible solutions and a robust point when one exists. Mutation-map requests rank possible latch substitutions against a requested free-energy change; those rankings are estimates, not mutation recommendations validated in the lab.
- **Design landscape and slices:** calculate a two-axis grid or a one-axis slice, with optional constraints and temperature-dependence inputs.
- **K_open scenarios and dosing:** the same supporting calculations used by the web interface.

The API is versioned through the request schemas rather than a separate version prefix. Use `/docs` as the source of truth for endpoint paths and accepted fields, since schemas can evolve with the tool.

## Interpreting results

These outputs help prioritize designs; they do not establish that a biosensor will fold, bind its target, or produce a given luminescence signal. The scanner's liability calibration is based on limited ECLIPSE reference data, helix confidence comes from sequence heuristics, and graft scores depend on the template and modeling assumptions. For experiments, confirm affinities and basal opening with measurements and validate candidate constructs in the assay conditions you plan to use.

The thermodynamic engine uses a shared three-state model. Its equations are general to the model, but defaults such as K_open, K_CK, lucKey, and calibration values are scaffold-specific. See [docs/README.md](docs/README.md) for the engine map and links to the thermodynamics, calibration, sequence-liability, and assembly notes.

## Help

```bash
lockr --help
lockr scan --help
lockr fc --help
lockr graft --help
lockr serve --help
```
