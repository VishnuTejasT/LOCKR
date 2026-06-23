# LOCKR Thermodynamics Design Tool

This software is a CPU-only thermodynamic engine for LOCKR biosensors. It essentially works by scanning a binding region in the intended latch, and identifies liabilities and problems that can restrict cage-key thermodynamic interactioons. It provideds a fix, and outputs a final K_CK (Cage-Key) score. Then, the user can use the K_CK score from the scanner in the Calculator section to calculate fold-change and other variables to determine whether their LOCKR design is a good fit. Both sections utlize the same physics, so they share the same thermodynamic engine.


# Install Instructions

```git clone <your gitlab repo url>
cd LOCKR
conda create -n igem python=3.10
conda activate igem
pip install -e .
```

# Using the Tool

Run the following:

```
conda activate igem
```