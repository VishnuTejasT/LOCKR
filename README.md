# LOCKR Thermodynamics Design Tool

This software is a CPU-only thermodynamic engine for LOCKR biosensors. It essentially works by scanning a binding region in the intended latch, and identifies liabilities and problems that can restrict cage-key thermodynamic interactioons. It provideds a fix, and outputs a final K_CK (Cage-Key) score. Then, the user can use the K_CK score from the scanner in the Calculator section to calculate fold-change and other variables to determine whether their LOCKR design is a good fit. Both sections utlize the same physics, so they share the same thermodynamic engine.


# Install Instructions

```
git clone https://github.com/VishnuTejasT/LOCKR.git
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

If you need help with running the software, please run this command:
'''
lockr --help
'''

Then, you have two options to actually run the program. 
   1. The first option is to run the tool through the web app. To do this, run the follwoing command:
   ```
   lockr serve
   ```
   Then, open the link that is provided

   2. The second option is to run the tool through the integrated CLI. There are 2 commands associated with this:
   
   "lockr scan" runs the scanner portion, while 'lockr fc" runs the fold-chnage calculation portion.
   ```
   lockr scan SEQUENCE
   lockr scan SEQUENCE --ph 7.4     #User can specify the pH of the hypothetical/inteded solution.
   lockr scan SEQUENCE --window 1:17 
   lockr scan SEQUENCE --preserve 1,2,11,12,15     #User can choose to preserve specific aa residues on the binder, to prevent their mutation.
   lockr scan SEQUENCE --policy conservative / neutralizing    # User can choose b/w conservation mutations or neutralizing mutations of liable residues in the binder.
   lockr scan SEQUENCE --suggest / --no-suggest
   lockr scan SEQUENCE --json    # Accepts FASTA files, raw files, or even a mix of both!
   lockr scan --file sequences.fasta      # Accepts FASTA files, raw files, or even a mix of both!
   lockr scan --file sequences.fasta --json     # Accepts FASTA files, raw files, or even a mix of both!
   ```
   AND 
   ```
   lockr fc --k-ck FLOAT --k-open FLOAT --pull FLOAT --luckey FLOAT
   lockr fc --k-ck FLOAT --k-open FLOAT --pull FLOAT --luckey FLOAT --k-target FLOAT --target FLOAT
   lockr fc --k-ck FLOAT --k-open FLOAT --pull FLOAT --luckey FLOAT --json
   ```

   Here are additional commands for more help:

   ```
   lockr --help
   lockr scan --help
   lockr fc --help
   lockr serve --help
   python3 -m lockr.cli          # Another entry point that also works
   ```
   

