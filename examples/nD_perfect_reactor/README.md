# Perfectly Stirred Reactor

Reference:
> G. B. Skinner and G. H. Ringrose, “Ignition Delays of a Hydrogen—Oxygen—Argon Mixture at Relatively Low Temperatures”, J. Chem. Phys., vol. 42, no. 6, pp. 2190–2192, Mar. 1965. Accessed: Oct. 13, 2024.

<img src="result.png" height="MAX_HEIGHT"/>

## Validation

After running the simulation, compare MFC species mass fractions and induction
time against a Cantera 0-D ideal-gas reactor reference:

```bash
python analyze.py
```

This reads the Silo output, runs an equivalent Cantera reactor, prints the
induction times (Skinner et al. / Cantera / (Che)MFC), and saves `plots.png`.
All dependencies are installed automatically by the MFC toolchain.
