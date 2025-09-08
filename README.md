# doppelganger

We are building new software tools for searching optical counterparts to high-energy transients in the LSST dataset, in real-time.
This repo contains a collection of code, scripts and notebooks produced in our research.
The project codename will likely change in the future. Work in progress.

## Installing

This is only needed if you intend to use the code in `doppelganger/`, or any of our scripts. You don't need to install the package to execute our notebooks: just copy them to the RSP and run them.

Clone this repository and `cd` into it. Create a new environment, activate it, and install with `pip install .`.
To install with optional developer dependencies, in an editable environment, use `pip install -e ".[dev]"`.

### Testing

Run tests from your favourite shell with:
```bash
python -m unittest
```
Don't forget to activate your environment first!

### Scripts

The `tunein` demo script will output parsed notices info from GCN, as they come in real-time.
Activate an environment where this package is installed, then launch `tunein` with:

```python
tunein fermi-gbm
```

For more info check the helper with `tunein --help`.
