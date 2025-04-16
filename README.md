# doppelganger

This is a repo collecting tools for searching LSST data for optical counterparts to high-energy events.
Work in progress.

## Installing

Clone this repository and `cd` into it. Create a new environment, activate it, and install with `pip install .`.
To install with optional developer dependencies, in an editable environment, use `pip install -e ".[dev]"`.

## Testing

Run tests from your favourite shell with:
```bash
python -m unittest
```
Don't forget to activate your environment first!

## Scripts

The `tunein` demo script will output parsed notices info from GCN, as they come in real-time.
Activate an environment where this package is installed, then launch `tunein` with:

```python
tunein fermi-gbm
```

For more info check the helper with `tunein --help`.