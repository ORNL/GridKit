# GridKit Cases

This directory contains reusable system datasets organized by modeling domain.
Each case package contains a `README.md`, its machine-readable case file, and
any case-specific diagrams or validation images. Case files are stored at:

```text
cases/<Domain>/<case>/<case>.case.json
```

Application study configurations, reference results, and executable examples
remain under `examples/`. A study references its case package rather than
keeping a second copy of the dataset.
