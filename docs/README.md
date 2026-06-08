# GridKit Documentation

GridKit documentation can be built two ways:

1. Read the Docs builds the Sphinx site from `docs/conf.py`.
2. CMake can build the existing Doxygen HTML target.

## Read the Docs Build

The Read the Docs proof of concept is configured by `.readthedocs.yaml`.
Before Sphinx runs, the build generates Markdown wrapper pages and Doxygen XML:

```sh
python docs/generate_model_docs.py
cd docs && doxygen Doxyfile
```

To test the same flow locally:

```sh
python -m pip install -r docs/requirements.txt
python docs/generate_model_docs.py
cd docs && doxygen Doxyfile
cd ..
sphinx-build -T -b html docs docs/_build/html
```

The generated Sphinx files, Doxygen XML, and HTML output are build artifacts and
should not be committed.

## CMake Doxygen Target

The existing CMake target is still available for standalone Doxygen HTML output:

```sh
cmake --build build -t GridKitDocs
```

The generated files are written under the build directory.
