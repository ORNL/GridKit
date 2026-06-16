# GridKit Documentation

GridKit documentation can be built two ways:

1. Read the Docs builds the Sphinx site from `docs/conf.py`.
2. CMake can build the existing Doxygen HTML target.

## Read the Docs Build

The Read the Docs proof of concept is configured by `.readthedocs.yaml`.
Static MyST wrapper pages under `docs/` include the repository Markdown
files directly. Before Sphinx runs, the build generates Doxygen XML:

```sh
rm -rf build/docs/doxygen
cd docs && doxygen Doxyfile
```

To test the same flow locally:

```sh
python -m pip install -r docs/requirements.txt
rm -rf build/docs/doxygen
cd docs && doxygen Doxyfile
cd ..
sphinx-build -T -b html docs docs/_build/html
```

Doxygen XML under `build/docs/doxygen`, generated API reference files, and HTML
output are build artifacts and should not be committed.

## CMake Doxygen Target

The existing CMake target is still available for standalone Doxygen HTML output:

```sh
cmake --build build -t GridKitDocs
```

The generated files are written under the build directory.
