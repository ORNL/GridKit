# GridKit Documentation

GridKit documentation can be built two ways:

1. Read the Docs builds the Sphinx site from `docs/conf.py`.
2. CMake can build the existing Doxygen HTML target.

## Read the Docs Build

The Read the Docs proof of concept is configured by `.readthedocs.yaml`.
Static MyST wrapper pages under `docs/` include the repository Markdown
files directly. Before Sphinx runs, the build generates Doxygen XML and
bundles the JavaScript for the interactive oneline diagrams:

```sh
rm -rf docs/xml
cd docs && doxygen Doxyfile
cd js && npm ci && npm run bundle
```

To test the same flow locally:

```sh
python -m pip install -r docs/requirements.txt
rm -rf docs/xml
cd docs && doxygen Doxyfile
cd js && npm ci && npm run bundle
cd ../..
sphinx-build -T -b html docs docs/_build/html
```

Doxygen XML under `docs/xml`, generated API reference files, the bundled
JavaScript under `docs/_static/js/latkit.bundle.js`, and HTML output are
build artifacts and should not be committed.

## Interactive Oneline Diagrams

Case pages can embed an interactive, WebGPU-rendered oneline with the
`oneline` directive (see `docs/_ext/oneline.py`), pointing at the case
JSON file:

````markdown
```{oneline} ../../../../../examples/PhasorDynamics/Large/WECC/wecc.json
:height: 480px
```
````

The directive converts the case's buses and branches to a render-ready
topology at build time. Buses are placed by `extension.longitude` and
`extension.latitude` when every bus has them, and by an automatic layout
otherwise. Browsers without WebGPU hide the diagram and show the page
unchanged.

## CMake Doxygen Target

The existing CMake target is still available for standalone Doxygen HTML output:

```sh
cmake --build build -t GridKitDocs
```

The generated files are written under the build directory.
