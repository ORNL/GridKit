# GridKit Documentation Target

We use CMake's built-in `FindDoxygen` module to generate a target that will
use Doxygen to build the project documentation.

## Building

The documentation target is excluded from the default set of targets built via
`make` (or `cmake --build .`). It must be designated explicitly with
```sh
make GridKitDocs
```
or
```sh
cmake --build . -t GridKitDocs
```
The generated files can be found in the build directory under `docs/html/` and
the main page is `docs/html/index.html`.

If you wish to inspect the Doxyfile itself, it is generated as
`docs/Doxyfile.GridKitDocs`.

## Notes
The reasoning behind taking this approach is as follows:

1. It is easier to maintain just the few options we need to customize rather
   than a whole Doxyfile (leave what can be generated out of the repo); if
   another option is needed, it's easy to look it up.
2. Since the documentation can be seen as a "product" or "artifact" of the code,
   it makes sense for it to be a buildable "target"
3. Generated files should not be placed in the source directory. The binary
   directory makes sense for this, and the CMake target makes this easy.
