# GridKit image

The runtime image contains the complete installed GridKit tree, built with
Apache Arrow so the `arrow` and `arrow_stream` monitor formats are available.
It has no GridKit entrypoint: select an installed executable after the image
name, or run the image interactively to open a shell.

Run all commands below from the repository root.

## Build

Build the development image first if it is not already available:

```bash
/usr/bin/podman build -t gridkit-dev:latest .devcontainer
```

Then build the application image, tagged for what it carries (`latest`, or a
dated tag such as `arrow-2026-09-09` that consumers can pin):

```bash
/usr/bin/podman build \
  --build-arg BUILDER_IMAGE=gridkit-dev:latest \
  --build-arg BUILD_JOBS=4 \
  -f docker/Dockerfile \
  -t gridkit:latest \
  .
```

## VS Code tasks

Install the recommended **VSCode Action Buttons** extension and reload VS Code.
Three status-bar buttons build `gridkit-dev:latest`, build `gridkit:<tag>`,
or confirm and push `gridkit:<tag>` to `ghcr.io/lukelowry/gridkit:<tag>`; the
build and publish buttons prompt for the tag, defaulting to `latest`.
Run **Refresh Action Buttons** from the Command Palette after changing their
configuration.

These tasks run on the machine hosting the VS Code workspace and require
`/usr/bin/podman`. Run them from a local or Remote-SSH window, not from inside
the development container. The publish action asks you to type `PUBLISH`
before pushing. Authenticate before publishing:

```bash
/usr/bin/podman login ghcr.io
```

## Run

Create a writable output directory:

```bash
mkdir -p docker/output
```

Run the two-bus example with `DynamicSimulation`:

```bash
/usr/bin/podman run --rm \
  --volume "$PWD/examples/PhasorDynamics/Tiny/TwoBus/Basic:/input:ro" \
  --volume "$PWD/docker/output:/work" \
  gridkit:latest \
  DynamicSimulation /input/TwoBusBasic.solver.json
```

Select a different installed application in the same way:

```bash
/usr/bin/podman run --rm \
  --volume "$PWD/examples/PhasorDynamics/Medium/NewEngland:/input:ro" \
  --volume "$PWD/docker/output:/work" \
  gridkit:latest \
  ContingencyAnalysis /input/newengland.solver.json
```

Input is mounted read-only at `/input`. Generated files are written under
`/work`, which is mapped to `docker/output` in these examples.

Open a shell to inspect the image or use several GridKit tools:

```bash
/usr/bin/podman run --rm --interactive --tty gridkit:latest
```

On SELinux hosts, add the appropriate `:z` or `:Z` bind-mount label. With
rootful Docker, add `--user "$(id -u):$(id -g)"` when host output files should
be owned by the invoking user.
