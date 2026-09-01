# GridKit development container

1. Install the VS Code **Dev Containers** extension. Also install **Remote -
   SSH** when the checkout is on another machine.
2. Clone GridKit with its submodules:

   ```bash
   git clone --recurse-submodules --branch develop https://github.com/ORNL/GridKit.git
   ```

3. Open the checkout in VS Code and run **Dev Containers: Reopen in
   Container**. The extension builds the image on first use.

When using Podman, set `dev.containers.dockerPath` to `/usr/bin/podman` in VS
Code's user settings before reopening the folder.

Configure, build, and test inside the container:

```bash
cmake --workflow --preset dev
```

For later builds or test runs:

```bash
cmake --build --preset dev
ctest --preset dev
```
