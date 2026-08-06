# This script gets executed in `make test_install`. 

# GridKit install prefix used for CMake dependency finding
$env:GridKit_DIR = $args[0]
echo "GridKit_DIR: ${env:GridKit_DIR}"

$env:COMPILER = $args[1]

# Locate source of the consumer test app
$env:INSTALL_BUILD_CONSUME = "${env:GridKit_DIR}/share/gridkit/examples/GridKitConsumer"
echo "Consumer directory: ${env:INSTALL_BUILD_CONSUME}"

# Create build directory
New-Item -ItemType Directory -Force -Path "${env:INSTALL_BUILD_CONSUME}/build" | Out-Null

Remove-Item -Recurse -Force "${env:INSTALL_BUILD_CONSUME}/build/*"

# Configure consumer project
cmake -B "${env:INSTALL_BUILD_CONSUME}/build" `
    -S "${env:INSTALL_BUILD_CONSUME}" `
    "-DCMAKE_CXX_COMPILER=${env:COMPILER}" `
    "-DGridKit_DIR=${env:GridKit_DIR}"

# Build and install
cmake --build "${env:INSTALL_BUILD_CONSUME}/build"

cmake --install "${env:INSTALL_BUILD_CONSUME}/build"

# Check
cd "${env:INSTALL_BUILD_CONSUME}/build"

# Get build configuration (Debug or Release)
$config = "Debug" # Default
$buildTypeLine = Get-Content "CMakeCache.txt" | Select-String "^CMAKE_BUILD_TYPE:STRING=(.*)$"
if ($buildTypeLine -and $buildTypeLine.Matches.Groups[1].Value) {
    $config = $buildTypeLine.Matches.Groups[1].Value
}

ctest -C $config --output-on-failure

exit $LASTEXITCODE