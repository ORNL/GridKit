{
  description = "gridkit development shell";

  inputs = {
    enzyme-src = {
      url = "github:EnzymeAD/Enzyme/v0.0.303";
      flake = false;
    };
    nlohmann-json-src = {
      url = "github:nlohmann/json/d70e46bc65a21a8eb1ef0e3e9ba707b7b9b78378";
      flake = false;
    };
    sundials-src = {
      url = "github:LLNL/sundials/v7.9.0";
      flake = false;
    };
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs =
    {
      self,
      nixpkgs,
      sundials-src,
      enzyme-src,
      nlohmann-json-src,
      ...
    }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
      ];
      forAllSystems = nixpkgs.lib.genAttrs systems;
    in
    {
      packages = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
          llvmPackages = pkgs.llvmPackages_21;
          enzyme = (pkgs.enzyme.override { inherit llvmPackages; }).overrideAttrs (old: {
            version = "0.0.303";
            src = enzyme-src;
            patches = (old.patches or [ ]) ++ [ ./nix/enzyme-optimization-levels.patch ];
          });
          sundials79 = pkgs.sundials.overrideAttrs (old: {
            version = "7.9.0";
            src = sundials-src;

            cmakeFlags = (lib.filter (flag: !lib.hasPrefix "-DEXAMPLES_INSTALL_PATH=" flag) old.cmakeFlags) ++ [
              (lib.cmakeBool "SUNDIALS_ENABLE_EXAMPLES_INSTALL" true)
              (lib.cmakeFeature "SUNDIALS_EXAMPLES_INSTALL_PATH" "${placeholder "examples"}/share/examples")
            ];
          });
          inherit (pkgs) lib;
          nlohmannJson = pkgs.nlohmann_json.overrideAttrs {
            version = "3.12.0-unstable-d70e46b";
            src = nlohmann-json-src;
          };
          gridkit = llvmPackages.stdenv.mkDerivation {
            pname = "gridkit";
            version = "unstable";
            src = self;

            nativeBuildInputs = [
              pkgs.cmake
              pkgs.gfortran
              pkgs.pkg-config
              llvmPackages.llvm
              enzyme
            ];

            buildInputs = [
              sundials79
              pkgs.ipopt
              (lib.getDev pkgs.ipopt)
              pkgs.cppzmq
            ];

            # link in what would otherwise be a submodule
            postPatch = ''
              mkdir -p third-party/magic-enum third-party/nlohmann-json
              ln -s ${pkgs.magic-enum}/include third-party/magic-enum/include
              ln -s ${nlohmannJson}/include third-party/nlohmann-json/include
              ln -s ${nlohmannJson}/include third-party/nlohmann-json/single_include
            '';

            cmakeBuildType = "Release";
            cmakeFlags = [
              (lib.cmakeBool "GridKit_ENABLE_ENZYME" true)
              (lib.cmakeFeature "ENZYME_DIR" "${enzyme}")
              (lib.cmakeBool "GridKit_ENABLE_SUNDIALS" true)
              (lib.cmakeBool "GridKit_ENABLE_IPOPT" true)
              (lib.cmakeFeature "SUNDIALS_DIR" "${sundials79}")
              (lib.cmakeFeature "IPOPT_DIR" "${pkgs.ipopt}")
              (lib.cmakeFeature "IPOPT_INCLUDE_DIR" "${lib.getDev pkgs.ipopt}/include/coin-or")
            ];

            enableParallelBuilding = true;
            doCheck = true;
            checkPhase = "ctest --output-on-failure";

            meta = with lib; {
              description = "Modeling framework for power systems simulations and analysis";
              homepage = "https://github.com/ORNL/GridKit";
              platforms = platforms.unix;
            };
          };
        in
        {
          default = gridkit;
          gridkit = gridkit;
          enzyme = enzyme;
        }
      );

      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
        in
        {
          default = (pkgs.mkShell.override { stdenv = pkgs.llvmPackages_21.stdenv; }) {
            inputsFrom = [ self.packages.${system}.gridkit ];
            inherit (self.packages.${system}.gridkit) cmakeFlags;
            packages = [
              pkgs.ninja
            ];
          };
        }
      );
    };
}
