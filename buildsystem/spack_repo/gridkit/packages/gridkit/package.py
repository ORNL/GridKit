# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack.package import *

class Gridkit(CMakePackage):
    """Experimental code for prototyping interfaces betwen numerical libraries and network models."""

    homepage = "https://github.com/ORNL/GridKit"
    git = "https://github.com/ORNL/GridKit.git"

    maintainers("nkoukpaizan", "pelesh")

    version("develop", submodules=True, branch="develop")

    variant("asan", default=False, description="Enable/Disable address sanitizer")
    variant("enzyme", default=False, description="Enable/Disable Enzyme")
    variant("ipopt", default=False, description="Enable/Disable Ipopt")
    variant("klu", default=False, description="Enable/Disable KLU")
    variant("resolve", default=False, description="Enable/Disable Re::Solve")
    variant("sundials", default=False, description="Enable/Disable SUNDIALS")
    variant("ubsan", default=False, description="Enable/Disable undefined behavior sanitizer")

    conflicts("+klu", when="~sundials")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    depends_on("enzyme", when="+enzyme")
    depends_on("ipopt", when="+ipopt")
    depends_on("sundials@develop+klu~mpi", when="+sundials+klu")
    depends_on("sundials@develop~klu~mpi", when="+sundials~klu")

    def cmake_args(self):
        args = []
        spec = self.spec

        args.extend(
            [
                self.define_from_variant("GRIDKIT_ENABLE_ASAN", "asan"),
                self.define_from_variant("GRIDKIT_ENABLE_ENZYME", "enzyme"),
                self.define_from_variant("GRIDKIT_ENABLE_IPOPT", "ipopt"),
                self.define_from_variant("GRIDKIT_ENABLE_RESOLVE", "resolve"),
                self.define_from_variant("GRIDKIT_ENABLE_SUNDIALS", "sundials"),
                self.define_from_variant("GRIDKIT_ENABLE_UBSAN", "ubsan"),
            ]
        )

        if self.spec.satisfies("+enzyme"):
            args.extend(
                [
                    self.define("CMAKE_C_COMPILER", f"{self.spec['llvm'].prefix}/bin/clang"),
                    self.define("CMAKE_CXX_COMPILER", f"{self.spec['llvm'].prefix}/bin/clang++"),
                ]
            )

        return args
