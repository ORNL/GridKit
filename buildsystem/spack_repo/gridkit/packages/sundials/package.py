from spack_repo.builtin.packages.sundials.package import Sundials as BuiltinSundials

from spack.package import *

class Sundials(BuiltinSundials):
    version("7.9.0", tag="v7.9.0", commit="312fc0f3684f27209ca9dc9249d194436eb41a7a")
    version("7.8.0", tag="v7.8.0", commit="aedc088437064dd55b35c000145f7f5db6ee49e3")
    version("partitioned-dae", branch="feature/partitioned-dae-stepper")

    # The branch registers ark_partitioned_dae.c without its .out file, which breaks install.
    # Remove once the branch includes examples/arkode/C_serial/ark_partitioned_dae.out.
    @when("@partitioned-dae")
    def patch(self):
        filter_file(
            r'^\s*"ark_partitioned_dae\.c\\;\\;develop"\s*$',
            "",
            join_path("examples", "arkode", "C_serial", "CMakeLists.txt"),
        )
