from spack_repo.builtin.packages.sundials.package import Sundials as BuiltinSundials

from spack.package import *

class Sundials(BuiltinSundials):
    version("7.9.0", tag="v7.9.0", commit="312fc0f3684f27209ca9dc9249d194436eb41a7a")
    version("7.8.0", tag="v7.8.0", commit="aedc088437064dd55b35c000145f7f5db6ee49e3")
