from spack_repo.builtin.packages.enzyme.package import Enzyme as BuiltinEnzyme

from spack.package import *

class Enzyme(BuiltinEnzyme):
    version("0.0.196", sha256="2b9cfcb7c34e56fc8191423042df06241cf32928eefbb113ac3c5199e3361cb2")
    version("0.0.186", sha256="125e612df0b6b82b07e1e13218c515bc54e04aa1407e57f4f31d3abe995f4714")
