from spack_repo.builtin.packages.resolve.package import Resolve as BuiltinResolve

from spack.package import *

class Resolve(BuiltinResolve):
    version("gridkit-pinned", commit="022b099623711486c420ea6609bbe4c36263e5bb")
