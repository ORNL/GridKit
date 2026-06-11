One way to build GridKit and its dependencies is by using Spack. 
We include Spack as a submodule to the Gridkit repository and use it for our CI builds. 
For more information on how to use Spack, please refer to their [documentation](https://spack-tutorial.readthedocs.io/en/latest).

Example steps to build using Spack with an environment created on the fly:
```
git submodule update --init buildsystem/spack
source ./buildsystem/spack/share/spack/setup-env.sh
spack env create gridkit
spack env activate -p gridkit
spack repo add buildsystem/spack_repo/gridkit
spack develop --path=$(pwd) gridkit@develop
spack compiler find
spack add gridkit+enzyme+ipopt+klu+sundials ^sundials@develop
spack concretize -f
spack install
spack env deactivate
```

As GridKit matures, we plan to expand on this build approach and provide Spack
configurations for different systems of interest.
