##################
### creates my basic conda env: run on compute node
### choose name and kfs

#!/bin/bash
echo "here we go"

#envs_dir=/scratch/isatkaus/flash/conda-envs/
### flash w 256K DoM
envs_dir=/home/isatkaus/flash-w-dom/conda-envs/

### package manager: conda or mamba
### mamba is much faster, unless conda is 23.x (uses libmamba as its default solver)
PKG_MGR=mamba

### everything from conda
env_name=h-py312-basic

#env_name=rev-min-cray-hdf5

if [ -d "$envs_dir$env_name" ]; then
    read -p "Environment '$envs_dir$env_name' exists. Overwrite? [y/N] " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 1
    else
        echo "Overwriting environment '$envs_dir$env_name'."
        rm -rf "$envs_dir$env_name" 
    fi
fi

cd $envs_dir

echo -e "\nCreating conda environment: $env_name"

if [[ "$PKG_MGR" == "mamba" ]]; then
    ml mamba
else
    ml conda
fi
echo "Using package manager: $PKG_MGR"
#ml cray-hdf5/1.12.2.9 (if using this, need to install h5py with pip)
#echo "using hdf5 module: $HDF5_DIR"

$PKG_MGR create --prefix ./$env_name python=3.12 numpy pandas pyyaml jupyter plotly h5py poppler pillow --yes
#$PKG_MGR create --prefix ./$env_name python=3.12 --yes
conda activate $envs_dir$env_name


### install additional:
# poppler is for pdftotext
#conda install -c conda-forge poppler
$PKG_MGR install -y scipy pyarrow


### install additional via pip
#echo $(which pip)
echo "installing via pip (pip version: $(pip --version))"
pip install tables matpowercaseframes==2.1.0

#pip install --no-cache-dir --no-binary=h5py h5py==3.14.0
#echo "lld'ing h5py: $(ldd $envs_dir$env_name/lib/python3.10/site-packages/h5py/h5.cpython-310-x86_64-linux-gnu.so | grep hdf5)"


#pip install nrel-rev==0.14.5 --no-cache-dir
#pip install darshan==3.4.7 --no-cache-dir

#### additional "max" packages:
#mamba install pandas pyyaml jupyter plotly --yes


echo -e "\ninstalled $env_name in: \n$envs_dir$env_name"
echo -e "\nTo activate this environment, use:"
echo -e "ml conda"
echo -e "conda activate $envs_dir$env_name"

