# To run this week's script.

## Install miniconda

### Windows user

Download [this](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh).
In terminal,
```sh
# **Use WSL 2**
chmod +x Mambaforge-Linux-x86_64.sh
./Mambaforge-Linux-x86_64.sh
```

### Mac users
Download
  - [Mac M1+](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh)
  - [Intel Mac](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh)

In terminal, replace the bracket with your version.
```sh
chmod +x Mambaforge-MacOSX-[arm64 if M1 or x86_64 if Intel].sh
./Mambaforge-MacOSX-[arm64 if M1 or x86_64 if Intel].sh
```

Enter yes for the license and install at the default directory.

**Restart your terminal before proceeding.**

## Create an environment and install packages

```sh
mamba create -y -n qmn python=3.10 jq
conda activate qmn
mamba install -y -c bioconda multiqc fastqc kallisto igv
pip install ffq gget
```

> mamba is the same as conda but much faster.

**You need to run `conda activate qmn` at the start of every terminal session to use the installed packages.**
