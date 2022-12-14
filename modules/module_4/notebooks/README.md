# To run this week's script

# CHECK IF YOU ALREADY HAVE CONDA

```sh
which conda
```

If something gets returned, you already have conda.
### Windows user: conda must be installed in WSL 2.


## If you have conda: Create an environment and install packages

```sh
source ~/.bashrc
conda create -y -n qmn -c conda-forge -c bioconda  python=3.10 multiqc fastqc kallisto igv jq
conda activate qmn
pip install ffq gget
```

### M1 users

`kallisto` won't install. Instead, use `brew` that we installed earlier to install kallisto.
Also install `samtools`, this is something we'll need.
```sh
brew install kallisto samtools
```

> mamba is the same as conda but much faster.

**You need to run `conda activate qmn` at the start of every terminal session to use the installed conda packages.**



## Install miniconda (if you do not already have conda)

### Windows user [use WSL 2]
Download [this](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh).
In terminal,
```sh
# If possible, try using /home/[YOUR USERNAME], the internal WSL drive, as the directory for scripts.
# There is significant disk read/write penalty when using anything in /mnt/c/.
# You can access your WSL drive from Windows Explorer at \\wsl$.
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

For Mac users, also install `homebrew`, a Mac package manager.
```sh
xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

# FAQs
- `kallisto` does not install.
  - If you're a Windows user, please use WSL 2.
  - If you're a Mac user, use Homebrew. Installation guide above.
