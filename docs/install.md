[Up to table of contents](contents.md)

### Installation methods

Partis can be installed via:
  - [pip and pipx](#installation-with-pip-or-pipx): recommended
  - [additional steps for simulation installation](#simulation)
  - [with Docker](#installation-with-docker)
  - [plotting dependencies](#plotting)

#### Installation with pip or pipx

First install system dependencies:
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install python3 python3-pip python-is-python3 pipx build-essential cmake libgsl-dev libyaml-cpp-dev scons mafft ncurses-base ncurses-bin

# macOS (with Homebrew)
brew install python3 pipx cmake gsl yaml-cpp scons mafft
```
Note, however, that macOS 14 or newer (ARM64) cannot run or compile the ig-sw dependency due to incompatibility with emmintrin/SSE2.
Until we've managed to resolve this, then, you'll unfortunately be limited to either Linux or macOS 13 or older.

Install from PyPI (recommended):
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: partis-env\Scripts\activate
pipx install partis-bcr
```
Alternatively, you can install from source if you think you'll want to modify the source code.
```bash
git clone https://github.com/psathyrella/partis.git
cd partis

# Initialize essential submodules
git submodule update --init packages/ham packages/ig-sw

# Create virtual environment and install
python -m venv .venv
source .venv/bin/activate  # On Windows: partis-env\Scripts\activate
pip install -e .
```

##### Testing the installation

```bash
# Run quick test
partis-test.py --quick

# Run more complete test (still without simulation, since the base/pip install above doesn't install simulation requirements)
partis-test.py --no-simu
partis-test.py --paired --no-simu
```

#### Simulation

A variety of partis simulation options will work using the base install described above; however for full functionality you'll need to install R and some R packages, as well as potentially compile bio++.
An overview of simulation options is [here](subcommands.md#simulate), and further information is available by running `partis simulate --help`.
Briefly, partis simulation can generate trees using the R packages TreeSim and TreeSimGM, which requires installing them (along with R).
Alternatively, if you have a list of your own newick trees in a file you can pass this in with `--input-simulation-treefname`, in which case you can avoid R installation entirely.
To install R, one option is:

```
# Ubuntu/Debian
sudo apt install r-base

# macOS (with Homebrew)
brew install r
```
To then install the R packages run:
```
R --slave -e 'dir.create(c(Sys.getenv("R_LIBS_USER")), recursive=TRUE); install.packages(c("TreeSim", "TreeSimGM", "geiger", "MASS"), repos="http://cran.rstudio.com/", dependencies=TRUE, lib=c(Sys.getenv("R_LIBS_USER")))'
```
This command creates a directory for local libraries in the user's home  in order to avoid system-wide install.

The full/default mutation model also requires compilation of an updated (development) version of bpp whose source code is in `packages/bpp-newlik/`.
However, if you don't much care about mutation model accuracy (in particular, if you don't care about modeling different rates between different bases e.g. A->C vs A->T), you can set `--no-per-base-mutation` and it'll instead use an older, already-compiled bpp version in `packages/bpp`.
To use the full mutation model first install some dependencies:

```
# Ubuntu/Debian
sudo apt install libblas-dev liblapack-dev libboost-all-dev libeigen3-dev libsuitesparse-dev libfftw3-dev

# macOS (with Homebrew)
brew install openblas lapack boost eigen suite-sparse fftw

```
And then compile the bpp newlik version from source:
```
./bin/build.sh with-simulation  # this is usually slow, maybe 30 minutes
```

#### Installation with Docker

Docker can be installed using the instructions [here](https://docs.docker.com). If you're new to Docker, start with the official [quick start](https://docs.docker.com/get-started/) guide. Then,

```
sudo docker pull quay.io/matsengrp/partis
sudo docker run -it --user=root --name container-1 -v ~:/host/home quay.io/matsengrp/partis /bin/bash
```
The `sudo` may not be necessary for some systems. With `docker run`, we create a new container from (i.e. instance of) the partis image. The `-v` mounts your home directory on the host machine to the path `/host/home` inside the container, so we can easily extract results.
The `--user=root` is a temporary fix: mamba uses the `mambauser` user inside docker, but it is difficult to mount the host filesystem in a way such that it is writeable by the mambauser.

Docker's default escape key is ctrl-p, which unfortunately conflicts with common command line shortcuts. To change this put`{"detachKeys": "ctrl-q,q"}` into the file `~/.docker/config.json` which will swap it to ctrl-q.

When finished, you probably don't want to just exit and `docker run` again -- this will create a new container from scratch, whereas you usually want to re-enter the previous container.
One solution is to run Docker inside [tmux](https://hackernoon.com/a-gentle-introduction-to-tmux-8d784c404340?gi=70388a0228fb) in which case you can just leave the container open.
Alternatively, you can detach from the container with `ctrl-q q` (or whatever you've remapped ctrl-p q to).
To then reattach to this running container, run `docker attach container-1`.
If you don't plan on reattaching to a container, you should run with `--rm` so it is removed when you exit; otherwise you'll use up a lot of disk space with accumulated old containers (view with `sudo docker ps -a`).

#### Plotting

In order to turn on plotting (typically by setting `--plotdir`) you'll have to install some extra dependencies:
```bash
pip install partis-bcr[plotting]
```
or if installing from source:
```bash
pip install -e .[plotting]
```

##### MDS Plotting with R

In order to make the [MDS plots](plotting.md#partition-plots), you need R and the bios2mds package installed:

```bash
apt-get install xorg libx11-dev libglu1-mesa-dev r-cran-rgl  # may not be necessary, depending on existing R install, but is at least necessary in docker
conda install -y -cr r-rgl r-essentials
unset R_LIBS_SITE
R --vanilla --slave -e 'install.packages(c("class", "bios2mds"), dependencies=T, repos="http://cran.rstudio.com/")'  # this may be all the deps: c("picante", "pspline", "deSolve", "igraph", "TESS", "fpc", "pvclust", "corpcor", "phytools", "mvMORPH", "geiger", "mvtnorm", "glassoFast", "Rmpfr")
```

