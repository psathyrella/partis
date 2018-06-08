### Installation

Partis can be installed either [with Docker](#installation-with-docker) or [from scratch](#installation-from-scratch).

[Using Docker](#installation-with-docker) makes the installation process much easier, because it installs specific versions of each dependency in a controlled environment. The process of running, however, is then a bit less convenient since you're inside Docker (the use of batch systems, for instance, will be difficult).

Installing [without Docker](#installation-from-scratch), on the other hand, will require more effort to install (or more accurately, a higher chance of descent into dependency hell), but once installed it'll be easier to use. The closer your system is to the latest Ubuntu, the easier this process will be (RHEL variants and macOS typically work without too much trouble).

#### Installation with Docker

Docker can be installed using the instructions [here](https://docs.docker.com). If you're new to Docker, start with the official [quick start](https://docs.docker.com/get-started/) guide. Then,

```
sudo docker pull psathyrella/partis  # pull partis image from docker hub
sudo docker run -it -v ~:/host/home psathyrella/partis /bin/bash  # start up a container from this image
```
The `sudo` may not be necessary for some systems. With `docker run`, we create a new container from (i.e. instance of) the partis image. The `-v` mounts your home directory on the host machine to the path `/host/home` inside the container, so we can easily extract results.

Docker's default escape key is ctrl-p, which unfortunately conflicts with common command line shortcuts. To change this put`{"detachKeys": "ctrl-q,q"}` into the file `~/.docker/config.json` which will swap it to ctrl-q.

When finished, make sure to detach using `ctrl-q q` (or whatever you've remapped ctrl-p q to), rather than exiting (e.g. with ctrl-d or `exit`), otherwise you'll have to make a new container and recompile.
To pick up again from were you left off, you typically want to reattach (if you instead exit and `docker run` again, you'll make a new container).
To reattach to a container from which you've detached, list running containers with `docker ps`, get the right container id, then `docker attach <CONTAINER_ID>`.
This is enabled by the `-it` and `/bin/bash` options we used above for `docker run`: these allocate a pseudo-tty, keep STDIN open, and run bash instead of the default command.

#### Installation from scratch

The docker-free approach is frequently a lot faster (if you already have numpy and scipy installed, for instance, you don't have to wait for docker to compile and install them from scratch).
You also don't have to deal with the additional complications of being inside docker, perhaps most importantly docker's likely incompatibility with your batch queueing system.

You're basically installing the dependencies listed in the [Dockerfile](https://github.com/psathyrella/partis/blob/master/Dockerfile), plus a few extras.
We have had good experiences using both [conda](#conda) (recommended), and [pip]((https://github.com/psathyrella/partis/blob/master/docs/pip-install.md).
Both are described below.

##### conda

If you don't have conda installed, follow the full installation instructions for your system [here](https://docs.anaconda.com/anaconda/install/), e.g. for [linux](https://docs.anaconda.com/anaconda/install/linux).

Any time you're using conda, it needs to be in your path, for instance:
```
export PATH=<path_to_conda>:$PATH
```
If you've used pip in the past, but you won't need it in the future, you should also completely remove `~/.local`.
If you might need pip in the future, you should expect some difficulty with having both conda and pip on the same system.
In most cases you can prevent conda from finding the packages in `~/.local` (and consequently breaking) by setting `export PYTHONNOUSERSITE=True` (see [this issue](https://github.com/conda/conda/issues/448) for some context).
You may also need to `unset LD_LIBRARY_PATH`.

Then make a conda environment for partis:

```
conda create -y -n partis
source activate partis  # do this _every_ time you start a new terminal
conda install -y biopython scikit-learn cmake gsl openblas pandas psutil pysam r-essentials scons seaborn pyyaml
conda install -y -c biocore mafft
```

Then start a new terminal, and:

```
source activate partis
pip install colored-traceback dendropy==3.12.3
unset R_LIBS_SITE
R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM", "bios2mds"), repos="http://cran.rstudio.com/")'
```

Then clone and compile:

```
git clone git@github.com:psathyrella/partis
cd partis
./bin/build.sh
```
