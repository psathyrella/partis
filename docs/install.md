### Installation

Partis can be installed either [with Docker](#installation-with-docker) or [from scratch](#installation-from-scratch).

[Using Docker](#installation-with-docker) makes the installation process much easier, because it installs specific versions of each dependency in a controlled environment. The process of running, however, is then a bit less convenient since you're inside Docker (the use of batch systems, for instance, will be difficult).

Installing [without Docker](#installation-from-scratch), on the other hand, will require more effort to install (or more accurately, a higher chance of descent into dependency hell), but once installed it'll be easier to use. The closer your system is to the latest Ubuntu, the easier this process will be (RHEL variants and macOS typically work without too much trouble).

If you want to make simulated samples, you'll also need to [install some R packages](#simulation).

#### Installation with Docker

Docker can be installed using the instructions [here](https://docs.docker.com). If you're new to Docker, start with the official [quick start](https://docs.docker.com/get-started/) guide. Then,

```
sudo docker pull psathyrella/partis  # pull partis image from docker hub
sudo docker run -it --name container-1 -v ~:/host/home psathyrella/partis /bin/bash  # start up a container from this image
```
The `sudo` may not be necessary for some systems. With `docker run`, we create a new container from (i.e. instance of) the partis image. The `-v` mounts your home directory on the host machine to the path `/host/home` inside the container, so we can easily extract results.

Docker's default escape key is ctrl-p, which unfortunately conflicts with common command line shortcuts. To change this put`{"detachKeys": "ctrl-q,q"}` into the file `~/.docker/config.json` which will swap it to ctrl-q.

When finished, you _don't_ want to just exit and `docker run` again -- this will create a new container.
One solution is to run Docker inside [tmux](https://hackernoon.com/a-gentle-introduction-to-tmux-8d784c404340?gi=70388a0228fb) in which case you can just leave the container open.
Alternatively, you can detach from the container with `ctrl-q q` (or whatever you've remapped ctrl-p q to).
To then reattach to this running container, run `docker attach container-1`.

#### Installation from scratch

To install without Docker, you first want to read the [Dockerfile](../Dockerfile) carefully to figure out the dependencies you'll need to install.
We currently recommend using conda, as in the Dockerfile, but [pip]((docs/pip-install.md) has also worked well in the past.

If you don't have conda installed, follow the full installation instructions [here](https://docs.anaconda.com/anaconda/install/), e.g. for [linux](https://docs.anaconda.com/anaconda/install/linux).

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
```

Then run the conda install commands from the [Dockerfile](../Dockerfile).

Finally, open a new terminal, and:

```
source activate partis
pip install colored-traceback dendropy==3.12.3
git clone git@github.com:psathyrella/partis
cd partis
./bin/build.sh
```

#### Simulation

If you want to run the `simulate` action, partis needs some tree-related R packages, which can be installed with

```
conda install -y r-essentials
unset R_LIBS_SITE
R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM", "bios2mds"), repos="http://cran.rstudio.com/")'
```

note that the above ^ process will likely screw up environment variables related to where your system looks for gcc, so after running this in a terminal you'll likely need to specify full paths.
