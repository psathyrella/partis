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

To install without Docker, you basically just run the commands in the [Dockerfile](../Dockerfile) (ignoring the `FROM` and `COPY` lines).
If you're on a debian variant, run the apt-get install as written.
On other distros, or on macOS, you'll have to figure out the equivalent package names, but after that just swap yum/brew for apt-get.
If you don't care about running simulation, you can ignore the parts that are indicated as simulation-specific.

If you don't have conda installed, follow the full installation instructions [here](https://docs.anaconda.com/anaconda/install/).
Any time you're using conda, it needs to be in your path, for example with `export PATH=<path_to_conda>:$PATH`.
While conda typically works better, you can also install the python packages with pip (you should probably user the `--user` option).

Because conda and pip do not play well with each other, if you've used pip in the past but are going to now use conda, and you won't need pip in the future, you should completely remove `~/.local`.
If you might need pip in the future, you should expect some difficulty with having both conda and pip on the same system.
In most cases you can prevent conda from finding the packages in `~/.local` (and consequently breaking) by setting `export PYTHONNOUSERSITE=True` (see [this issue](https://github.com/conda/conda/issues/448) for some context).
You may also need to `unset LD_LIBRARY_PATH`.

Once conda is installed, run the rest of the commands in the Dockerfile whose lines are marked with `RUN` or `CMD`, except substitute the line `WORKDIR /partis` with `cd partis/`.

In order to avoid polluting your environment, we do not automatically add partis to your path.
Several methods of accomplishing this are described [here](subcommands.md#subcommands).

#### Plotting

For plotting to work with the `partition` action (i.e. if you're setting `--plotdir`), you need to install R (skip if you already have it installed), along with some packages

```
apt-get install xorg libx11-dev libglu1-mesa-dev  # may not be necessary, depending on existing R install, but is at least necessary in docker
apt-get install r-cran-rgl
conda install -y -cr r-rgl r-essentials
unset R_LIBS_SITE
R --vanilla --slave -e 'install.packages(c("bios2mds"), repos="http://cran.rstudio.com/")'

# RPANDA stuff (sorry this is hackey, but at least a.t.m. we need to use this modified src [i.e. not what it is in CRAN])
mkdir -p packages/RPANDA/lib
R CMD INSTALL -l packages/RPANDA/lib packages/RPANDA/
```
