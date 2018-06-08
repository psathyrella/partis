
##### pip

For Ubuntu (tested on 16.04):

```
sudo apt-get install python-pip cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft r-base
pip install --user numpy scipy scikit-learn matplotlib pandas biopython dendropy==3.12.3 pysam pyyaml seaborn colored_traceback psutil
R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM", "bios2mds"), repos="http://cran.rstudio.com/")'
```
If you're running an old RHEL variant, with a (very) old gcc, you may need to either update gcc, or swap out `-Ofast` for something else in several places.

For macOS, you pretty much swap out the `apt-get`/`pip` calls for a combination of `brew` and `pip`.
For instance, depending what is already installed on your system (the Ubuntu list above is more comprehensive) this is likely sufficient:

```
pip install --user colored-traceback pysam cmake pyyaml psutil  # adding --egg may help if it fails
brew install gsl scons
```

Then clone and compile:

```
git clone git@github.com:psathyrella/partis
cd partis
./bin/build.sh
```
