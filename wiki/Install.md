Currently, StochHMM source code is available by cloning the git repository (Github or Google Code).   The git repositories are updated in tandem.  

StochHMM source code is also available as a tar-gzipped archive. 


#Getting StochHMM

## Downloading from Google Code
####[Download StochHMM Archive v0.34](https://stochhmm.googlecode.com/files/StochHMM_v034.tar.gz)

After downloading follow the compiling instructions


###Cloning from Github
Change to a directory were you would like to clone the StochHMM repository.   It will create a folder called StochHMM in your current directory
```
$ git clone git://github.com/KorfLab/StochHMM.git
```

###Cloning from Google Code
```
$ git clone https://lottpaul@code.google.com/p/stochhmm/ 
```


# Compiling StochHMM using GNU make(Linux, Mac OS X, Windows)
StochHMM uses the GNU automake tools to generate the make file and compile the library and executable
If you downloaded the archive, please extract the files.

```
$ cd [StochHMM directory]
$ ./configure
$ make
```
This uses the makefile to compile StochHMM binary and static library. 

The StochHMM executable will be compiled in the projects root directory. On Linux and Mac OS X systems, the executable will be ```stochhmm```.  On Windows systems, the executable will be ```stochhmm.exe```.  The library libstochhmm.a will be located in the src/ directory.

StochHMM has been compiled and tested successfully on Ubuntu 12.04, Mac OS X, and Windows 


[Next: Running StochHMM](Running-StochHMM)