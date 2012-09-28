SMCTC 1.0: Sequential Monte Carlo Template Class
----------------------------------------------------
Last modified:       2009-04-15
Author / Maintainer: Adam M. Johansen (a.m.johansen@warwick.ac.uk)
Website:             
http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic/johansen/smctc/



This file provides a very brief overview of the SMCTC software. It is not
intended to replace the documentation available both in the form of the
automatically generated reference guide or the accompanying article.

1. Overview:

Sequential Monte Carlo methods are a very general class of Monte Carlo methods
for sampling from sequences of distributions. Simple examples of these
algorithms are used very widely in the tracking and signal processing
literature. Recent developments illustrate that these techniques have much
more general applicability, and can be applied very effectively to statistical
inference problems. Unfortunately, these methods are often perceived as being
computationally expensive and difficult to implement. This software and the
accompanying article (SMCTC: Sequential Monte Carlo in C++) seeks to
address both of these problems.

A C++ template class library for the efficient and convenient
implementation of very general Sequential Monte Carlo algorithms is
provided. Two example applications are provided: a simple particle filter for
illustrative purposes and a state-of-the-art algorithm for rare event
estimation.

2. Installation:

2.1 Prerequisites
Use of this software requires a working C++ implementation and that the GNU
Scientific Library (GSL) is installed and working. It is necessary to either have
a working make environment or to make use of the proprietary project
and solution file approach provided by Microsoft Visual Studio. If the libraries are not
available in the directories searched by the compiler then it may be
necessary to modify the Makefile so that the appropriate search paths are
supplied to the compiler and linker.

2.2 Installation Procedure
This section assumes that the Makefile is being used from a command line
interpreter of some description. Windows users may prefer to make use of the
Visual studio solution and project files which are provided. See the main
documentation and the project file itself for details of this alternative.

In a sense, installation is a minimal procedure. The software should be placed
in an appropriate directory and the command
make libraries
within the top-level directory will compile the binary portion of the library
and copy it to the lib subdirectory.

Several optional steps will make it easier to make use of the library:

After compiling the library there will be a static library named libsmctc.a
within the lib subdirectory. This should either be copied to your preferred
library location  (typically /usr/lib on a Linux system) or its locations must
be specified every time the library is used.

The header files contained within the include subdirectory should be copied to
a system-wide include directory (such as /usr/include) or it will be necessary
to specify the location of the SMCTC include directory whenever a file which
makes use of the library is compiled.

3. Using SMCTC:

In order to use the library one simply needs to write C++ code which makes use
of the templated-objects provided by the library and then to link the
resulting object files against the SMCTC library as well as the GSL. The
example programs can be compiled by entering the command:
make examples
in the top-level directory.

The use of a Makefile is not essential at this stage; compilation proceeds
rather simply in the case of self-contained applications. The following
commands, for example, are sufficient to compile the example program described
in section 5.1:

g ++ -I ../../ include -c pfexample . cc pffuncs . cc
g ++ pfexample . o pffuncs . o -L ../../ lib - lsmctc - lgsl - lgslcblas - opf

It is, of course, advisable to include some additional options to encourage
the compiler to optimise the code as much as possible once it has been debugged.

Please see the associated documentation for more details.
