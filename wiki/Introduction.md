StochHMM is an easy to use and flexible hidden Markov model (HMM) application and C++ library. It implements HMMs from a simple user-defined text file, thereby allowing researchers to focus on the design and optimization of the model rather than spending precious time implementing HMM algorithms. 

StochHMM implements standard HMM and hidden semi-Markov model architectures and algorithms. It grants researchers the power to integrate additional datasets in their HMM to improve predictions. Finally, it adapts HMM algorithms to provide stochastic decoding giving researchers the ability to explore and rank sub-optimal predictions.

**Why StochHMM?**

***
	
*  Current HMM libraries/toolboxes lack the abilities to integrate additional datasets and to implement stochastic decoding algorithms. [1,2,3,4,5] StochHMM extends standard HMM architecture to allow integration of a priori/related data at multiple points within the HMM framework.

*  Provides a HMM application and C++ library that allows HMMs to be rapidly developed, flexible, and accessible to researchers regardless of programming experience.

*  Handles ambiguous characters in a state-dependent manner.

*  Implements stochastic decoding algorithms to allow for sub-optimal predictions.



**Integrating Additional Data**

***

StochHMM allows the integration of a priori or related supporting data in multiple ways:

1. Implementing multiple emission states. Multiple datasets of supporting evidence can be integrated either as joint or independent emission distributions in the analysis.

2. Linking the HMM framework to user-defined functions or existing applications. Such linking provides the ability to integrate results from existing function, databases and/or applications into either the emissions or transitions of the states. This allows researchers to incorporate, for example, a Pfam or BLAST database query into the HMM to improve the prediction.

3. Weighting or explicit definition of state path. This allows researchers to restrict or weight predictions within regions of the sequence by what is already known. 



**Stochastic Sampling**

***

StochHMM implements standard HMM decoding algorithms (Viterbi, Posterior, and Nth-best Viterbi) and extends them to provide stochastic decoding (stochastic Viterbi, stochastic Posterior, and stochastic Forward). 

Stochastic decoding allows exploration of suboptimal paths. This provides researchers the ability to develop a gene-finder, which can predict and rank multiple splice isoforms. Additionally, StochHMM gives researchers the ability to create heat maps or simplified weighted directed acyclic graphs from the stochastic traceback table. 


**Additional Features**

***

We have added many additional features to increase its utility in diverse research projects. StochHMM supports: explicit state durations, mixed higher-order emissions, lexical transitions (Mealy machine), user-defined ambiguous character definition, and emission-dependent handling of ambiguous characters.


**StochHMM Project**

***

We have successfully used StochHMM in a diverse set of research projects.  We developed SkewR, a HMM to predict R-loop formation in the human genome in collaboration with Paul Ginno.[6] In addition, we worked with Diane Schroeder to rapidly develop HMMs to identify large-scale methylation domains in cell lines (SH-SY5Y and IMR90)[7] and full-term human placenta tissue.[8] 


**Future Direction**

***

Currently, we are optimizing StochHMM to improve speed and memory usage. We are also working to provide additional examples, documentation, multiple training tools, and an optional GUI interface. We have plans on adding SWIG language bindings to allow StochHMM library utilization in additional languages (Perl, Java, Python). To assist researchers, a forum has been setup on Google Groups to provide advice and assistance in the development of their own HMMs using StochHMM.


**Availability**

***

We are providing StochHMM as a free standalone application and C++ library to give researchers the ability to rapidly develop HMMs.

StochHMM is provided as source code and compiles on GNU/CLANG compilers under Windows, Mac OS X, and Linux. We are providing StochHMM under the MIT open source license to increase accessibility and to give researchers the ability to use it in derivative works without restrictions.

**References**

***

1. Schliep, A., Georgi, B., Rungsarityotin, W., Costa, I. G. & Schonhuth, A. The general hidden markov model library: Analyzing systems with unobservable states. Proceedings of the Heinz-Billing-Price 121–136 (2004).
2. Sand, A., Pedersen, C., Mailund, T. & Brask, A. HMMlib: A C++ Library for General Hidden Markov Models Exploiting Modern CPUs. CORD Conference Proceedings 126–134 (2010).doi:10.1109/PDMC-HiBi.2010.24
3. Lunter, G. HMMoC--a compiler for hidden Markov models. BIOINFORMATICS 23, 2485–2487 (2007).
4. Himmelman, L. R - HMM package.  (2010) Available from http://cran.r-project.org/web/packages/HMM/HMM.pdf
5. Murphy, K. Hidden Markov Model (HMM) Toolbox for Matlab. (2004) Available from http://www.cs.ubc.ca/~murphyk/Software/HMM/hmm.html
6. Ginno, P. A., Lott, P. L., Christensen, H. C., Korf, I. & Chédin, F. R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters. Mol. Cell 45, 814–825 (2012).
7. Schroeder, D. I., Lott, P., Korf, I. & LaSalle, J. M. Large-scale methylation domains mark a functional subset of neuronally expressed genes. Genome Res 21, 1583–1591 (2011).
8. Schroeder D.I., Blair J.D., Lott P., Yu H.O., Hong D., Crary F., Ashwood P., Walker C., Korf I., Robinson W.P., Lasalle J.M. The human placenta methylome. Proc Natl Acad Sci 15, 6037-42 (2013).