#StochHMM
##A flexible Hidden Markov Model Framework in C++

**News:** Check out how Matthew Porter used StochHMM to plug in custom functions to state emissions in order to  quickly create a HMM to annotate copy number from SAM files. [Ploidamatic](https://github.com/KorfLab/Ploidamatic)  *Manuscript submitted for review at Oxford Bioinformatics.*

1. [Introduction](wiki/Introduction)
2. **NEW** [Comparison to Mamot / R-HMM / HMMoc](wiki/comparison)
3. [Installing/Compiling StochHMM](wiki/Install)
4. [Running StochHMM](https://github.com/KorfLab/StochHMM/wiki/Running-StochHMM)
	1. [Running Examples](wiki/Running-Examples)
5. [Developing a HMM using StochHMM](wiki/Developing-HMMs-using-StochHMM)
6. File Formats
	1. [Model File](wiki/Model-File)
		1. [Model Information](wiki/Model-Information)
		2. [Track Symbol Definitions](wiki/Track-Symbol-Definitions)
			1. [Track Functions](wiki/Track-Functions)
		3. [Ambiguous Symbol Definitions](wiki/Ambiguous-Symbol-Definitions)
		4. [Scaling Values](wiki/Scaling-Functions)
		5. [Templated States](wiki/Templated_States)
		6. [State Definitions](wiki/State-Definitions)
			1. [Transitions](wiki/State-Transitions)
			2. [Emissions](wiki/State-Emissions)
			3. [Additional User-defined State Functions](wiki/State-Functions)
				1. [Pre-defined Univariate and Multivariate PDF Functions](wiki/Predefined-PDFs)
	2. [Sequence Files](wiki/Sequence-Files)
		1. [External Definitions](wiki/External-Definitions)
	3. [Models Files](wiki/Multiple-Model-File)
7. [Examples Model and Sequences](wiki/Example-Models-and-Sequences)
	1. [Simple CpG](wiki/CpG)
	2. [CpG with Ambiguous](wiki/Cpg_ambiguous)
	3. [Simple Dice](wiki/Dice)
	4. [User-Defined Emissions Dishonest Casino](wiki/Dice-Continuous)
	5. [Multiple Independent Emissions](wiki/Multi-Independent)
	6. [Multiple Joint Emissions](wiki/Multi-Joint)

8. [C++ Library](wiki/C---Library)
	1. [Example Code](wiki/Example-code)
	2. [StochHMM main explained](wiki/stochhmm-main)

9. [Additional Resources](wiki/Additional-Resources)