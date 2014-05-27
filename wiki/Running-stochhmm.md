StochHMM implements standard HMM, HMM with duration (for viterbi algorithm) and in the future
will implement hidden semi-Markov model architectures and algorithms. It
grants researchers the power to integrate additional datasets in their
HMM to improve predictions. Finally, it adapts HMM algorithms to provide
stochastic decoding giving researchers the ability to explore and rank
sub-optimal predictions.

To run the StochHMM command-line application you'll need a model file and a sequence.

On command-line type:
```
stochhmm -model <model file> -seq <sequence file> [options]
```

##Command-line options:

###Required command-line options and files
```
-model <model file>	import model file
-seq <sequence file>	import sequence file in fasta format
```

###Non-stochastic Decoding:  Different algorithms available for decoding
```
-viterbi	performs viterbi traceback
-posterior	Calculates posterior probabilities
		If no output options are supplied, this will return the posterior scores
		for all of the states.

	-threshold <score>: Return only the States with a GFF_DESC, if they are
			greater than or equal to the threshold amount.

-nbest <number of paths> 		performs n-best viterbi algorithm
```

###Stochastic Decoding:
```
-stochastic <Type> 
	Types:
		forward		performs stochastic traceback using forward algorithm
		viterbi		performs stochastic traceback using modified-viterbi algorithm
		posterior	performs stochastic traceback using posterior algorithm
-rep  <number of tracebacks to sample>
```

###Output options:
```
-gff			prints path in GFF format
-path			prints state path according to state number
-label			prints state path as labels
-hits			prints hit table for multiple tracebacks for all states at each position
```

##Example Models
A couple example model files have been provided.
```
3_16_Eddy.hmm - GC rich model from Problem 3.16 of 
	"Problems and Solutions in Biological Sequence Analysis. M. Borodovsky and S. Ekisheva. Cambridge
	 Press, UK (2006)"

Dice.hmm - Dishonest Casino Dice model from pg 65 of 
	"Biological Sequence Analysis: Probabilistic models of proteins and nucleic acids. R Durbin, S.
	 Eddy, A. Krogh, and G. Mitchison. Cambridge Press, UK (1998)"

GC-skew.hmm - SkewR model for predicting R-loop forming regions in the human genome. See 
	"Ginno,P.A. et al. (2012) R-loop formation is a distinctive characteristic of unmethylated human
	 CpG island promoters. Mol. Cell, 45, 814–825."
```
[->Goto Models](Example-Models-and-Sequences)

[Next: Running Examples](Running-Examples)