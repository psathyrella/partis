##StochHMM - Stochastic Hidden Markov Model Framework

##Usage: StochHMM -model <model_file>  -seq <seq_file> [options]
Example: ./StochHMM -model Dice.hmm -seq Dice.fa -viterbi -gff

##Command Line options:

###Files: StochHMM requires a sequence file and a model file
	-model <model file>		import model file
	-seq <sequence file>		import sequence file in fasta format

###Non-stochastic Decoding:  Different algorithms available for decoding
	-viterbi		performs viterbi traceback
	-nbest <int>		performs the N-best Viterbi and return # number of paths

###Posterior Decoding:
	-posterior		Calculates posterior probabilities

###Stochastic Decoding
	-stochastic <type> : Performs stochastic decoding using algorithm type. Must also set number of  repetitions
		Types:
			viterbi		 Uses viterbi scores to calculate stochastic traceback probabilities
			forward		 Uses forward scores to calculate stochastic traceback probabilities
			posterior	 Uses posterior scores to calculate stochastic traceback probabilities

	-repetitions <int>	Number of stochastic tracebacks to perform

###Output options:
	-gff		prints path in GFF format
	-path		prints state path according to state number
	-label		prints state path as labels
	-trellis	print the trellis values

Written by Paul Lott at University of California, Davis
Please direct any questions, suggestions or bugs reports to Paul Lott at plott@ucdavis.edu

[Next: Running Examples](Examples)