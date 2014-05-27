Scores are reported as natural log values. (Posterior and Stochastic scores not reported)
***


##Dice Model Examples
Each example provides the command used and the output from StochHMM

####Print PATH_LABEL using Viterbi algorithm
```
#Print Viterbi traceback as State Path Label
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -viterbi -label
>>Eddy Dice TRACK_NAME:TRACK1	Score: -539.062
F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F 
L L L L L L L L L L L L L L L L L L F F F F F F F F F F F F L L L L L L L L L L L L L L L L L L 
L L L L L L L L L L L L L L L L F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F 
F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F L L L L L L L L L L L L L 
F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F 
F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F L L L L L L L L L L L L L L L L L L 
L F F F F F F F F F F F
```

####Print GFF using Viterbi algorithm
```
#Print Viterbi traceback as GFF (Only states with GFF_DESC will be output)
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -viterbi -gff
#Score: -539.062
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	1	48	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	49	66	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	67	78	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	79	112	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	113	179	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	180	192	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	193	270	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	271	289	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	290	300	.	+	.
```

####Print state number using Viterbi algorithm
```
#Print Viterbi traceback as state position in model
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -viterbi -path
>>Eddy Dice TRACK_NAME:TRACK1	Score: -539.062
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 0 0 0 0 0 0 0 0 0 0 0 
```

####Print Posterior scores
```
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -posterior

Posterior Probabilities Table
Model:	CASINO DICE MODEL
Sequence:	>Eddy Dice TRACK_NAME:TRACK1
Probability of Sequence from Forward: Natural Log'd	-516.544812
Probability of Sequence from Backward:Natural Log'd	-516.544812
Position	FAIR	LOADED
1	0.812	0.188
2	0.849	0.151
3	0.861	0.139
4	0.850	0.150
5	0.814	0.186
6	0.738	0.262
...(output truncated)
```

####Print Posterior decoding as GFF
```
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -posterior -gff
#Score: 0
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	1	47	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	48	66	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	67	78	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	79	95	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	96	104	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	105	112	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	113	129	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	130	138	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	139	179	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	180	192	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	193	201	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	202	207	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	208	269	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	270	289	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	290	300	.	+	.
```


####Print GFF using Stochastic Viterbi algorithm 10 samples
Reports the number of times the specific traceback path has occurred during the sampling
```
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -stochastic viterbi -rep 10 -gff
Traceback occurred:	 1
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	1	50	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	51	68	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	69	89	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	90	112	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	113	131	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	132	143	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	144	154	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	155	155	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	156	183	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	184	196	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	197	203	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	204	210	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	211	275	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	276	288	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	289	300	.	+	.

Traceback occurred:	 1
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	1	48	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	49	63	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	64	65	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	66	66	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	67	81	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	82	85	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	86	86	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	87	117	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	118	131	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	132	142	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	143	167	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	168	168	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	169	179	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	180	209	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	210	230	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	231	231	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	232	270	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	LOADED	271	288	.	+	.
Eddy Dice TRACK_NAME:TRACK1	StochHMM	FAIR	289	300	.	+	.

...(output truncated)
```

####Print Hits table for Stochastic Viterbi sampling
```
$ stochhmm -model ../examples/Dice.hmm -seq ../examples/Dice.fa -stochastic viterbi -rep 10 -hits
Position	FAIR	LOADED
1	10	0
2	10	0
3	10	0
4	10	0
5	10	0
6	8	2
7	6	4
... (truncated output)
```
