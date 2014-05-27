StochHMM implements 3 traditional way of decoding (Viterbi, N-Best Viterbi, Posterior)

##Viterbi
Viterbi decoding produces the most probable path of states through the complete sequence (Overall Best path).

```
//-------- Import Model ----------//
StateFuncs default_functions;  //Create a StateFuncs with default functions 

//Create model
model hmm;

//Import model
hmm.import(filename,&default_functions);


//------- Import Sequence --------//
//seqTracks handles importing all the sequences for a model
seqTracks jobs;
jobs.loadSeqs(hmm,seq_filename,FASTA);
seqJob *job=jobs.getJob();


//------- Decoding -------//
//Create trellis and initialize with model and sequence
trellis trell(hmm,job->getSeqs());

//Perform viterbi
trell.viterbi();

//Create traceback path
traceback_path path(hmm)

//Get traceback
trell.traceback(path);

//Output traceback as State index
path.print_path();

//Output traceback as State label
path.print_label();

//Output traceback as GFF
path.print_gff();
```

##N-Best Viterbi
N-best Viterbi allows the user to get the nth best path.  However, it is very memory intensive, because you have to store (n) viterbi scores and the traceback for the score for each state. 

##Posterior
Posterior decoding produces the traceback that is the most likely state at every position of the sequence.  
Where the Viterbi is the best overall, the posterior is the best at each position.   It may produce a path that is not valid according the topology of the model.




#Stochastic Sampling
##Stochastic Forward


##Stochastic Viterbi


##Stochastic Posterior