This section is ordered as you would to create functions using the StochHMM C++ library.

1. Defining User-defined Functions
2. Importing Model
3. Importing Sequence
4. Decoding
5. Outputting Results

***

##1. User-defined Functions
Four categories of user-defined State functions are supported by StochHMM. See [State Functions](State-Functions)

1. Transition Function

2. Emission Function

3. PDF Function

4. Multivariate PDF Function


* Track Function which can be called upon import of sequence and return a std::vector<floats>. (Still in preliminary implementation)

* Model Attribute Function which can be called on a sequence and returns a double.  This is used to determine select from multiple models.(Still in preliminary implementation.)

Transition and Emission Functions are functions that will be called when the emission or transition are evaluated.  They are called in addition to calculating the emission or transition.  Whereas, the PDF functions can replace the emission, such that the user-can define their own continuous/discrete functions.   See [Continuous dice model](Dice-Continuous] for an example of model file where the PDF functions replace the discrete distribution defined in the example Dice model.


###Transition State Function
Transition functions are called when a transition is evaluated. The function is of the form:
``` 
//Transition State Function
//Parameters:
//	const std::string*	= Pointer to undigitized sequence
//	const size_t		= Position in sequence being evautated
//	const std::string*	= Pointer to combined/edited sequence based on tracebacks. See State-Functions section
//	const size_t		= Size of combined/edited sequence
double function_name (const std::string*, const size_t, const std::string*, const size_t);
```

###Emission State Function
Emission State functions are called when an emission is evaluated. The function is of the form:
``` 
//Emission Function
//Parameters:
//	const std::string*	= Pointer to undigitized sequence
//	const size_t		= Position in sequence being evaluated
double function_name (const std::string*, const size_t);
```

###PDF State Function
PDF State functions are called as an emission. The function is of the form:
``` 
//PDF Function
//Parameters:
//	const double		= Emission value
//	const std::vector<double>*= Vector of PDF parameters as defined in model's emission. See Dice-Continuous
double function_name (const double, const std::vector<double>*);
```

###Multivariate PDF State Function
Multivariate PDF State functions are called as an emission. The function is of the form:
``` 
//PDF Function
//Parameters:
//	const std::vector<double>*= Vector of all emission values
//	const std::vector<double>*= Vector of PDF parameters as defined in model's emission. See Dice-Continuous
double function_name (const std::vector<double>*, const std::vector<double>*);
```



##Registering User-defined Functions
Registering the fucntion allows StochHMM to associate the function defined in the model file with the correct function. Names of functions must be unique.
```
//Create state functions class
StateFuncs functions;  //Create a StateFuncs with default functions

//Add a Transition State Function
functions.assignTransitionFunction("My_Transition", *transition_function);

//Add a Emission State Function
functions.assignEmissionFunction("My_Emission", *emission_function);

//Add a PDF State Function
functions.assignPDFFunction("My_PDF", *pdf_function);

//Add a multivariate PDF State Function
functions.assignPDFFunction("My_multi_PDF", *multi_pdf_function);
```

###Example of [Dice-Continuous](Dice-Continuous)
```
//Defined PDF function for fair di emissions. Return probability
double fair(const double& val, const std::vector<double>* param){
	return 0.167;
}

//Loaded di function  Return probability
double loaded(const double& val, const std::vector<double>* param){
	if (val == 6.0){
		return 0.5;
	}
	
	return 0.1;
}


//Define State PDF Functions names- Correspond with those found in the model
functions.assignPDFFunction("FAIR", *fair);
functions.assignPDFFunction("LOADED", *loaded);
```

***

##2. Import Model
```
//Create model
model hmm;

//Import filename of model using user-defined functions--
hmm.import(filename,&default_functions);

//Import filename of model without user-defined functions
hmm.import(filename);

```


##3. Import Sequences
```
//seqTracks handles importing all the sequences for a model
seqTracks jobs;
jobs.loadSeqs(hmm,seq_filename,FASTA);

//For single sequence
seqJob *job=jobs.getJob();
sequences* temp = job->getSeqs()

//For multiple sequences in fasta file
seqJob *job=jobs.getJob();
while (job!=NULL){

	//Perform decoding Here

	job = jobs.getJob() //Get next sequence
}
```

***

##4. Decoding
###Initialize Trellis
```
//Create trellis and initialize with model and sequence
trellis trell(hmm,job->getSeqs());
```

###Perform analysis

####Viterbi
```
//Perform viterbi
trell.viterbi();

//Create traceback path class
traceback_path path(hmm)

//Get decoded traceback path
trell.traceback(path);
```

####Posterior
```
//Perform posterior
trell.posterior();

//Perform Posterior decoding
traceback_path path(hmm);
trell.traceback_posterior(path);

//Get posterior table to process scores
trell.getPosteriorTable();
```

####Nth-best Viterbi
```
//Perform N-best (get N top paths)
trell.naive_nth_viterbi(N);
//Get the N tracebacks and output them
for(size_t i=0;i<nth;i++){
	traceback_path path(hmm);
	trell.traceback_nth(path, i); //ith path
	//print_output here
}
```

####Stochastic Viterbi
```
//Perform stochastic viterbi
trell.stochastic_viterbi();
multiTraceback paths;
trell.stochastic_traceback(paths, repetitions);

```

####Stochastic Forward
```
//Perform stochastic forward
trell.stochastic_forward();
multiTraceback paths;
trell.stochastic_traceback(paths, repetitions);

```

***

## 5. Outputting Results
###Printing traceback_path class
```
//Output traceback as State index
path.print_path();

//Output traceback as State label
path.print_label();

//Output traceback as GFF
path.print_gff();
```