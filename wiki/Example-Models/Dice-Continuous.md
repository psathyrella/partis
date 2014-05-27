##Occasionally Dishonest Casino Model

[Durbin, R., Eddy, S., Krogh, A., Mitchison, G., Biological Sequence Analysis: Probabilistic models of proteins and nucleic acids. Cambride University Press, 1998, pg 54]

This model uses two PDF emission functions which are defined at compile time. 

## PDF Function code for the LOADED and FAIR states

```
//Defined PDF function for fair dice emissions. Return log of probability
// 0.167 for everything.   Could be expanded to check for value and 
// return a different value if we see a number outside the range of [1,6]
double fair(const double val, const std::vector<double>* param){
	return log(0.167);
}

//Loaded dice function  Return probability of log(0.5), if dice roll is 6.
//Returns log(0.1) otherwise.   
double loaded(const double val, const std::vector<double>* param){
	if (val == 6.0){
		return log(0.5);
	}
	
	return log(0.1);
}

int main(int argc, const char * argv[])
{	
	srand(time(NULL));

	model hmm; // HMM model 
	StateFuncs functions;  //State functions


	//Define State PDF Functions names- Correspond with those found in the model
	//Assign pointer to functions we defined above;
	functions.assignPDFFunction("FAIR", *fair);
	functions.assignPDFFunction("LOADED", *loaded);

//Check that model argument is defined and import the model	
	std::string model_file = "Dice_Continuous.hmm";
	std::string seq_file   = "Dice_real.fa";

	//Import model (Pass State functions to model import)
	hmm.import(model_file,&functions);

	//Import sequences
	jobs.loadSeqs(hmm, seq_file, FASTA);
    
	seqJob *job=jobs.getJob();
    
	//Create Trellis
	trellis test_trellis(job->getModel(),job->getSeqs());
    
	//Perform Viterbi
	test_trellis.viterbi();
	
	//Get traceback and print the labels
	traceback_path test(job->getModel());
	test_trellis.traceback(test);
	test.print_label();
	
	return 0;
}
```

##MODEL FILE

**Note**: The parameters in this model are not used and are only supplied to illustrate how they are passed to the function.  The fair function will receive a vector composed of [1,2,3,4,5].  In a model, the parameters could be used to describe the parameters of a probability distribution function, the index to a table to use for emissions, or any other possibility.


```
#STOCHHMM MODEL FILE
<MODEL INFORMATION>
======================================================
MODEL_NAME:	CASINO DICE MODEL
MODEL_DESCRIPTION:	Taken from CH3 Durbin/Eddy
MODEL_CREATION_DATE:	August 28,2009

<TRACK SYMBOL DEFINITIONS>
======================================================
DICE:	REAL_NUMBER

<STATE DEFINITIONS>
#############################################
STATE:	
	NAME:	INIT
TRANSITION:	STANDARD: P(X)
	FAIR:	0.5
	LOADED:	0.5
#############################################
STATE:	
	NAME:	FAIR
	PATH_LABEL:	F
	GFF_DESC:	FAIR
TRANSITION:	STANDARD: P(X)
	FAIR:	0.95
	LOADED:	0.05
	END:	1
EMISSION:	DICE: CONTINUOUS
	PDF: FAIR	PARAMETERS: 1,2,3,4,5
#############################################
STATE:
	NAME:	LOADED
	PATH_LABEL:	L
	GFF_DESC:	LOADED
TRANSITION:	STANDARD: P(X)
	FAIR:	0.1
	LOADED:	0.9
	END:	1
EMISSION:	DICE: CONTINUOUS
	PDF: LOADED	PARAMETERS: 6,7,8,9,10
#############################################
//END
```


##Sequence File

```
>Dice
3,1,5,1,1,6,2,4,6,4,4,6,6,4,4,2,4,5,3,1,1,3,2,1,6,3,1,1,6,4,1,5,2,1,3,3,6,2,5,1,4,4,5,4,3,6,3,1,6,5,6,6,2,6,5,6,6,6,6,6,
6,5,1,1,6,6,4,5,3,1,3,2,6,5,1,2,4,5,6,3,6,6,6,4,6,3,1,6,3,6,6,6,3,1,6,2,3,2,6,4,5,5,2,3,6,2,6,6,6,6,6,6,2,5,1,5,1,6,3,1,
2,2,2,5,5,5,4,4,1,6,6,6,5,6,6,5,6,3,5,6,4,3,2,4,3,6,4,1,3,1,5,1,3,4,6,5,1,4,6,3,5,3,4,1,1,1,2,6,4,1,4,6,2,6,2,5,3,3,5,6,
3,6,6,1,6,3,6,6,6,4,6,6,2,3,2,5,3,4,4,1,3,6,6,1,6,6,1,1,6,3,2,5,2,5,6,2,4,6,2,2,5,5,2,6,5,2,5,2,2,6,6,4,3,5,3,5,3,3,3,6,
2,3,3,1,2,1,6,2,5,3,6,4,4,1,4,4,3,2,3,3,5,1,6,3,2,4,3,6,3,3,6,6,5,5,6,2,4,6,6,6,6,2,6,3,2,6,6,6,6,1,2,3,5,5,2,4,5,2,4,2
```