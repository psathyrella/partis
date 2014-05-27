#State Functions
State functions are user-defined functions that are called within a state to get a probability to apply to the scoring.   The functions can be called at either evaluation of an emission or transition.

The function passed to assigned in a StateFunc class object and passed to the model import function.   Only a single StateFunc can be passed to the model import.  If multiple functions need to be passed to the model, the functions need to be assigned to a single StateFunc object.

(Note: State function values are returned in log space.   If a tag function or a function is part of a multiple emission, it will be combined with the corresponding type by multiplication (addition in log space).)


###Creating and Assigning StateFunc

```
StateFunc stFunc;

//Assigning a transition function
stFunc.assignTransitionFunction( const char* <Name of Function>, <pointer to transition function>);

//Assigning a emission function
stFunc.assignTransitionFunction( const char* <Name of Function>, <pointer to emission function>);

//Assigning a univariate probability density function
stFunc.assignPDFFunction( const char* <Name of Function>, <pointer to probability density function>);

//Assigning a multivariate probability density function
stFunc.assignMultivariatePdfFunction( const char* <Name of Function>, <pointer to multivariate probability density function>);

```

##Emission Functions
Emission functions can be defined and called on any State emission.  The return value of this function will be applied to the standard emission probability as a natural logarithmic value.  This allows the user to either use a function to calculate the emission probability or use a function to weight the emissions based on some other criteria. 

###Definition
Emission functions need to return a double and take the parameters (const std::string*, size_t).
The std::string the function will receive is the sequence for the defined track.   The size_t will be the position within the string being evaluated by the emission.

Any function of the form ```double myFunctionName (const std::string*, const size_t)``` can be used for an emission function.

###Assigning Function to StateFunc object

```
StateFunc stFunc;

//Assigning a emission function
stFunc.assignEmissionFunction( "MyFunction", &myFunctionName);
```

###Calling within the Model
Within the model the function will be called at any emission where Emission Tag is placed.

####Example
```
EMISSION:	SEQ: P(X)	[FUNCTION: MyFunction	TRACK: SEQ	SCALE: 0]
	ORDER:	0
@A	C	G	T	
0.1	0.1	0.1	0.1
```

When evaluating this emission, the function "MyFunction" will be called and passed the track labeled "SEQ".   The returned score will be scaled by a value 0 or any be the corresponding SCALING Function that was defined in the model.

Note:  The ```SCALE``` is optional in Emission TAG.  If ```SCALE``` is left out the value will be passed, directly to the model and not scaled by any value, which is equivalent to ```SCALE: 0```


##Transition Functions
A transition function is called anytime the corresponding transition is evaluated. The return value of this function will be applied to the standard transition probability as a natural logarithmic value.  When the transition is evaluated the function can traceback through the trellis and pass the corresponding sequence to the function for evaluation.

###Definition
The transition function returns a double and must have the form ```double MyFunctionName (const std::string*, const size_t, const std::string*, const size_t)```.   The first string corresponds to the complete sequence.  The first size_t corresponds to the current position within the sequence that is being evaluated.  The second std::string is determined by a traceback and combining paramenters set in the tag.  The second size_t corresponds to the traceback lenght of the traceback.


###Transition TAG

FUNCTION:  Function Name to use

TRACK:  Name of Track to use

SCALE: (optional) Scaling parameter to use to scale the score

TRACEBACK LABELS:  see valid Traceback Labels below

COMBINE TAG: see valid Combine tags below

####Traceback Labels:  (Alternate TAG separated by "/")

TO_LABEL / TB->LABEL: Traceback through trellis until the specified State PATH_LABEL is reached

TO_GFF / TB->GFF: Traceback through trellis until the specified State GFF_DESC is reached

TO_STATE / TB->STATE: Traceback through trellis until the specified State NAME is reached

TO_START / TB->START: Traceback through trellis until the start is reached

DIFF_STATE: Traceback through trellis until a different State is reached.


####Combine Tags:  Combine sequence according to specified criteria

COMBINE_LABEL:  Combine the sequence where the traceback State's PATH_LABEL is equal to specified value

COMBINE_GFF: Combine the sequence where the traceback State's GFF_DESC is equal to specified value

COMBINE_STATE: Combine the sequence where the traceback State's NAME is equal to specified value

NO_COMBINE:  Don't combine the sequence in any manner.   Will return the sequence corresponding to the traceback.

###Assigning Function to StateFunc object

```
StateFunc stFunc;

//Assigning a transition function
stFunc.assignTransitionFunction( "MyFunction", &myFunctionName);
```


###Calling within the Model
Within the model the function will be called at any transition where Transition Tag is placed.  Tags are surrounded by square brackets [...]

####Example

```
TRANSITION:	STANDARD:	P(X)
ACCEPTOR2: 1	[FUNCTION:	HMMER	TRACK:	SEQ	COMBINE_LABEL:	E	TO_LABEL:	N	SCALE:	0.4	]
```

Anytime the ACCEPTOR2 transition was evaluated, it would also traceback to any State with a ```PATH_LABEL``` of "N". Then the sequence from track SEQ which corresponded to any States with the ```PATH_LABEL``` "E" would be combined. The function HMMER would then be called and be passed the complete SEQ sequence, the current position within the trellis being evaluated, the combined sequence and the total length of the traceback.  
The returned value of the function would be scaled and that value would be applied(added in log space) to the set probability for transition to ACCEPTOR2 state of P(X) = 1; 


####Computation Costs
With the traditional HMM the emissions and transitions are based on a simple table lookup.  Within StochHMM the Lexical emissions and transitions are retrieved from a table.   When including a user-defined function the computation complexity can greatly increase depending upon you function.

The basic lexical emissions retrieval routine calls a *lexical* class to retrieve the associated value from a table based upon some emission (0th or higher order).  In the case of Lexical emission, the emission is a character or word, with it's associated dependencies. Upon calculating the proper indices for the work and its dependencies, the *Lexical* class can retrieve them form the associated tableUpon retrieving the value for the emission the *lexical class* returns the value to *emission* class, which promptly combines this with any other emissions (by multiplication) (addition in log space.)   The same emission function calls the user-defined functions, which are expected to return a value in log space.

Based upon the call to user-defined function, it can be expected that the complexity of the call could be less than the baseline or greater.   If for example, the user defined function is using the position and a table of emission values, it could be reduce the amount of computation below that used by the *Lexical* class.   However, if the user-defined function calls a function that has an larger computational complexity  the additional cost will be accrued at each call to the function.

To theoretically calculate the savings or cost of your user-defined function for a worst-case scenario, you'll need to understand that the HMM could call your function approximately (length of sequence * number of states) times.

Cost = (length of sequence * states * user-defined function cost)


To minimize calls to the emissions and transition functions, StochHMM has implemented some low-cost optimizations into the algorithms.  During the calculation of each position, StochHMM keeps track of what states are possible in the next step and which are not.   StochHMM only calls the function for those states which have a possible transition from previous state with a non-zero score. This does come at a slight cost and we have made every effort to minimize the computational cost.   

For example, in the Occasionally Dishones Casino model, StochHMM will still be keeping track of what states are possible.  Because there are only 2 states and they are completely connected, there isn't any improvement in speed using this optimization.

However, when working with a less connected model such as a 200-state gene model, where some paths are impossible for example (1st Codon to 3rd Codon), we can use the possible transitions to inform and limit our calls to the emissions.    Also, not all paths are possible because of emissions, for example if we've seen a *T* in the previous state, we don't need to visit the a state that describes the second position in the Start codon.  (Case where accepted start codon is ATG).
This allows us to reduce the number of times that emissions or user-defined functions are called in less connected models.

If your user-defined function comes at a high cost, you possibly address the cost by making sure that the function is called infrequently (this could be done by placement in the model).  Taking advantage of the fact that StochHMM won't evaluate the state unless it's predecessor is valid, the function would only be evaluated when it was necessary.

If the function is costly and the model is such that it will be evaluated frequently, then it may be advantageous to precompute the values and integrate them using a simple user-defined function or emission/transitions.  



   