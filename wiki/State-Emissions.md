#Emissions
A state can emit a single value or multiple values at each time point(see Multi-Emission States below).    Emissions are defined by the type of track that the emission uses.  Emissions can be either Lexical, Lexical Function, Real_Number, or Continuous.

Emissions must be defined after the Transition has been defined.

The emissions under an HMM will output a value between 0 and 1.   However, emission values will not be checked.

Multiple emissions will be combined in log space via addition in the emission subroutines.  (Note: Joint emissions are stored and retrieved as a single value for the emissions)  

##Lexical Emission


Lexical emissions provide emissions for each character P(X=Char|Word).  They can have a higher order.  Higher orders increase the probability table by the power of ORDER+1.  For 4 characters in 5th order, the emissions table is 4^6 = 4096 cells.

Emission values can be provided as counts, probabilities P(X), or natural log(P(x)).

To specify the VALUE_TYPE include:
	COUNTS for counts;
	P(X) for probabilities
	LOG for natural log probabilities
	
When COUNTS are defined, rows will be converted to frequencies.
P(X) and LOG values are not checked to see if the row total value is equal to 1.

By providing P(X) or LOG values, the user can set probabilities that are more inline with a CRF rather than an HMM.
	
###Formatting:
```
EMISSION:	TRACK_NAME:	VALUE_TYPE
	ORDER:	ORDER_NUMBER
VALUES....
```

###Example 
```
EMISSION:	SEQ:	COUNTS
	ORDER:	0
1	0	0	0
```

The emission uses the SEQ track and the values provided are counts.   The order of dependence is set to zero.   The columns of the table following the emission should correlate with the order of the TRACK SYMBOLS.   Column 1 = A, Column 2 = C ....

For higher order emissions the rows correlate with words.

```
EMISSION:	SEQ:	COUNTS
	ORDER:	2
1	0	0	0	#Row AA
1	0	0	0	#Row AC
1	0	0	0	#Row AG
1	0	0	0	#Row AT
1	0	0	0	#Row CA
1	0	0	0	#Row CC
1	0	0	0	#Row CG
	...
1	0	0	0	#Row TT
```

Lexical tables can also have the word for easier reference, but the order of words must be maintained.  (They aren’t checked).   By default, StochHMM will reprint models in this format.

Example of including row word definitions:

```
EMISSION:	SEQ:	COUNTS
	ORDER:	2
@A	C	G	T
@AA	1	0	0	0	#Row AA
@AC	1	0	0	0	#Row AC
@AG	1	0	0	0	#Row AG
@AT	1	0	0	0	#Row AT
@CA	1	0	0	0	#Row CA
@CC	1	0	0	0	#Row CC
@CG	1	0	0	0	#Row CG
	...
@TT	1	0	0	0	#Row TT
```

###Example Models

##Lexical Ambiguity

Lexical Emissions also support ambiguity.  For example in nucleotide sequence, N could be either A,C,G,T/U.   If ambiguity is not defined for an emission, then all ambiguous characters will be set to P(X)=0.

Ambiguous characters can be scored by either the emission: MIN, MAX, AVG, or assigned a constant value;

###Format:

```
EMISSION:	TRACK_NAME:	VALUE_TYPE
	ORDER:	ORDER_NUMBER	AMBIGUOUS: [SCORING_PARAMETERS]
VALUES....
```


###Example:

```
EMISSION:	SEQ:	COUNTS
	ORDER:	5	AMBIGUOUS:	P(X): 0.0002
...
```

OR 

```
EMISSION:	SEQ:	COUNTS
	ORDER:	5	AMBIGUOUS:	LOG: -INF
```

OR 

```
EMISSION:	SEQ:	COUNTS
	ORDER:	5	AMBIGUOUS:	AVG
```

OR 

```
EMISSION:	SEQ:	COUNTS
	ORDER:	5	AMBIGUOUS:	MIN
```

OR 

```
EMISSION:	SEQ:	COUNTS
	ORDER:	5	AMBIGUOUS:	MAX
```

###Example Models


##Lexical Function Emissions
A user can define a function to compute the emission probability and link this to the state emission.  The PWM_START function will be called to calculate the natural log probability.

The PWM_START function will be passed a std::string of the sequence for the defined Track and a size_t, which is the position currently being evaluated.

###Format

```
EMISSION:	TRACK_NAME:		FUNCTION:	FUNCTION_NAME
```

###Example

```
EMISSION:	SEQ:		FUNCTION:	PWM_START
```


###Example Models


##REAL Number Emissions
A REAL number emissions is a given real number value that the state emits, corresponding to a REAL_NUMBER track.  Values are considered to be natural log based when imported or created by Track functions. 

###Format

```
EMISSION:	TRACK_NAME:		REAL_NUMBER
```

###Example

```
EMISSION:	PWM:		REAL_NUMBER
```

Where PWM is previously defined as a REAL_NUMBER track.


###REAL Number Track Complement
The complement (1-P(X)) of a REAL_Number track can be used in the emission by specifiying 1-P(X) or COMPLEMENT

###Example

```
EMISSION:	SIDD:		REAL_NUMBER	COMPLEMENT
```

###Example Models


##Continuous PDF Emission

A PDF emission uses a real number track and output probabilities based on a function.   That function could be a probability distribution function or any user-defined function. 

StochHMM has a few pre-defined probability distribution functions. For a list of pre-defined distributions, ```PDF_FUNCTION_NAME``` and the required parameters, see [Probability Distribution Functions](Predefined-PDFs).  The order of parameters must match the PDF function order. If more parameters than necessary are defined, these will be ignored.


###Format
```
EMISSION:	TRACK_NAME: CONTINUOUS
	PDF: PDF_FUNCTION_NAME	PARAMETERS: PARAMETERS COMMA SEPARATED
```

###Example
```
EMISSION:	DICE: CONTINUOUS
	PDF: FAIR	PARAMETERS: 1,2,3,4,5
```

###Example Model
See [Dice_Continuous.hmm]()

###Compliment of PDF
####Example
```
EMISSION:	DICE: CONTINUOUS	COMPLEMENT
	PDF: FAIR	PARAMETERS: 1,2,3,4,5
```
OR

```
EMISSION:	DICE: CONTINUOUS	1-P(X)
	PDF: FAIR	PARAMETERS: 1,2,3,4,5
```




***


##Multi-emissions states

A state can emit multiple symbols (individually or jointly).  To define a multiple emission state in which the emissions are individually output, all that needs to be done is to include multiple emission definitions.

###Independent Emissions
Each emission probability is considered independent of the other.  So the final emission for the state will be the product of all emission probabilities.

###Example


```
EMISSION:	SEQ:	COUNTS
	ORDER:	0
1	0	0	0

EMISSION:	USER:	COUNTS
	ORDER:	0
1	0	0	0

EMISSION:	SIDD:		REAL_NUMBER
```

Here the state outputs three emissions.


##Joint Emission
Joint emission probabilities are used when the joint probabilities are not independent.  A Joint emission can be defined using 2 or more different tracks. Each track may have the same or different ORDER.

###Example
To create a joint emission for SEQ and USER tracks.   
SEQ track has characters A,C,G,T
USER track has characters 1,2,3,4
Then the joint emission probability for a state outputting A and 2 would be P(A and 2) 

```
EMISSION:	SEQ,USER:	COUNTS
	ORDER:	0,0
@A1	A2	A3	A4	C1	C2	C3	C4	G1	G2	G3	G4	T1	T2	T3	T4
1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
```

The order for both emissions is set to 0 (SEQ track) , 0 (USER track).  The resulting table would be 1 row x 16 columns.   The columns correspond to the joint word.

The order of emissions can also differ for each track used.  For example, we may want a first order emission for SEQ and a zero order emission for USER

 
```
EMISSION:	SEQ,USER:	COUNTS
	ORDER:	1,0
@A1	A2	A3	A4	C1	C2	C3	C4	G1	G2	G3	G4	T1	T2	T3	T4
@A	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
@C	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
@G	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
@T	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
```

##Multivariate Continuous PDF Emission
Multivariate Continuous PDF functions are defined similar to Univariate Continuous PDF emissions. Currently, there are no multivariate PDF functions loaded by default into StochHMM.   In order to use Multivariate Continuous PDF, the function need to be created at compile-time and the function needs to be added to the StateFuncs index, which is then passed to HMM import function.  (See StateFuncs section)


###Format
```
EMISSION:	<TRACK_NAMES (COMMA-SEPARATED>: MULTI_CONTINUOUS
	PDF: <PDF_FUNCTION_NAME>	PARAMETERS: <PARAMETERS (COMMA-SEPARATED)>
```

####Example
```
EMISSION:	TRACK1,TRACK2: MULTI_CONTINUOUS
	PDF: NORMAL	PARAMETERS: 1,2,0.1,0.2,0.1,0.2
```

#Emission State Function
An Emission State Function is simply an additional function that can be linked with an emission.  When ever the emission is checked or calculated the fucntion will be called and return a value.

Emission State Function takes a pointer to a std::string to a track sequence, and a size_t value for current position.   This function could be used for linking the emission to a position weight matrix score or some upstream analysis that couldn’t be performed before hand.   In this way the HMM emission can be link to a long-range scanning function.

The Emission State Function should be placed on the same line as the emission definitioin and is surrounded by square brackets.

The functions output can be scaled by a constant value or by <SCALING VALUE> definition.

###Format
```
EMISSION DEFINITION	[ FUNCTION:	FUNCTION_NAME	TRACK: TRACK_NAME	SCALE: VALUE OR SCALING FUNCTION]
```

###Example
```
EMISSION:	SEQ:	COUNTS	[FUNCTION: TestPWM		Track: SEQ	Scale: 0.5]
```