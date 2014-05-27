#State Transitions
A state must define at least one transition for each state.  Transitions can be one of three 3 types (STANDARD, LEXICAL, DURATION).  STANDARD transitions are a simple probability associated with the states transition from one state to the next.  LEXICAL transitions are transitions that use the sequence to figure out the probability (Mealy Machine).  DURATION transitions define a survival distribution for the length that a state is likely to emit.

##GENERAL FORMAT
```
TRANSITION:	TRANSITION_TYPE:	VALUE_TYPE
	STATE_NAME:	VALUE OR TRACK_NAME OR TB_INFO
```

VALUE_TYPE:  P(X) or LOG (natural log)

STATE_NAME: Name of state that transition describes


##STANDARD
Standard transitions have a fixed probability associated with the transition.  We define standard transition as shown below.

###Format
```
TRANSITION:	STANDARD:		VALUE_TYPE
	STATE_NAME:	VALUE
```

###EXAMPLE
```
TRANSITION:	STANDARD:		P(X)
	INTER:	1.0
```

or

```
TRANSITION:	STANDARD:		LOG
	INTER:	0.0
```

Describes a transition to INTER state with a probability 1.0.

We can also include multiple transitions, under a single STANDARD transition definition.

```
TRANSITION:	STANDARD:		P(X)
	INTER:	0.5
	START:	0.5
```

If it is valid to end the sequence with the state then a STANDARD transition to END state must be defined.

```
TRANSITION:	STANDARD:		P(X)
	INTER:	0.5
	START:	0.5
	END:	1.0
```
 
 
##LEXICAL
Lexical transition use the sequence to determine the likelihood of transitioning from one state to the next.   So, like the emissions, we must provide sequence counts, frequencies P(X), or ln(P(X)) for emitted observations.   Because they are sequence dependent, we’ll have a table corresponding to a sequence track and track alphabet.

###FORMAT
```
TRANSITION:	LEXICAL:	VALUE_TYPE
	STATE_NAME:	TRACK_NAME
		ORDER:	ORDER_VALUE
SEQUENCE_P(X)_TABLE
```

Rows and columns are defined the same as in the emissions.  If we have an alphabet for SEQ of "A,C,G,T" then Column 1 corresponds to "A", Column 2="C", Column 3="G", and Column 4="T".   Row 1 corresponds to "given A", Row 2= "given C", Row 3="given G", Row 4="given T".

###EXAMPLE
```
TRANSITION:	LEXICAL:		P(X)
	INTER:	SEQ
		ORDER:	1
0.996096	0.998151	0.998046	0.997835
0.997860	0.997702	0.998058	0.997916
0.997821	0.997921	0.997690	0.997879
0.997816	0.997893	0.997766	0.996080
```

This example defines the transition to a state called "INTER" as 1st-order dependence upon the sequence.
If the next state is emitting an A and the current state is emitting a T, then the transition probability P(A|T) = 0.997816.


OR 

```
TRANSITION:	LEXICAL:		P(X)
	INTER:	SEQ
		ORDER:	1
@A	C	G	T
@A	0.996096	0.998151	0.998046	0.997835
@C	0.997860	0.997702	0.998058	0.997916
@G	0.997821	0.997921	0.997690	0.997879
@T	0.997816	0.997893	0.997766	0.996080
```

The transition is lexical with values provided as probabilities.   The transition is to the INTER state and we’ll use the SEQ track to define the transition probability.   The order of dependence is first order, meaning we’ll look at the probability P(X | Y).

If our sequence was ACGT and we are currently at C position in our sequence, we’d compute P(C | A) = 0.998046


##Ambiguous Characters
If we want to define how the ambiguous characters are treated, we need to type AMBIGUOUS

Scoring of ambiguous can be AVG, MAX, MIN, P(X) or LOG  See [Emissions](State-Emissions)

###EXAMPLE
```
TRANSITION:	LEXICAL:		P(X)
	INTER:	SEQ
		ORDER:	1	AMBIGUOUS: AVG
0.996096	0.998151	0.998046	0.997835
0.997860	0.997702	0.998058	0.997916
0.997821	0.997921	0.997690	0.997879
0.997816	0.997893	0.997766	0.996080
```

##LEXICAL FUNCTION (Compile Time)
**Requires definition of Lexical Function at Compile time.
 
Lexical transition can also link to an emission function, which uses a lexical table and sequence to calculate the transition probability.  See [Lexical Function Emission](State-Emissions)   

###EXAMPLE
```
TRANSITION:	LEXICAL:		FUNCTION:  PWM2
	INTER:	SEQ
```

This defines a Transition to the state INTER using the SEQ track.   The transition probability is calculated using the PWM2 function



##DURATION Transition (Preliminary)
Duration transitions are dependent upon how long has been spent in the state.  The initial length will be calculated by tracing back until a given condition. That length will be used to calculate the transition probability.

The duration transition code is still very preliminary.   The forward/backward algorithm have not yet been implemented for duration transitions.

###Traceback Options 
DIFF_STATE: Traceback until a different state is encountered
TO_STATE:  Traceback until a state with given name is encountered
TO_LABEL: Traceback until a state with given label is encountered
TO_GFF: Traceback until a state with a given GFF description is encountered
TO_START: Traceback to the start of the sequence


###FORMAT
```   
TRANSITION:	DURATION:	VALUE_TYPE
	STATE_NAME:	TRACEBACK_OPTION
		LENGTH	VALUE
		LENGTH	VALUE
		...
```


```   
TRANSITION:	DURATION:	P(X)
	ENTER:	DIFF_STATE
		3	1
		4	0.995
		5	0.50
		10	0.10
		100	0.01
		101	0
```
When initially encountered, the transition to ENTER will traceback until the label of states is different than the current state.  The length will then be used to calculate the transition probability.

This definition would define the transition to the ENTER state.  Anything less than 4 would be considered probability 1.0. And anything more than 100 would have a probability of zero. 

If the position index of the distribution begins with a number that number will be extended to the extremes.  For example using the above transition, 2 would be 1 also.  

If position 101, ended with 0.1 then 102-> ∞ would also be 0.1

If you want the probability of less than 3 to be zero, then you should define 2 as P(X) = 0.

To describe the outgoing transition we’d need to define a reciprocal transition

```
TRANSITION:	DURATION:	P(X)
	ENTER:	DIFF_STATE
		3	1
		4	0.995
		5	0.50
		10	0.10
		100	0.01
		101	0
TRANSITION:	DURATION:	P(X)
	NEXT:	DIFF_STATE
		3	0
		4	0.05
		5	0.50
		10	0.90
		100	0.99
		101	1
```

##END TRANSITIONS
The ending state doesn’t contain any emissions or transitions from the state.   However, every state that has a possible transition to the end of the state must define a STANDARD transition to the END state.