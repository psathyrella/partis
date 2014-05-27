MotifScore imports a position weight matrix and scores the matrix across a sequence.

It is compatible with a simple motif (One that has no variable spacer region),  Undefined Spacer (Where there is no scoring of the spacer regions, but the spacer sizes can be defined), Defined Spacer (where each spacer position has it's own weight matrices).

The position weights can use zeroth order or higher order emission tables to calculate the matching scores.




##Track symbol Definitions - Required
Track symbol definitions are required.   Ambiguous definitions are optional
```
<TRACK SYMBOL DEFINITIONS>
======================================================
TRACK1:	A,C,G,T
```

## Background Definitions - Optional
A background definition is optional.

```
<BACKGROUND DEFINITION>
======================================================
EMISSION:	TRACK1:	COUNTS
	ORDER:	0
25	25	25	25
```


## Threshold Definitions - (Optional)
Creates a preliminary threshold definition for scoring the matrix.
If a position weight defines a separate threshold that threshold will
be used from that state onward.  
```
<THRESHOLD DEFINITION>
======================================================
-10
```

## Undefined Spacer Definition
If you want to define a variable length undefined spacer region, you need to identify the positions which transitions to or from the spacer.

```
<SPACER DEFINITIONS>
======================================================
H7-SPACER-N1: 12,20,22
```
The above transition creates a spacer of 12,20, or 22 bp between H7 position and the N1 position  (Where H7 and N1 are the names assigned to a particular position weight definition.

##Position Weight Definitions
Position weight definitions define the name of the position.  The can optionally set a threshold to use from that position onward.  The also need to identify which position comes next.  The emissions are the weight that will be used for that position.   See emissions for different types of values allowed and also allowing ambiguous characters.


```
<POSITION WEIGHT DEFINITIONS>
##################################################
NAME:	H1
TRANSITION:	H2
THRESHOLD: -100
EMISSION:	TRACK1:	COUNTS
	ORDER:	0
0	562	0	0
##################################################
...  ##### Additional Positions #####>
##################################################
//END
```