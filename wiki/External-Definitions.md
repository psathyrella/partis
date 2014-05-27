#External Definitions
[Previous: Sequence Files](Sequence-Files)

An external definitions is a definition of a weight or explicit emission state for a particular sequence.   It can be a weighting of a states within a region of sequence.  Or a explicit definition of the states that will be used.   

Weightings will be applied to the emissions of the states at the current position.   If the path of the model doesn't pass through the given state at the positions defined the weight will not be applied to those positions. 

**Currently, the weights must be appended on to the end of a fasta entry.

##Weighting Regions
States within a region can be weighted up or down based on the state's name, label, or GFF label.  Weighted states can be identified by different classifiers: **STATE_NAME**, **STATE_LABEL**, or **STATE_GFF**.

**STATE_NAME**:	Name of state defined in the model

**STATE_LABEL**:	Label of state(s) defined in the model

**STATE_GFF**:	GFF Label of state(s) defined in the model


###Format
To weight a state(s) in the given position(s):

```
[EXDEF:	WEIGHTED	START: <start_position>	END: <end_position>	<CLASSIFIER_TYPE>:	<classifier>
VALUE:	<value>	VALUE_TYPE:	<value_type>]
```

Weighting by **State Name**
```
[EXDEF: WEIGHTED	START: <start position>	END: <end position>	STATE_NAME: <state_name>	
VALUE:	<value>	VALUE_TYPE:	<value type>]
```

Weighting by **State Label**
```
[EXDEF: WEIGHTED	START: <start position>	END: <end position>	STATE_LABEL: <state_label>
VALUE:	<value>	VALUE_TYPE:	<value type>]
```

Weighting by **State GFF Label**
```
[EXDEF: WEIGHTED	START: <start position>	END: <end position>	STATE_GFF: <state_GFF>
VALUE:	<value>	VALUE_TYPE:	<value type>]
```

##Example

The following will up weight the **H** state at position 100 to 120 by a value of 100. 
```
[EXDEF: WEIGHTED	START:	100	END:	120	STATE_NAME:	H	VALUE: 100	VALUE_TYPE:P(X)]
```

The following will down weight the **L** state at position 100 to 120 by a value of 0.001.
```
[EXDEF: WEIGHTED	START:	100	END:	120	STATE_NAME:	L	VALUE: 0.001	VALUE_TYPE:P(X)]
```

Using the **STATE_LABEL**, we can weight multiple states:

```
[EXDEF:	WEIGHTED	START:	1	END:	5	STATE_LABEL: H	VALUE:100	VALUE_TYPE:LOG]
```

This up weights all the states with a label **H** in position 1 to 5 by a natural log value of 100;


Using the **STATE_GFF**, we can weight multiple states:

```
[EXDEF: WEIGHTED	START:	50	END:	55	STATE_GFF:	LOW	VALUE: 100	VALUE_TYPE:P(X)]
```

This will up weight states with a GFF label of **LOW** by a probability of 100.


##Explicit definitions of regions states
If you already know the path of the states through a sequence region, you can explicitly define the path of the model through the sequence.   This will restrict the model through this region to the states defined.  The algorithm will only compute the probabilities for these states. 

***Note:***If the path is invalid according to model's grammar, you can end without a valid path of the model through the sequence.

State names must match state names defined in the model.   The number of states included in the **PATH OF STATE NAMES** must match the size of region defined by **START** and **END** positions.  This will be be equal to (START-END)+1  

###Format
```
[EXDEF:	ABSOLUTE	START:	<start position>	END:	<end position>	TRACE: <PATH OF STATE NAMES>	]
```


###Example

```
[EXDEF:	ABSOLUTE	START:	31	END:	40	TRACE:	H,H,H,H,H,H,H,H,L,L]
```

This explicitly defines the state path of the sequence through positions 31-40 to be: H,H,H,H,H,H,H,H,L,L .



##Fasta with External Definitions
```
>MFN2
GTAGGCGGGGCGAGCCGGCTGGGCTCAGGGTCCACCAGCTCACCCGGGTCGAGGGGCAAT
CTGAGGCGACTGGTGACGCGCTTATCCACTTCCCTCCTCCCGCCTCCCCCTGGGGTGGCG
CTCGCTGGTGACGTAGTGAGTGTGATGGCCGCCGCGAGGCCGGGAAGGTGAAG
[EXDEF:	ABSOLUTE	START:	31	END:	40	TRACE:	H,H,H,H,H,H,H,H,L,L]
[EXDEF:	ABSOLUTE	START:	11	END:	15	TRACE:	H,H,H,L,L]
[EXDEF:	WEIGHTED	START:	1	END:	5	STATE_LABEL: H	VALUE:100	VALUE_TYPE:LOG]
[EXDEF: WEIGHTED	START:	50	END:	55	STATE_GFF:	LOW	VALUE: 100	VALUE_TYPE:P(X)]
[EXDEF: WEIGHTED	START:	100	END:	120	STATE_NAME:	H	VALUE: 100	VALUE_TYPE:P(X)]
[EXDEF: WEIGHTED	START:	100	END:	120	STATE_NAME:	L	VALUE: 0.001	VALUE_TYPE:P(X)]
```