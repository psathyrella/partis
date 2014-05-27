##State Definitions (Required)
Defines each state within the model.  Each model is required to include an initial state (INIT = q0) and transitions to an ending state(END). Each state should have a name, a label, a transition, and an emission, excluding the INIT and END states.

The state definition section should begin with <STATE DEFINITIONS> line

###Format
```
<STATE DEFINITIONS>
#########################################################
```

After the definitions line begins the definitions for each state, beginning with INIT.

##Initial State (INIT)
Each model STATE DEFINITIONS section should begin with an initial state named INIT. This state doesn't have an emission.  It only defines valid transitions from the start position to other states.  
	
###FORMAT

```
<STATE DEFINITIONS>
#########################################################
STATE:
     NAME:	INIT
TRANSITION:	STANDARD:	VALUE_TYPE
     STATE_NAME:	VALUE
     ...
```

###EXAMPLE


```
<STATE DEFINITIONS>
#########################################################
STATE:
     NAME:	INIT
TRANSITION:	STANDARD:	P(X)
     UNDER:	0.5
     OVER:	0.5
```
	
Here we have the STATE DEFINITIONS  header and we’ve also defined a state named INIT.   INIT is the initial state of the model.  It essentially says at the beginning of the sequence, where is it possible to start the model.   We’ve defined transitions to the states named UNDER and OVER, meaning it is only possible to start with the state UNDER or OVER states.

Only STANDARD transitions are valid for the INIT state.  However, the values can be provided as either probabilities P(X) or LOG (natural log value).

## END State
In addition to the INIT state, each model is required to define which states can transition to an END state.   At least one state, excluding the INIT state must contain an transition to END state.   END state probabilities should be considered independent of other transitions and must be a STANDARD transition. [see Transitions](State-Transitions)
There should not be any State definition for the END state, only transitions to END.

If no transitions to END state are defined, viterbi decoding can't find a valid path through the sequence.


## States
Each state is divided by a line of ###... and must begin with STATE:
Following the state, we define the identifiers and essential parts of the state

###Format

```
#############################################
STATE:
	NAME:	NAME OF STATE
	PATH_LABEL:	SINGLE_CHARACTER_LABEL
	GFF_DESC:	GFF_DESCRIPTOR_NAME
	...
```

NAME: (REQUIRED)
Each state must be given a unique name.  This name will be used in defining transition from and to different states.

GFF_DESC: (OPTIONAL)
If you want to output the traceback as a GFF, you must include a GFF description for all the states you want to include in the output.   GFF descriptions don’t need to be unique.  If there are more than one state with the same descriptor the GFF coordinates will combine the states into one GFF descriptor.

PATH_LABEL: (REQUIRED)
Single character label that doesn’t need to be unique.  For example for an exon, we may give all the states that describe the exon a label of “E”.


###EXAMPLE
Excludes transition and emission definitions for the states.   See [Transitions](State-Transitions) and [Emissions](State-Emissions)

```
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
...
#############################################
STATE:
	NAME:	LOADED
	PATH_LABEL:	L
	GFF_DESC:	LOADED
...
#############################################
//END
```

This model defines two states besides INIT and END.

[Next: Transitions](State-Transitions)
