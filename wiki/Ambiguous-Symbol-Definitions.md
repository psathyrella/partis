[Previous: Track Symbol Definitions](Track-Symbol-Definitions)
##Ambiguous Track Symbol Definitions (Optional)
If ambiguous characters or symbols are not defined then StochHMM will interpret any undefined symbol as an unknown character.  This will cause the model to not be able to find a valid path through the sequence.  

If ambiguous characters are defined, then the first defined ambiguous symbol will be used as the default ambiguous character. This means if an unknown symbol is present but doesn't match any of the other characters, then it will be assigned the default value. If we encountered an “Q” in the example below, it would be interpreted as an N.

In addition to defining the ambiguous characters, each state must have an **AMBIGUOUS** tag in the emission definition in order for ambiguous characters to be emitted from the state.[See State Emissions](State-Emissions)

If **AMBIGUOUS SYMBOL DEFINITION** is missing from the model, then ambiguous characters will not be allowed, even if an individual state is set to allow Ambiguous Characters.

###Format

```
<AMBIGUOUS SYMBOL DEFINITIONS>
======================================================
TRACK_NAME: AMBIGUOUS_CHARACTER[UNAMBIGUOUS CHARACTERS SEPARATED BY COMMAS]
```


###EXAMPLE
If SEQ track has A,C,G,T, then we can define N to represent A or C or G or T.  
R = A or G
Y = G or C
W = A or T
K = G or T
M = A or C
 so on...


```
<AMBIGUOUS SYMBOL DEFINITIONS>
======================================================
SEQ: N[A,C,G,T],R[A,G],Y[C,T],S[G,C],W[A,T],K[G,T],M[A,C],B[C,G,T]
```

The cost of using ambiguous characters is accrued at the time the model is imported.  The model will calculate emission values for all of the ambiguous characters at the time the model is imported.  Defining N,R,Y,S,W,K,M,B will greatly increase the size of the emission table.   

Therefore:

1. If the sequences your are analyzing dont have all the ambiguous characters, don't define them all.  

2. If you need to treat all unknown character as an ambiguous character, define one ambiguous characters.  By design the first ambiguous character will be used for all unknowns. 

[Next: Scaling Functions](Scaling-Functions)
