[Previous: Model Information](https://github.com/KorfLab/StochHMM/wiki/Model-Information)

##Track Symbol Definitions (Required)
StochHMM digitizes the character or word to an integer value.  Track symbol definitions requires a name and symbols(characters or strings) used by the track. A track can have up to 255 discrete symbols, all symbols must be explicitly defined.  Characters or Strings of a track are not limited to 1 character, but can be any size and variable size.  As in the USER track below.

###Note:  StochHMM considers uppercase and lowercase to be unique.   If you define uppercase and then use lowercase sequence, your lowercase will be either not accepted as valid alphabet or be assigned as the default ambiguous character.

**If variable size words are used, please see [Sequence Files](Sequence-Files) for formatting requirements***

If there are ambiguous characters, [see Ambigious Symbol Definitions.](Ambiguous-Symbol-Definitions)

If a track is a real number track, instead of symbols use REAL_NUMBER.  [See Sequence-Files](Sequence-Files) for formatting requirements.


##Format
```
<TRACK SYMBOL DEFINITIONS>
======================================================
TRACK_NAME: TRACK_SYMBOLS
TRACK_NAME: REAL_NUMBER
...
```

##EXAMPLE
For a single alphanumerical track called SEQ, with characters A,C,G,T
```
<TRACK SYMBOL DEFINITIONS>
======================================================
SEQ:	A,C,G,T 
```

For a single REAL_NUMBER track called PWM
```
<TRACK SYMBOL DEFINITIONS>
======================================================
PWM:	REAL_NUMBER 
```

For multiple Tracks,
```
<TRACK SYMBOL DEFINITIONS>
======================================================
SEQ:	A,C,G,T 
USER:	HOW,I,ARE,WE   
PWM:	REAL_NUMBER 
```



##TRACK FUNCTION (Compile Time Function)
A Track Function can create a REAL_NUMBER track from a alphanumeric track, but must be defined at compile-time rather than at runtime.  

##FORMAT
```
<TRACK SYMBOL DEFINITIONS>
======================================================
SEQ:	A,C,G,T      
TRACK_NAME:	TRACK_FUNCTION:	TRACK_TO_USE	SCALE:	SCALING_FACTOR
```

##EXAMPLE
```
TRACK SYMBOL DEFINITIONS
======================================================
SEQ:		A,C,G,T      
SIDD:		siddAnalysis:	SEQ	SCALE:	0.9
```
SEQ is an alphanumerical track with four characters A,C,G, and T.

SIDD is a REAL_NUMBER track but is defined by a user-defined function called siddAnalysis.     `siddAnalysis: SEQ` means that the siddAnalysis function will use the SEQ sequence data to determine the SIDD track values.  The values provided by siddAnalysis will multiplied by the SCALE number.  Here we’ll multiply every value by 0.9.

If a user defined weighting function is defined then the weights section must be defined before the track section of the model. [See Scaling Values](Scaling-Values)

At run-time StochHMM will import the SEQ sequence and use siddAnalysis with SEQ to determine the SIDD track values.  

[Next: Ambiguous Symbol Definition](Ambiguous-Symbol-Definitions)