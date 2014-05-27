##CpG Model

[Borodovsky, M., Ekisheva, S., Problems and Solutions in Biological Sequence Analysis. Cambride University Press, 2006, pg 80]


##Model File
```
#STOCHHMM MODEL FILE

<MODEL INFORMATION>
======================================================
MODEL_NAME:	CpG
MODEL_DESCRIPTION:	Taken from CH3 Durbin/Eddy 3.16, 3.17
MODEL_CREATION_DATE:	August 28,2009


<TRACK SYMBOL DEFINITIONS>
======================================================
TRACK1:	A,C,G,T

<STATE DEFINITIONS>
######################################################
STATE:	
	NAME:	INIT
TRANSITION:	STANDARD:	P(X)
	HIGH:	0.5
	LOW:	0.5

####################################################
STATE:	
	NAME:	HIGH
	GFF_DESC:	HIGH
	PATH_LABEL:	H
TRANSITION:	STANDARD:	P(X)
	HIGH:	0.5
	LOW:	0.5
	END:	1
EMISSION:	TRACK1: P(X)
	ORDER:	0
0.2	0.3	0.3	0.2

###################################################
STATE:
	NAME:	LOW
	GFF_DESC:	LOW
	PATH_LABEL:	L
TRANSITION:	STANDARD:	P(X)
	HIGH:	0.4
	LOW:	0.6
	END:	1
EMISSION:	TRACK1: P(X)
	ORDER:	0
0.3	0.2	0.2	0.3

##################################################
//END
```

##Sequence File
```
>Example 3.16 Example
GGCACTGAA
```