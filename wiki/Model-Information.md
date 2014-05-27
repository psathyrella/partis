[Previous: Model File](Model-File)
## Model Information (Optional)

Each model can contain a name, description, creation date information, creation commands, author, numerical attribute or range.  Any of the categories are optional.

If a numerical attribute is assigned, a model function may be defined to choose the model based on the numerical attribute or Upper and Lower range. 

###Format

```
MODEL INFORMATION
================================================
NAME:	VALUE
DESCRIPTION:	VALUE
CREATION_DATE:	VALUE
CREATION_COMMAND:	VALUE
AUTHOR:	VALUE
NUM_ATTRIB:	VALUE
UPPER:	VALUE
LOWER:	VALUE
```

*NAME*: Name of the model  

*DESCRIPTION*: Description of the model  

*CREATION_DATE*: Date that the model was created  

*CREATION_COMMAND*: Command or script used to create the model  

*AUTHOR*: Who created the model  

*NUM_ATTRIB*: Numerical value associated with the model.   When using multiple models, StochHMM could select the model based on the value associated with the model.  For example, if you have multiple models trained for sequences with different GC percentages.  StochHMM could select the model based on the GC percentage of the sequence.  

*UPPER*:  Range value, similar to NUM_ATTRIB except defines the inclusive upper range of values associated with the model.  *If used LOWER must also be supplied.*  

*LOWER*:  Range value, similar to NUM_ATTRIB except defines the inclusive lower range of values associated with the model.  *If used UPPER must also be supplied.*  

###Example

```
#STOCHHMM MODEL FILE

MODEL INFORMATION
================================================
NAME:	Simplified Human Gene
DESCRIPTION:	Created using perl script
CREATION_DATE:	August 12,2009
CREATION_COMMAND:	./creation_script.pl fasta.fa....
AUTHOR:	Paul Lott
NUM_ATTRIB:	20
UPPER:	29
LOWER:	20
```


[Next: Track Symbol Definitions](Track-Symbol-Definitions)