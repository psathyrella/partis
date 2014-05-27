Scaling Functions section allows you define more complex scaling function for scores.

Section is define in the header section of the model before the STATE definitions.

Each **SCALE** defines a function that can be used throughout the model to scale values received by external functions.


**SCALE:** Definition of scaling function name (must be unique)

**VALUE:** Definition of x values

**SCALED:** Definition of transformed scores f(x)


Optional:

**MIN_VALUE:** Define the minimum value allowed.  Any value lower will use the minimum

**MIN_SCALED:** Define the score to be used when MIN_VALUE is passed.

**MAX_VALUE:** Define the maximum value allowed.  Any value higher will use this value.

**MAX_SCALED:** Define the score to be used when MAX_VALUE is passed.


If values fall in between x values used to define the function, the value will be interpolated.   Values beyond the first or last values provided use linear extrapolation to calculate the scaled score.   On the other hand, if a max and min values are defined, then the minimum/maximum scaled score will be used if the values is out of range of the min and max.


##Format
```
SCALING FUNCTIONS
====================================================================================================
SCALE:	<function name>
	VALUE:	<value_type	[<List of x Values>]
	SCALED:	<value_type>	[<List of f(x) Values>]
	MIN_VALUE:	<value_type	<value>
	MIN_SCALED:	<value_type>	<value>
	MAX_VALUE:	<value_type>	<value>
	MAX_SCALED:	<value_type>	<value>
```

##Example
This example creates two scaling function for use in the model. HMMER defines a max and min value to be used by the function.  Whereas, SIDD does not.   So values higher or lower than what has been chosen will use linear extrapolation to calculate the scaled score. 

```
SCALING FUNCTIONS
====================================================================================================
SCALE:	HMMER
	VALUE:	P(X)	[1,2,3,4,5,6]
	SCALED:	LOG	[2,4,6,8,10,12]
	MIN_VALUE:	P(X)	0
	MIN_SCALED:	LOG	0
	MAX_VALUE:	P(X)	7
	MAX_SCALED:	LOG	14
SCALE:	SIDD
	VALUE:	LOG	[0,-0.2,-0.4,-0.6,-0.8,-1]
	SCALED:	P(X)	[0,1,2,3,4,5]
```
