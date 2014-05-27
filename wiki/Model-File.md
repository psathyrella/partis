The StochHMM model file contain the complete model topology (States, transitions) and values of the model. It can define which states call additional user functions.   However, user-functions must be defined before the executable is compiled.    

The model is composed of multiple sections (Model Header Line, Model Information, Track Symbol Definitions, Ambiguous Symbol Definitions, Scaling Values, Templated States, and State Definitions)

Divider lines shown as “====” or “####” are optional, but you’re encouraged to use them to make the model easier to read.

***The first line of a model file must contain the header.***

#Model Header (Required)

##Format

```
#STOCHHMM MODEL FILE
```


[Next: Model Information](Model-Information)
