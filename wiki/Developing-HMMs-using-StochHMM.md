The StochHMM application allows you to quickly implement a HMM using a sequence and model file. Regardless of using the C++ library or application the first place to start is creating a text model file.

The model file includes metadata about the model, the alphabet, any ambiguous characters that will be in the data, the state definitions (transitions, emissions, associated functions)

Currently, the training of the models is handled by external scripts that label and count training sequence sets. The Baum-Welch, Viterbi Training, and other forms of unsupervised training are not supported by StochHMM. We do have plans to implement these in the future. The framework for the Baum-Welch code is implemented, but needs additional testing and optimization.
