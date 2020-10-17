#include "flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;

// Parse neural network and layer from a text file as classes
// Please provide the get and set function for each member in the two classes.

class NeuralNetwork{

    //
    protected:
        int num_of_inputs;
        int num_of_outputs;
        int num_of_hidden_layers;
        // use interval type for offset and scale_factor
        Interval offset;
        Interval scale_factor;
        // including hidden layers and output layer
        Layer layers[num_of_hidden_layers];

}

class Layer{
    protected:
        // hidden layer and output layer (input layer excluded)
        int neuron_number_last_layer;
        int neuron_number_this_layer;
        // activation of this layer: can be 'ReLU' or 'tanh' or 'sigmoid'
        string activation;
        // even though weight and bias are real matrix, we use interval to describe the access of each matrix for convenience
        Matrix<Interval> weight(neuron_number_this_layer)(neuron_number_last_layer);
        Matrix<Interval> bias(neuron_number_this_layer)(1);
}