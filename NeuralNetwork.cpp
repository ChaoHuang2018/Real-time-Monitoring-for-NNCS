#include <fstream>
#include <iostream>
#include "NeuralNetwork.h"


Layer::create_layer(int last_layer_dim, int dim, string act, Matrix<Interval> w, Matrix<Interval> b){
    neuron_number_last_layer = last_layer_dim;
    neuron_number_this_layer = dim;
    activation = act;
    weight = w;
    bias = b;
}


NeuralNetwork::NeuralNetwork(string filename, string act){
    std::ifstream input(filename);
    std::string line;

    // Parse the structure of neural networks
    getline(input, line);
    num_of_inputs = stoi(line);
    getline(input, line);
    num_of_outputs = stoi(line);
    getline(input, line);
    num_of_hidden_layers = stoi(line);

    std::vector<int> network_structure (num_of_hidden_layers + 1, 0);
    for (int idx = 0; idx < num_of_hidden_layers; idx++){
        getline(input, line);
        network_structure[idx] = stoi(line);
    }
    num_of_outputs = network_structure.back();

    // Parse the input text file and store weights and bias

    // compute parameters of the input layer
    Matrix<Interval> weight0(network_structure[0])(num_of_inputs);
    Matrix<Interval> bias0(network_structure[0])(1);
    for (int i = 0; i < network_structure[0]; i++){
        for (int j = 0; j < num_of_inputs; j++){
            getline(input, line);
            double value = stod(line);
            Interval I(value, value)
            weight0[i][j] = I;
        }
        getline(input, line);
        double value = stod(line);
        Interval I(value, value);
        bias0[i][0] = I;
    }
    Layer input_layer;
    input_layer.create_layer(num_of_inputs, network_structure[0], act, weight0, bias0);
    layers.push_back(input_layer);

    // compute the parameters of hidden layers
    for (int layer_idx = 0; layer_idx < num_of_hidden_layers; layer_idx++){
        Matrix<Interval> weight(network_structure[layer_idx + 1])(network_structure[layer_idx]);
        Matrix<Interval> bias(network_structure[layer_idx + 1])(1);

        for (int i = 0; i < network_structure[layer_idx + 1]; i++){
            for (int j = 0; j < network_structure[layer_idx]; j ++){
                getline(input, line);
                double value = stod(line);
                Interval I(value, value)
                weight0[i][j] = I;
            }
            getline(input, line);
            double value = stod(line);
            Interval I(value, value);
            bias0[i][0] = I;
        }
        Layer hidden_layer;
        hidden_layer.create_layer(network_structure[layer_idx], network_structure[layer_idx + 1], act, weight, bias);
        layers.push_back(hidden_layer);
    }
    // Affine mapping of the output
    getline(input, line);
    double value = stod(line);
    Inverval I(value, value);
    offset = I;
    getline(input, line);
    double value = stod(line);
    Inverval I(value, value);
    scale_factor = I;
}
