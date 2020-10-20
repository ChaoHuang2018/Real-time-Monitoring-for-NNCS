#include "NNTaylor.h"

using namespace flowstar;
using namespace std;

NNTaylor::NNTaylor()
{
}

NNTaylor::NNTaylor(NeuralNetwork nn)
{
    this->nn = nn;
}

void NNTaylor::set_taylor_linear(vector<string> state_vars, vector<Interval> network_input_box)
{
    Matrix<double> center(this->nn.get_num_of_inputs(), 1);
    for (int j = 0; j < network_input_box.size(); j++)
    {
        center[j][0] = network_input_box[j].midpoint();
    }

    Matrix<double> half_len(this->nn.get_num_of_inputs(), 1);
    for (int j = 0; j < network_input_box.size(); j++)
    {
        half_len[j][0] = network_input_box[j].width();
    }

    // layer_info_all_layer stores all the information (value/range of output/jocobian/hessian) of all the layers
    vector<vector<Neuron>> layer_info_all_layer;

    // process the input layer
    vector<Neuron> layer_info;

    for (int i = 0; i < this->nn.get_num_of_inputs(); i++)
    {
        Matrix<Interval> first_order_der_value(this->nn.get_num_of_inputs(), 1);
        first_order_der_value[i][0] = Interval(1, 1);
        Matrix<Interval> first_order_der_range(this->nn.get_num_of_inputs(), 1);
        for (int j = 0; j < first_order_der_range.rows(); j++)
        {
            first_order_der_range[j][0] = Interval(1, 1);
        }

        Matrix<Interval> second_order_der_value(this->nn.get_num_of_inputs(), this->nn.get_num_of_inputs());
        Matrix<Interval> second_order_der_range(this->nn.get_num_of_inputs(), this->nn.get_num_of_inputs());

        Neuron neuron(this->nn.get_num_of_inputs());
        neuron.set_input_value(Interval(center[i][0], center[i][0]));
        neuron.set_input_range(network_input_box[i]);
        neuron.set_first_order_der_value(first_order_der_value);
        neuron.set_first_order_der_range(first_order_der_range);
        neuron.set_second_order_der_value(second_order_der_value);
        neuron.set_second_order_der_range(second_order_der_range);
        neuron.set_activation_info("Affine", "taylor");

        layer_info.push_back(neuron);
    }
    layer_info_all_layer.push_back(layer_info);

    // start to process hidden layers and output layer
    // may be incorrect since I did not define copy constructor
    vector<Neuron> last_layer_info = layer_info;
    for (int s = 0; s < this->nn.get_num_of_inputs() + 1; s++)
    {
        vector<Neuron> this_layer_info;
        Layer layer = this->nn.get_layers()[s];
        Matrix<Interval> weight = layer.get_weight();
        Matrix<Interval> bias = layer.get_bias();

        for (int i = 0; i < layer.get_neuron_number_this_layer(); i++)
        {
            Matrix<Interval> weight_i(layer.get_neuron_number_last_layer(), 1);
            weight.getRowVec(weight_i, i);
            Matrix<Interval> bias_i_matrix(1, 1);
            bias.getRowVec(bias_i_matrix, i);
            Interval bias_i = bias_i_matrix[0][0];

            Neuron neuron(this->nn.get_num_of_inputs());
            neuron.set_input_value(last_layer_info, weight_i, bias_i);
            neuron.set_input_range(last_layer_info, weight_i, bias_i);
            neuron.set_first_order_der_value(last_layer_info, weight_i);
            neuron.set_first_order_der_range(last_layer_info, weight_i);
            neuron.set_second_order_der_value(last_layer_info, weight_i);
            neuron.set_second_order_der_range(last_layer_info, weight_i);
            neuron.set_activation_info(layer.get_activation(), "taylor");
            this_layer_info.push_back(neuron);
        }

        layer_info_all_layer.push_back(this_layer_info);
        last_layer_info = this_layer_info;
    }

    // process the scalar and offset by constructing an addtional virtual single-neuron layer
    Interval scale_factor = this->nn.get_scale_factor();
    Interval offset = this->nn.get_offset();
    vector<Neuron> virtual_layer_info;
    Neuron virtual_neruon(this->nn.get_num_of_inputs());
    Matrix<Interval> virtual_weight(1, 1);
    virtual_weight[0][0] = scale_factor;
    Interval virtual_bias = -offset * scale_factor;
    virtual_neruon.set_input_value(last_layer_info, virtual_weight, virtual_bias);
    virtual_neruon.set_input_range(last_layer_info, virtual_weight, virtual_bias);
    virtual_neruon.set_first_order_der_value(last_layer_info, virtual_weight);
    virtual_neruon.set_first_order_der_range(last_layer_info, virtual_weight);
    virtual_neruon.set_second_order_der_value(last_layer_info, virtual_weight);
    virtual_neruon.set_second_order_der_range(last_layer_info, virtual_weight);
    virtual_neruon.set_activation_info("Affine", "taylor");

    Interval nn_output = virtual_neruon.get_input_value();
    cout << "output on center: " << nn_output.inf() << endl;
    Matrix<Interval> jacobian_value = virtual_neruon.get_first_order_der_value();
    cout << "output_der on center: ";
    for (int i = 0; i < jacobian_value.rows(); i++)
    {
        cout << jacobian_value[i][0].inf() << ", ";
    }
    cout << endl;

    string expression = to_string(nn_output.inf());
    for (int i = 0; i < state_vars.size(); i++)
    {
        expression = expression + " + " + to_string(jacobian_value[i][0].inf()) + " * " + state_vars[i];
    }
    this->taylor_linear_expression = expression;

    Matrix<Interval> hessian_range = virtual_neruon.get_second_order_der_range();
    Matrix<double> hessian_max(this->nn.get_num_of_inputs(), this->nn.get_num_of_inputs());
    for (int i = 0; i < hessian_max.rows(); i++)
    {
        for (int j = 0; j < hessian_max.cols(); j++)
        {
            hessian_max[i][j] = max(abs(hessian_range[i][j].inf()), abs(hessian_range[i][j].sup()));
        }
    }

    double error = 0.5 * inf_norm(half_len) * inf_norm(hessian_max);
    this->taylor_linear_remainder = Interval(-error, error);
}

void NNTaylor::set_range_by_IBP(vector<Interval> network_input_box)
{
    Matrix<double> center(this->nn.get_num_of_inputs(), 1);
    for (int j = 0; j < network_input_box.size(); j++)
    {
        center[j][0] = network_input_box[j].midpoint();
    }

    Matrix<double> half_len(this->nn.get_num_of_inputs(), 1);
    for (int j = 0; j < network_input_box.size(); j++)
    {
        half_len[j][0] = network_input_box[j].width();
    }

    // layer_info_all_layer stores all the information (value/range of output/jocobian/hessian) of all the layers
    vector<vector<Neuron>> layer_info_all_layer;

    // process the input layer
    vector<Neuron> layer_info;

    for (int i = 0; i < this->nn.get_num_of_inputs(); i++)
    {
        Neuron neuron(this->nn.get_num_of_inputs());
        neuron.set_input_value(Interval(center[i][0], center[i][0]));
        neuron.set_input_range(network_input_box[i]);
        neuron.set_activation_info("Affine", "IBP");

        layer_info.push_back(neuron);
    }
    layer_info_all_layer.push_back(layer_info);

    // start to process hidden layers and output layer
    // may be incorrect since I did not define copy constructor
    vector<Neuron> last_layer_info = layer_info;
    for (int s = 0; s < this->nn.get_num_of_inputs() + 1; s++)
    {
        vector<Neuron> this_layer_info;
        Layer layer = this->nn.get_layers()[s];
        Matrix<Interval> weight = layer.get_weight();
        Matrix<Interval> bias = layer.get_bias();

        for (int i = 0; i < layer.get_neuron_number_this_layer(); i++)
        {
            Matrix<Interval> weight_i(layer.get_neuron_number_last_layer(), 1);
            weight.getRowVec(weight_i, i);
            Matrix<Interval> bias_i_matrix(1, 1);
            bias.getRowVec(bias_i_matrix, i);
            Interval bias_i = bias_i_matrix[0][0];

            Neuron neuron(this->nn.get_num_of_inputs());
            neuron.set_input_value(last_layer_info, weight_i, bias_i);
            neuron.set_input_range(last_layer_info, weight_i, bias_i);
            neuron.set_activation_info(layer.get_activation(), "IBP");
            this_layer_info.push_back(neuron);
        }

        layer_info_all_layer.push_back(this_layer_info);
        last_layer_info = this_layer_info;
    }

    // process the scalar and offset by constructing an addtional virtual single-neuron layer
    Interval scale_factor = this->nn.get_scale_factor();
    Interval offset = this->nn.get_offset();
    vector<Neuron> virtual_layer_info;
    Neuron virtual_neruon(this->nn.get_num_of_inputs());
    Matrix<Interval> virtual_weight(1, 1);
    virtual_weight[0][0] = scale_factor;
    Interval virtual_bias = -offset * scale_factor;
    virtual_neruon.set_input_value(last_layer_info, virtual_weight, virtual_bias);
    virtual_neruon.set_input_range(last_layer_info, virtual_weight, virtual_bias);
    virtual_neruon.set_activation_info("Affine", "IBP");

    this->output_range_IBP = virtual_neruon.get_input_range();
}

string NNTaylor::get_taylor_expression()
{
    return this->taylor_linear_expression;
}

Interval NNTaylor::get_taylor_remainder()
{
    return this->taylor_linear_remainder;
}

Interval NNTaylor::get_range_by_IBP()
{
    return this->output_range_IBP;
}