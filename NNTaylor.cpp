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
        half_len[j][0] = network_input_box[j].width() / 2.0;
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
        first_order_der_range[i][0] = Interval(1, 1);

        // Matrix<Interval> second_order_der_value(this->nn.get_num_of_inputs(), this->nn.get_num_of_inputs());
        // for (int p = 0; p < this->nn.get_num_of_inputs(); p++)
        // {
        //     for (int q = 0; q < this->nn.get_num_of_inputs(); q++)
        //     {
        //         second_order_der_value[p][q] = Interval(0);
        //     }
        // }
        Matrix<Interval> second_order_der_range(this->nn.get_num_of_inputs(), this->nn.get_num_of_inputs());
        // for (int p = 0; p < this->nn.get_num_of_inputs(); p++)
        // {
        //     for (int q = 0; q < this->nn.get_num_of_inputs(); q++)
        //     {
        //         second_order_der_range[p][q] = Interval(0);
        //     }
        // }

        Neuron neuron(this->nn.get_num_of_inputs());
        neuron.set_input_value(Interval(center[i][0], center[i][0]));
        // cout << neuron.get_input_value() << endl;
        neuron.set_input_range(network_input_box[i]);
        // cout << neuron.get_input_range() << endl;
        neuron.set_first_order_der_value(first_order_der_value);
        // cout << "Input layer, Neuron " << i << ", first order value: " << neuron.get_first_order_der_value() << endl;
        neuron.set_first_order_der_range(first_order_der_range);
        // cout << "Input layer, Neuron " << i << ", first order range: " << neuron.get_first_order_der_range() << endl;
        // neuron.set_second_order_der_value(second_order_der_value);
        neuron.set_second_order_der_range(second_order_der_range);
        // cout << "123: " << neuron.get_second_order_der_value() << endl;
        // cout << "456: " << neuron.get_second_order_der_range() << endl;
        neuron.set_activation_info("Affine", "taylor");
        // cout << "Input layer, Neuron " << i << ", first order value after activation: " << neuron.get_activation_info().get_de() << endl;
        // cout << "Input layer, Neuron " << i << ", first order range after activation: " << neuron.get_activation_info().get_de_range() << endl;

        layer_info.push_back(neuron);
    }
    layer_info_all_layer.push_back(layer_info);

    // start to process hidden layers and output layer
    // may be incorrect since I did not define copy constructor
    vector<Neuron> last_layer_info = layer_info;
    for (int s = 0; s < this->nn.get_num_of_hidden_layers() + 1; s++)
    {
        // cout << "s: " << s << endl;

        vector<Neuron> this_layer_info;
        Layer layer = this->nn.get_layers()[s];
        Matrix<Interval> weight = layer.get_weight();
        Matrix<Interval> bias = layer.get_bias();

        for (int i = 0; i < layer.get_neuron_number_this_layer(); i++)
        {
            Matrix<Interval> weight_i(1, layer.get_neuron_number_last_layer());
            weight.getRowVec(weight_i, i);
            // cout << "Layer " << s << ", Neuron " << i << ", weight: " << weight_i << endl;
            Matrix<Interval> bias_i_matrix(1, 1);
            bias.getRowVec(bias_i_matrix, i);
            Interval bias_i = bias_i_matrix[0][0];

            Neuron neuron(this->nn.get_num_of_inputs());
            neuron.set_input_value(last_layer_info, weight_i, bias_i);
            neuron.set_input_range(last_layer_info, weight_i, bias_i);
            neuron.set_first_order_der_value(last_layer_info, weight_i);
            // cout << "Layer " << s << ", Neuron " << i << ", first order value: " << neuron.get_first_order_der_value() << endl;
            neuron.set_first_order_der_range(last_layer_info, weight_i);
            // cout << "Layer " << s << ", Neuron " << i << ", first order range: " << neuron.get_first_order_der_range() << endl;
            // neuron.set_second_order_der_value(last_layer_info, weight_i);
            neuron.set_second_order_der_range(last_layer_info, weight_i);
            neuron.set_activation_info(layer.get_activation(), "taylor");
            // cout << "Layer " << s << ", Neuron " << i << ", first order value after activation: " << neuron.get_activation_info().get_de() << endl;
            // cout << "Layer " << s << ", Neuron " << i << ", first order range after activation: " << neuron.get_activation_info().get_de_range() << endl;

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
    // virtual_neruon.set_second_order_der_value(last_layer_info, virtual_weight);
    virtual_neruon.set_second_order_der_range(last_layer_info, virtual_weight);
    virtual_neruon.set_activation_info("Affine", "taylor");

    Interval const nn_output = virtual_neruon.get_input_value();
    this->output = nn_output.inf();
    cout << "output on center: " << nn_output.inf() << endl;
    Matrix<Interval> jacobian_value = virtual_neruon.get_first_order_der_value();
    cout << "jocobian: " << jacobian_value << endl;
    cout << "output_der on center: ";
    for (int i = 0; i < jacobian_value.rows(); i++)
    {
        cout << jacobian_value[i][0].inf() << ", ";
    }
    cout << endl;
    for (int i = 0; i < jacobian_value.rows(); i++)
    {
        this->jacobian.push_back(jacobian_value[i][0].inf());
    }

    string expression = to_string(nn_output.inf());
    for (int i = 0; i < state_vars.size(); i++)
    {
        expression = expression + " + (" + to_string(jacobian_value[i][0].inf()) + ") * (" + state_vars[i] + " - (" + to_string(center[i][0]) + "))";
    }
    this->taylor_linear_expression = expression;

    Matrix<Interval> hessian_range = virtual_neruon.get_second_order_der_range();
    cout << "hessian_range: " << hessian_range << endl;
    Matrix<double> hessian_max(this->nn.get_num_of_inputs(), this->nn.get_num_of_inputs());
    for (int i = 0; i < hessian_max.rows(); i++)
    {
        for (int j = 0; j < hessian_max.cols(); j++)
        {
            hessian_max[i][j] = max(abs(hessian_range[i][j].inf()), abs(hessian_range[i][j].sup()));
        }
    }
    cout << inf_norm(half_len) << endl;
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
        // cout << neuron.get_input_value() << endl;
        neuron.set_input_range(network_input_box[i]);
        // cout << neuron.get_input_range() << endl;
        neuron.set_activation_info("Affine", "IBP");

        layer_info.push_back(neuron);
    }
    layer_info_all_layer.push_back(layer_info);

    // start to process hidden layers and output layer
    // may be incorrect since I did not define copy constructor
    vector<Neuron> last_layer_info = layer_info;
    // cout << "num_of_hidden_layers: " << this->nn.get_num_of_hidden_layers() << endl;
    for (int s = 0; s < this->nn.get_num_of_hidden_layers() + 1; s++)
    {
        // cout << "s: " << s << endl;

        vector<Neuron> this_layer_info;
        Layer layer = this->nn.get_layers()[s];
        Matrix<Interval> weight = layer.get_weight();
        Matrix<Interval> bias = layer.get_bias();

        // cout << "this layer neuron number : " << layer.get_neuron_number_this_layer() << endl;

        // cout << "Layer " << s << ", Weight: " << weight << endl;
        // cout << "Layer " << s << ", Bias: " << bias << endl;

        for (int i = 0; i < layer.get_neuron_number_this_layer(); i++)
        {
            Matrix<Interval> weight_i(1, layer.get_neuron_number_last_layer());
            weight.getRowVec(weight_i, i);
            // cout << "Layer " << s << ", Neuron " << i << ", weight: " << weight_i << endl;
            Matrix<Interval> bias_i_matrix(1, 1);
            bias.getRowVec(bias_i_matrix, i);
            Interval bias_i = bias_i_matrix[0][0];
            // cout << "Layer " << s << ", Neuron " << i << ", bias: " << bias_i << endl;

            Neuron neuron(this->nn.get_num_of_inputs());
            neuron.set_input_value(last_layer_info, weight_i, bias_i);
            // cout << "Layer " << s << ", Neuron " << i << ", input value: " << neuron.get_input_value() << endl;
            neuron.set_input_range(last_layer_info, weight_i, bias_i);
            // cout << "Layer " << s << ", Neuron " << i << ", input range: " << neuron.get_input_range() << endl;
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

    Interval const nn_output_value = virtual_neruon.get_input_value();
    cout << "output on center: " << nn_output_value.inf() << endl;
    Interval nn_output_range = virtual_neruon.get_input_range();
    cout << "output range: " << nn_output_range << endl;
}

string NNTaylor::get_taylor_expression()
{
    return this->taylor_linear_expression;
}

Interval NNTaylor::get_taylor_remainder()
{
    return this->taylor_linear_remainder;
}

double NNTaylor::get_output()
{
    return this->output;
}

vector<double> NNTaylor::get_jacobian()
{
    return this->jacobian;
}

Interval NNTaylor::get_range_by_IBP()
{
    return this->output_range_IBP;
}