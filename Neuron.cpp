#include "flowstar-template/Continuous.h"
#include "Neuron.h"

using namespace flowstar;
using namespace std;

Neuron::Neuron()
{
}

Neuron::Neuron(int nn_input_dim)
{
    this->nn_input_dim = nn_input_dim;
}

int Neuron::get_nn_input_dim()
{
    return this->nn_input_dim;
}

Activation Neuron::get_activation_info()
{
    return this->activation_info;
}

Interval Neuron::get_input_value()
{
    return Interval(this->input_value);
}

Interval Neuron::get_input_range()
{
    return Interval(this->input_range);
}

Matrix<Interval> Neuron::get_first_order_der_value()
{
    return Matrix<Interval>(this->first_order_der_value);
}

Matrix<Interval> Neuron::get_first_order_der_range()
{
    return Matrix<Interval>(this->first_order_der_range);
}

Matrix<Interval> Neuron::get_second_order_der_value()
{
    return Matrix<Interval>(this->second_order_der_value);
}

Matrix<Interval> Neuron::get_second_order_der_range()
{
    return Matrix<Interval>(this->second_order_der_range);
}

void Neuron::set_input_value(Interval input_value)
{
    this->input_value = Interval(input_value);
}

void Neuron::set_input_range(Interval input_range)
{
    this->input_range = Interval(input_range);
}

void Neuron::set_first_order_der_value(Matrix<Interval> first_order_der_value)
{
    this->first_order_der_value = Matrix<Interval>(first_order_der_value);
}

void Neuron::set_first_order_der_range(Matrix<Interval> first_order_der_range)
{
    this->first_order_der_range = Matrix<Interval>(first_order_der_range);
}

void Neuron::set_second_order_der_value(Matrix<Interval> second_order_der_value)
{
    this->second_order_der_value = Matrix<Interval>(second_order_der_value);
}

void Neuron::set_second_order_der_range(Matrix<Interval> second_order_der_range)
{
    this->second_order_der_range = Matrix<Interval>(second_order_der_range);
}

void Neuron::set_input_value(vector<Neuron> last_layer_info, Matrix<Interval> weight, Interval bias)
{
    Interval input_value(bias);
    for (int j = 0; j < last_layer_info.size(); j++)
    {
        // cout << weight[0][j] << endl;
        // cout << last_layer_info[j].activation_info.get_value() << endl;
        input_value = input_value + weight[0][j] * last_layer_info[j].activation_info.get_value();
        // cout << "iterative value: " << input_value << endl;
    }
    this->input_value = input_value;
}

void Neuron::set_input_range(vector<Neuron> last_layer_info, Matrix<Interval> weight, Interval bias)
{
    Interval input_range(bias);
    for (int j = 0; j < last_layer_info.size(); j++)
    {
        // cout << weight[0][j] << endl;
        input_range = input_range + weight[0][j] * last_layer_info[j].activation_info.get_output_range();
        // cout << "iterative range: " << input_range << endl;
    }
    this->input_range = input_range;
}

void Neuron::set_first_order_der_value(vector<Neuron> last_layer_info, Matrix<Interval> weight)
{
    Matrix<Interval> first_order_der_value(this->nn_input_dim, 1);
    for (int i = 0; i < this->nn_input_dim; i++)
    {
        Interval first_order_der_value_i(0);
        for (int j = 0; j < last_layer_info.size(); j++)
        {
            first_order_der_value_i = first_order_der_value_i + (weight[0][j] * last_layer_info[j].activation_info.get_de()) * last_layer_info[j].get_first_order_der_value()[i][0];
        }
        first_order_der_value[i][0] = first_order_der_value_i;
    }
    this->first_order_der_value = first_order_der_value;
}

void Neuron::set_first_order_der_range(vector<Neuron> last_layer_info, Matrix<Interval> weight)
{
    Matrix<Interval> first_order_der_range(this->nn_input_dim, 1);
    for (int i = 0; i < this->nn_input_dim; i++)
    {
        Interval first_order_der_range_i(0);
        for (int j = 0; j < last_layer_info.size(); j++)
        {
            first_order_der_range_i = first_order_der_range_i + (weight[0][j] * last_layer_info[j].activation_info.get_de_range()) * last_layer_info[j].get_first_order_der_range()[i][0];
        }
        first_order_der_range[i][0] = first_order_der_range_i;
    }
    this->first_order_der_range = first_order_der_range;
}

void Neuron::set_second_order_der_value(vector<Neuron> last_layer_info, Matrix<Interval> weight)
{
    Matrix<Interval> second_order_der_value(this->nn_input_dim, this->nn_input_dim);
    for (int i = 0; i < this->nn_input_dim; i++)
    {
        for (int s = i; s < this->nn_input_dim; s++)
        {
            Interval second_order_der_value_is(0);
            for (int j = 0; j < last_layer_info.size(); j++)
            {
                second_order_der_value_is = second_order_der_value_is + weight[0][j] * (last_layer_info[j].activation_info.get_de2() * *(last_layer_info[j].get_first_order_der_value()[i]) * *(last_layer_info[j].get_first_order_der_value()[s]) + last_layer_info[j].get_activation_info().get_de() * (last_layer_info[j].get_second_order_der_value()[i][s]));
            }
            second_order_der_value[i][s] = second_order_der_value_is;
            second_order_der_value[s][i] = second_order_der_value_is;
        }
    }
    this->second_order_der_value = second_order_der_value;
}

void Neuron::set_second_order_der_range(vector<Neuron> last_layer_info, Matrix<Interval> weight)
{
    Matrix<Interval> second_order_der_range(this->nn_input_dim, this->nn_input_dim);
    for (int i = 0; i < this->nn_input_dim; i++)
    {
        for (int s = i; s < this->nn_input_dim; s++)
        {
            Interval second_order_der_range_is(0);
            for (int j = 0; j < last_layer_info.size(); j++)
            {
                // cout << "111111" << endl;
                // cout << weight[0][j] << endl;
                // cout << "222222" << endl;
                // cout << last_layer_info[j].activation_info.get_de2_range() << endl;
                // cout << "333333" << endl;
                // cout << *(last_layer_info[j].get_first_order_der_range()[i]) << endl;
                // cout << "444444" << endl;
                // cout << *(last_layer_info[j].get_first_order_der_range()[s]) << endl;
                // cout << "555555" << endl;
                // cout << last_layer_info[j].get_activation_info().get_de_range() << endl;
                // cout << "666666" << endl;
                // cout << (last_layer_info[j].get_second_order_der_range()[i][s]) << endl;
                // cout << "777777" << endl;
                second_order_der_range_is = second_order_der_range_is + weight[0][j] * (last_layer_info[j].activation_info.get_de2_range() * *(last_layer_info[j].get_first_order_der_range()[i]) * *(last_layer_info[j].get_first_order_der_range()[s]) + last_layer_info[j].get_activation_info().get_de_range() * (last_layer_info[j].get_second_order_der_range()[i][s]));
            }
            second_order_der_range[i][s] = second_order_der_range_is;
            second_order_der_range[s][i] = second_order_der_range_is;
        }
    }
    this->second_order_der_range = second_order_der_range;
}

void Neuron::set_activation_info(string activation_type, string approach)
{
    this->activation_info = Activation(activation_type, this->input_value, this->input_range, approach);
}