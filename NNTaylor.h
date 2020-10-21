#include "flowstar-template/Continuous.h"
#include "Neuron.h"
#include "NeuralNetwork.h"

using namespace flowstar;
using namespace std;

class NNTaylor
{
protected:
    NeuralNetwork nn;
    string taylor_linear_expression;
    Interval taylor_linear_remainder;
    Interval output_range_IBP;

public:
    NNTaylor();
    NNTaylor(NeuralNetwork nn);

    void set_taylor_linear(vector<string> state_vars, vector<Interval> network_input_box);
    void set_range_by_IBP(vector<Interval> network_input_box);

    string get_taylor_expression();
    Interval get_taylor_remainder();
    Interval get_range_by_IBP();

    static double inf_norm(Matrix<double> m)
    {
        double norm = 0;
        for (int i = 0; i < m.rows(); i++)
        {
            double row_abs_sum = 0;
            for (int j = 0; j < m.cols(); j++)
            {
                row_abs_sum = row_abs_sum + abs(m[i][j]);
            }
            if (norm <= row_abs_sum)
            {
                norm = row_abs_sum;
            }
        }
        return norm;
    }
};