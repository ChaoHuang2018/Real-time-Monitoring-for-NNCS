#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"
#include "domain_computation.h"
#include "dynamics_linearization.h"

using namespace flowstar;
using namespace std;

// void test()
// {
//     string nn_name = "systems_with_networks/Benchmark1/nn_1_sigmoid_new";
//     string act_name = "sigmoid";
//     NeuralNetwork nn(nn_name, act_name);

//     // input layer -> 1st hidden layer
//     Matrix<Interval> weight_0 = nn.get_layers()[0].get_weight();
//     Matrix<Real> W1(weight_0.rows(), weight_0.cols());
//     for (int i = 0; i < weight_0.rows(); i++)
//     {
//         for (int j = 0; j < weight_0.cols(); j++)
//         {
//             W1[i][j] = weight_0[i][j].sup();
//         }
//     }
//     // cout << W1 << endl;

//     Matrix<Interval> bias_0 = nn.get_layers()[0].get_bias();
//     Matrix<Real> b1(bias_0.rows(), bias_0.cols());
//     for (int i = 0; i < bias_0.rows(); i++)
//     {
//         for (int j = 0; j < bias_0.cols(); j++)
//         {
//             b1[i][j] = bias_0[i][j].sup();
//         }
//     }
//     // cout << b1 << endl;

//     // 1st hidden layer -> 2nd hidden layer
//     Matrix<Interval> weight_1 = nn.get_layers()[1].get_weight();
//     Matrix<Real> W2(weight_1.rows(), weight_1.cols());
//     for (int i = 0; i < weight_1.rows(); i++)
//     {
//         for (int j = 0; j < weight_1.cols(); j++)
//         {
//             W2[i][j] = weight_1[i][j].sup();
//         }
//     }

//     Matrix<Interval> bias_1 = nn.get_layers()[1].get_bias();
//     Matrix<Real> b2(bias_1.rows(), bias_1.cols());
//     for (int i = 0; i < bias_1.rows(); i++)
//     {
//         for (int j = 0; j < bias_1.cols(); j++)
//         {
//             b2[i][j] = bias_1[i][j].sup();
//         }
//     }

//     // 2nd hidden layer -> output layer

//     Matrix<Interval> weight_2 = nn.get_layers()[2].get_weight();
//     Matrix<Real> W3(weight_2.rows(), weight_2.cols());
//     for (int i = 0; i < weight_2.rows(); i++)
//     {
//         for (int j = 0; j < weight_2.cols(); j++)
//         {
//             W3[i][j] = weight_2[i][j].sup();
//         }
//     }

//     Matrix<Interval> bias_2 = nn.get_layers()[2].get_bias();
//     Matrix<Real> b3(bias_2.rows(), bias_2.cols());
//     for (int i = 0; i < bias_2.rows(); i++)
//     {
//         for (int j = 0; j < bias_2.cols(); j++)
//         {
//             b3[i][j] = bias_2[i][j].sup();
//         }
//     }

//     // the input set is an interval
//     Interval I(0.7, 0.9), J(0.7, 0.9);
//     vector<Interval> input_set;
//     input_set.push_back(I);
//     input_set.push_back(J);

//     // setting the overall max order
//     unsigned int maxOrder = 10;
//     Global_Computation_Setting g_setting;
//     g_setting.prepareForReachability(maxOrder);

//     // the order in use
//     unsigned int order = 5;
//     Interval cutoff_threshold(-1e-12, 1e-12);

//     // translating the interval input set to a Taylor model vector
//     vector<Interval> tmv_domain;
//     TaylorModelVec<Real> tmv_input(input_set, tmv_domain), tmvTemp;

//     // displaying the Taylor model vector
//     stateVars.declareVar("y1");
//     stateVars.declareVar("y2");

//     tmVars.declareVar("t");
//     tmVars.declareVar("x1");
//     tmVars.declareVar("x2");

//     // reachability of input layer -> 1st hidden layer
//     TaylorModelVec<Real> tmv_1 = W1 * tmv_input;
//     tmv_1 += b1;

//     // for (int i = 0; i < tmv_1.tms.size(); i++)
//     // {
//     //     tmv_1.tms[i].output(cout, tmVars);
//     //     cout << endl;
//     // }

//     // applying the activation function
//     tmv_1.sigmoid_taylor(tmvTemp, tmv_domain, order, cutoff_threshold, g_setting);
//     tmv_1 = tmvTemp;

//     // reachability of 1st hidden layer -> 2nd hidden layer
//     TaylorModelVec<Real> tmv_2 = W2 * tmv_1;
//     tmv_2 += b2;

//     // cout << "Second layer:" << endl;
//     // for (int i = 0; i < tmv_2.tms.size(); i++)
//     // {
//     //     tmv_2.tms[i].output(cout, tmVars);
//     //     cout << endl;
//     // }

//     // applying the activation function
//     tmv_2.sigmoid_taylor(tmvTemp, tmv_domain, order, cutoff_threshold, g_setting);
//     tmv_2 = tmvTemp;

//     // reachability of 2nd hidden layer -> output layer
//     TaylorModelVec<Real> tmv_output = W3 * tmv_2;
//     tmv_output += b3;

//     // cout << "Third layer:" << endl;
//     // for (int i = 0; i < tmv_output.tms.size(); i++)
//     // {
//     //     tmv_output.tms[i].output(cout, tmVars);
//     //     cout << endl;
//     // }

//     // applying the activation function
//     tmv_output.sigmoid_taylor(tmvTemp, tmv_domain, order, cutoff_threshold, g_setting);
//     tmv_output = tmvTemp;

//     Matrix<Real> offset(1, 1);
//     offset[0][0] = -nn.get_offset().sup();
//     tmv_output += offset;

//     Matrix<Real> scalar(1, 1);
//     scalar[0][0] = nn.get_scale_factor().sup();
//     tmv_output = scalar * tmv_output;

//     tmv_input.output(cout, stateVars, tmVars);

//     tmv_output.output(cout, stateVars, tmVars);
//     cout << endl;

//     cout << endl;

//     tmv_output.output(cout, stateVars, tmVars);

//     vector<Interval> result;

//     tmv_output.intEval(result, tmv_domain);

//     cout << result[0] << endl;
// }

int main(int argc, char *argv[])
{
    clock_t begin, end;
    begin = clock();

    string str(argv[1]);
    string nn_name = "systems_with_networks/" + str;
    string act_name(argv[2]);
    NeuralNetwork nn(nn_name, act_name);

    vector<string> state_vars;
    state_vars.push_back("x0");
    state_vars.push_back("x1");

    NNTaylor nn_taylor(nn);
    /* The follow is used to test nn parser */
    vector<Interval> network_input_box;
    network_input_box.push_back(Interval(0.79, 0.81));
    network_input_box.push_back(Interval(0.79, 0.81));
    nn_taylor.set_taylor_linear(state_vars, network_input_box);
    nn_taylor.set_range_by_IBP(network_input_box);
    // nn_taylor.set_range_by_IBP(network_input_box);
    cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
    cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;
    cout << "Output Range of nn controller: " << nn_taylor.get_range_by_IBP() << endl;
    /* The above is used to test nn parser */

    /* The follow is used to test TMP */
    cout << endl;
    NNTaylor nn_taylor_1(nn);

    unsigned int maxOrder = 15;
    Global_Computation_Setting g_setting;
    g_setting.prepareForReachability(maxOrder);

    // the order in use
    unsigned int order = 5;
    Interval cutoff_threshold(-1e-12, 1e-12);
    unsigned int bernstein_order = 8;
    unsigned int partition_num = 2000;

    vector<Interval> tmv_domain;

    // TaylorInfo ti(g_setting, order, cutoff_threshold, tmv_domain);

    // nn_taylor_1.set_tm_by_TMP(state_vars, network_input_box, ti);

    // stateVars.declareVar("y");

    // tmVars.declareVar("t");
    // for (int i = 0; i < state_vars.size(); i++)
    // {
    //     tmVars.declareVar("x_" + to_string(i));
    // }

    // TaylorModelVec<Real> tmv_input = nn_taylor_1.get_input_tmv();
    // tmv_input.output(cout, stateVars, tmVars);

    // cout << endl;

    // TaylorModel<Real> tm_output = nn_taylor_1.get_output_tm();
    // vector<TaylorModel<Real>> temp_tm_v;
    // temp_tm_v.push_back(tm_output);
    // TaylorModelVec<Real> tmv_output(temp_tm_v);

    // tmv_output.output(cout, stateVars, tmVars);
    /* The above is used to test TMP */

    /* The following is used to test another function */
    cout << endl;

    // displaying the Taylor model vector
    stateVars.declareVar("y1");
    stateVars.declareVar("y2");

    tmVars.declareVar("t");
    tmVars.declareVar("x1");
    tmVars.declareVar("x2");

    begin = clock();

    TaylorInfo ti_1(g_setting, order, bernstein_order, partition_num, cutoff_threshold);
    TaylorModelVec<Real> tmv_input(network_input_box, tmv_domain);
    TaylorModelVec<Real> tmv_output_2;
    nn_taylor_1.get_output_tmv(tmv_output_2, tmv_input, ti_1, tmv_domain);
    tmv_output_2.output(cout, stateVars, tmVars);

    end = clock();

    cout << "Computition time for TMP: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

    vector<Interval> result;

    tmv_output_2.intEval(result, tmv_domain);

    cout << result[0] << endl;
    /* The above is used to test another function */

    // cout << endl;
    // test();

    return 0;
}
