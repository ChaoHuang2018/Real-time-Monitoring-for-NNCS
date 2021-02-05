#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"
#include "domain_computation.h"
#include "dynamics_linearization.h"

using namespace flowstar;
using namespace std;

int main()
{
    clock_t begin, end;
    begin = clock();

    string nn_name = "nn_test";
    string act_name = "sigmoid";
    NeuralNetwork nn(nn_name, act_name);

    vector<string> state_vars;
    state_vars.push_back("x");

    NNTaylor nn_taylor(nn);

    /* The follow is used to test nn parser */
    vector<Interval> network_input_box;
    network_input_box.push_back(Interval(2.5, 2.5));
    nn_taylor.set_taylor_linear(state_vars, network_input_box);
    nn_taylor.set_range_by_IBP(network_input_box);
    // nn_taylor.set_range_by_IBP(network_input_box);
    cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
    cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;
    cout << "Output Range of nn controller: " << nn_taylor.get_range_by_IBP() << endl;
    /* The above is used to test nn parser */

    return 0;
}
