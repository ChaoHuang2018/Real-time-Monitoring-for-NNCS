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

	string nn_name = "nn/nn_1_sigmoid";
	string act_name = "sigmoid";
	NeuralNetwork nn(nn_name, act_name);

	vector<string> state_vars;
	state_vars.push_back("x0");
	state_vars.push_back("x1");

	NNTaylor nn_taylor(nn);

	/* The follow is used to test nn parser */
	vector<Interval> network_input_box;
	network_input_box.push_back(Interval(8.550138163576985e-01, 9.320589519204974e-01));
	network_input_box.push_back(Interval(2.827279999999999e-01, 1.202728000000000e+00));
	nn_taylor.set_taylor_linear(state_vars, network_input_box);
	// nn_taylor.set_range_by_IBP(network_input_box);
	cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
	cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;
	/* The above is used to test nn parser */

	return 0;
}
