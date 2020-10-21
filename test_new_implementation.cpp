#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"

using namespace flowstar;
using namespace std;

int main()
{
	string nn_name = "nn/nn_1_sigmoid";
	string act_name = "sigmoid";
	NeuralNetwork nn(nn_name, act_name);

	vector<string> state_vars;
	state_vars.push_back("x0");
	state_vars.push_back("x1");

	vector<Interval> network_input_box;
	network_input_box.push_back(Interval(0.4, 0.5));
	network_input_box.push_back(Interval(0.4, 0.5));

	NNTaylor nn_taylor(nn);

	clock_t begin, end;
	begin = clock();
	// nn_taylor.set_taylor_linear(state_vars, network_input_box);
	nn_taylor.set_range_by_IBP(network_input_box);
	end = clock();
	// cout << "Linear Taylor Expression: " << nn_taylor.get_taylor_expression() << endl;
	// cout << "Linear Taylor Remainder: " << nn_taylor.get_taylor_remainder() << endl;
	// cout << "Time for generating error: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

	return 0;
}
