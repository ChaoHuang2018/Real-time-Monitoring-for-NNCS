#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"
#include "domain_computation.h"
#include "dynamics_linearization.h"

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

	NNTaylor nn_taylor(nn);

	clock_t begin, end;
	begin = clock();

	Matrix<Interval> remainder(2, 1), stateSpace(2, 1), u(1, 1), d(2, 1);
	Interval x0(-2, 2), x1(-2, 2), control(-1, 1), disturbance(-0.1, 0.1);

	stateSpace[0][0] = x0;
	stateSpace[1][0] = x1;

	u[0][0] = control;
	d[0][0] = d[1][0] = disturbance;

	// set the current state
	Matrix<Interval> x_current(2, 1);
	x_current[0][0] = 0.5;
	x_current[1][0] = 0.2;

	// get the primary guess of the reachable set in 0.1 second by 2 iterations
	Matrix<Interval> domain(2, 1);
	remainderEval_benchmark1(remainder, stateSpace, u, d, 0.1);
	int N = 2;
	for (int i = 0; i < N; ++i)
	{
		stateSpace = domain;
		remainderEval_benchmark1(remainder, stateSpace, u, d, 0.1);
		linearizationDomainEval_benchmark1(domain, x_current, remainder, u, d, 0.1);
	}
	cout << "Initial guess of the reachable set:" << domain << endl;

	// get the linear taylor model and the corresponding remainder of the dynamics over the primary guess
	Matrix<Interval> linearization_remainder(2, 1);
	Matrix<string> dynamics(2, 1);
	dynamics_linear_taylor_benchmark1(dynamics, domain);
	remainder_linear_taylor_benchmark1(linearization_remainder, domain);
	cout << "linear taylor expansion of dynamics: " << dynamics[0][0] << ", " << dynamics[1][0] << endl;
	cout << "linear taylor Remainder of dynamics: " << linearization_remainder << endl;

	// get the linear taylor model and the corresponding remainder of the nn controller over the primary guess
	vector<Interval> network_input_box;
	network_input_box.push_back(domain[0][0]);
	network_input_box.push_back(domain[1][0]);
	nn_taylor.set_taylor_linear(state_vars, network_input_box);
	// nn_taylor.set_range_by_IBP(network_input_box);
	end = clock();
	cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
	cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;
	cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

	return 0;
}
