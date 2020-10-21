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
	Matrix<string> dynamics_string(2, 1);
	Matrix<double> dynamics_coeff(2, 2);
	Matrix<double> dynamics_const_term(2, 1);
	dynamics_linear_taylor_benchmark1(dynamics_string, dynamics_coeff, dynamics_const_term, domain);
	remainder_linear_taylor_benchmark1(linearization_remainder, domain);
	cout << "linear taylor expansion of dynamics: " << dynamics_string[0][0] << ", " << dynamics_string[1][0] << endl;
	cout << "linear taylor Remainder of dynamics: " << linearization_remainder << endl;

	// get the linear taylor model and the corresponding remainder of the nn controller over the primary guess
	vector<Interval> network_input_box;
	network_input_box.push_back(domain[0][0]);
	network_input_box.push_back(domain[1][0]);
	nn_taylor.set_taylor_linear(state_vars, network_input_box);
	// nn_taylor.set_range_by_IBP(network_input_box);
	cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
	cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;

	Matrix<double> nn_coeff(2, 2);
	Matrix<double> nn_const_term(2, 1);
	nn_coeff[1][0] = nn_taylor.get_jacobian()[0];
	nn_coeff[1][1] = nn_taylor.get_jacobian()[1];
	nn_const_term[1][0] = nn_taylor.get_output();

	// A is the coeff of linear term
	Matrix<double> A = dynamics_coeff + nn_coeff;
	// B is the coeff of constant term
	Matrix<double> B = dynamics_const_term + nn_const_term;
	// C is the remainder, including dynamics remainder, nn controller remainder, disturbace
	Matrix<Interval> C(2, 1);
	C[0][0] = linearization_remainder[0][0] + d[0][0];
	C[1][0] = linearization_remainder[1][0] + nn_taylor.get_taylor_remainder() + d[1][0];

	cout << "A: " << A << endl;
	cout << "B: " << C + B << endl;
	// cout << "C: " << C << endl;

	// compute reachable set for LTI
	// declare the number of variables
	unsigned int numVars = 2;

	int x0_id = stateVars.declareVar("x0");
	int x1_id = stateVars.declareVar("x1");

	// define the dynamics
	Linear_Time_Invariant_Dynamics dynamics(A, B);

	// set the reachability parameters
	Computational_Setting setting;

	// set the stepsize and the order
	setting.setFixedStepsize(0.02, 4);

	// set the time horizon
	setting.setTime(5);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// print out the computation steps: Can be turned off using printOff()
	setting.printOn();

	// call this function when all of the parameters are defined
	setting.prepare();

	// define the initial set which is a box
	Interval init_x0(x_current[0][0]), init_x1(x_current[0][1]);

	vector<Interval> box(numVars);
	box[x0_id] = init_x0;
	box[x1_id] = init_x1;

	Flowpipe initialSet(box);

	// unsafe set
	vector<Constraint> unsafeSet;
	Constraint constraint("-x0 + 2.2"); // x0 >= 2.2, change the constraint to x0 >= 2.1 will produce UNKNOWN
	unsafeSet.push_back(constraint);

	/*
	 * The structure of the class Result_of_Reachability is defined as below:
	 * nonlinear_flowpipes: the list of computed flowpipes
	 * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
	 * fp_end_of_time: the flowpipe at the time T
	 */
	Result_of_Reachability result;

	// run the reachability computation
	dynamics.reach(result, setting, initialSet, unsafeSet);

	switch (result.status)
	{
	case COMPLETED_SAFE:
		printf("Safe\n");
		break;
	case COMPLETED_UNSAFE:
		printf("Unsafe\n");
		break;
	case COMPLETED_UNKNOWN:
		printf("Unknown\n");
		break;
	default: // never happen to linear systems
		printf("Fail to compute flowpipes.\n");
	}

	end = clock();
	cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

	return 0;
}
