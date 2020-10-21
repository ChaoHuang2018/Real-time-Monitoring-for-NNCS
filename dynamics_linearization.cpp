#include "dynamics_linearization.h"

using namespace flowstar;
using namespace std;

// Benchmark 1
void remainder_linear_taylor_benchmark1(Matrix<Interval> &remainder, Matrix<Interval> &domain)
{
	Matrix<double> center(2, 1);
	for (int j = 0; j < 2; j++)
	{
		center[j][0] = domain[j][0].midpoint();
	}

	Real c1(-6);

	remainder[0][0] = c1 * domain[0][0] * (domain[0][0] - center[0][0]) * (domain[0][0] - center[0][0]);
	remainder[1][0] = 0;
}

void dynamics_linear_taylor_benchmark1(Matrix<string> &dynamics, Matrix<Interval> &domain)
{
	Matrix<double> center(2, 1);
	for (int j = 0; j < 2; j++)
	{
		center[j][0] = domain[j][0].midpoint();
	}

	double const_term_ode1 = center[1][0] - pow(center[0][0], 3.0);
	vector<double> jacobian;
	jacobian.push_back(-3 * pow(center[0][0], 2.0));
	jacobian.push_back(1);
	string ode1_linear = to_string(const_term_ode1) + " + (" + to_string(jacobian[0]) + ") * " + "x0" + " + (" + to_string(jacobian[1]) + ") * " + "x1";

	string ode2_linear = "0";

	dynamics[0][0] = ode1_linear;
	dynamics[1][0] = ode2_linear;
}