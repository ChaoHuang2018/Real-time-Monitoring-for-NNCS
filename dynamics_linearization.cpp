#include "../flowstar-template/Continuous.h"

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

	Real factor(1.0 / 24.0), c1(-6), c2(84), c3(-105), c4(-12), c5(21);

	remainder[0][0] = c1 * domain[0][0] * (domain[0][0] - center[0][0]) * (domain[0][0] - center[0][0]);
	remainder[1][0] = 0;
}

vector<string> dynamics_linear_taylor_benchmark1(Matrix<Interval> &domain)
{
	Matrix<double> center(2, 1);
	for (int j = 0; j < 2; j++)
	{
		center[j][0] = domain[j][0].midpoint();
	}

	vector<string> dynamics;

	double const_term_ode1 = center[1][0] - pow(center[0][0], 3.0);
	vector<double> jacobian;
	jacobian.push_back(-3 * pow(center[0][0], 2.0));
	jacobian.push_back(1);
	string ode1_linear = to_string(const_term_ode1) + " + " + to_string(jacobian[0]) + " * " + "x0" + " + " + to_string(jacobian[1]) + " * " + "x1";

	string ode2_linear = "0";

	dynamics.push_back(ode1_linear);
	dynamics.push_back(ode2_linear);

	return dynamics;
}

int main()
{

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	Matrix<Interval> remainder(4, 1), stateSpace(4, 1), u(1, 1), d(4, 1);
	Interval x(-2, 2), y(-2, 2), z(-2, 2), s(-2, 2), control(-1, 1), disturbance(-0.1, 0.1);
	stateSpace[0][0] = x;
	stateSpace[1][0] = y;
	stateSpace[2][0] = z;
	stateSpace[3][0] = s;

	u[0][0] = control;

	d[0][0] = d[1][0] = d[2][0] = d[3][0] = disturbance;

	remainderEval_benchmark6(remainder, stateSpace, u, d, 0.1);

	cout << "Initial:\n"
		 << remainder << endl;

	Matrix<Interval> x0(4, 1);
	x0[0][0] = 0.5;
	x0[1][0] = 0.2;
	x0[2][0] = 0.4;
	x0[3][0] = 0.8;

	Matrix<Interval> domain(4, 1);

	linearizationDomainEval_benchmark6(domain, x0, remainder, u, d, 0.1);

	cout << domain << endl;

	int N = 5;

	for (int i = 0; i < N; ++i)
	{
		stateSpace = domain;

		remainderEval_benchmark6(remainder, stateSpace, u, d, 0.1);

		linearizationDomainEval_benchmark6(domain, x0, remainder, u, d, 0.1);
	}

	cout << "\nAfter " << N << "-round refinement:\n"
		 << remainder << endl;
	cout << domain << endl;

	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	return 0;
}
