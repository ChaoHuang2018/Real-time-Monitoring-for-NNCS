//#include "../flowstar-template/TaylorModel.h"
#include "../flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	// input layer -> 1st hidden layer
	Matrix<Real> W1(2, 1), b1(2, 1);

	W1[0][0] = 1;
	W1[1][0] = 2;

	b1[0][0] = -1;
	b1[1][0] = 1;


	// 1st hidden layer -> 2nd hidden layer
	Matrix<Real> W2(2, 2), b2(2, 1);

	W2[0][0] = 1;
	W2[1][0] = 2;
	W2[0][1] = 3;
	W2[1][1] = 4;

	b2[0][0] = -1;
	b2[1][0] = 1;


	// 2nd hidden layer -> output layer
	Matrix<Real> W3(1, 2), b3(1, 1);

	W3[0][0] = 1;
	W3[0][1] = 2;

	b3[0][0] = 5;


	// the input set is an interval
	Interval I(0,5);
	vector<Interval> input_set;
	input_set.push_back(I);


	// setting the overall max order
	unsigned int maxOrder = 15;
	Global_Computation_Setting g_setting;
	g_setting.prepareForReachability(maxOrder);


	// the order in use
	unsigned int order = 5;
	Interval cutoff_threshold(-1e-12, 1e-12);

	unsigned int bernstein_order = 10;
	unsigned int partition_num = 20;

	// translating the interval input set to a Taylor model vector
	vector<Interval> tmv_domain;
	TaylorModelVec<Real> tmv_input(input_set, tmv_domain), tmvTemp;


	// reachability of input layer -> 1st hidden layer
	TaylorModelVec<Real> tmv_1 = W1 * tmv_input;
	tmv_1 += b1;

	// applying the activation function
	tmv_1.sigmoid_taylor(tmvTemp, tmv_domain, order, bernstein_order, partition_num, cutoff_threshold, g_setting);
	tmv_1 = tmvTemp;


	// reachability of 1st hidden layer -> 2nd hidden layer
	TaylorModelVec<Real> tmv_2 = W2 * tmv_1;
	tmv_2 += b2;

	// applying the activation function
	tmv_2.sigmoid_taylor(tmvTemp, tmv_domain, order, bernstein_order, partition_num, cutoff_threshold, g_setting);
	tmv_2 = tmvTemp;


	// reachability of 2nd hidden layer -> output layer
	TaylorModelVec<Real> tmv_output = W3 * tmv_2;
	tmv_output += b3;





	// displaying the Taylor model vector
	stateVars.declareVar("y");

	tmVars.declareVar("t");
	tmVars.declareVar("x");

	tmv_input.output(cout, stateVars, tmVars);

	cout << endl;

	tmv_output.output(cout, stateVars, tmVars);

	return 0;
}
