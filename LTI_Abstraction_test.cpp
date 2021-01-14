#include "LTI_Abstraction.h"

using namespace flowstar;
using namespace std;


int main()
{
	// declare the number of variables
	unsigned int numVars = 2;

	int x_id = stateVars.declareVar("x");
	int y_id = stateVars.declareVar("y");


	// define the matrix A, B, C such that B and C are identity matrices
	Matrix<Real> A(numVars, numVars), B(numVars), C(numVars), a(numVars, 1);

	A[0][0] = -2;
	A[0][1] = -4;
	A[1][0] = 4;
	A[1][1] = -2;


	Interval I1(-0.3,0.3), I2(-0.02,0.01), I3(-0.05,0.05), I4(-0.05,0.05);

	Matrix<Interval> Ux(numVars, 1), Uu(numVars, 1);

	Ux[0][0] = I1;
	Ux[1][0] = I2;

	Uu[0][0] = I3;
	Uu[1][0] = I4;



	// define the abstraction
	LTI_Abstraction abstraction(A, B, a, C, Ux, Uu);


	// set the reachability parameters
	Computational_Setting setting;

	// set the stepsize and the order
	setting.setFixedStepsize(0.02, 4);

	// set the time horizon
	setting.setTime(10);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// print out the computation steps
	setting.printOn();


	// call this function when all of the parameters are defined
	setting.prepare();


	// define the initial set which is a box
	Interval init_x(1.9,2.1), init_y(-0.1,0.1);

	vector<Interval> box(numVars);
	box[x_id] = init_x;
	box[y_id] = init_y;

	Flowpipe initialSet(box);


	// unsafe set
	vector<Constraint> unsafeSet;
	Constraint constraint("-x + 2.2"); // x >= 2.2, change the constraint to x >= 2.1 will produce UNKNOWN
	unsafeSet.push_back(constraint);

	/*
	 * The structure of the class Result_of_Reachability is defined as below:
	 * nonlinear_flowpipes: the list of computed flowpipes
	 * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
	 * fp_end_of_time: the flowpipe at the time T
	 */
	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	abstraction.reach(result, setting, initialSet, unsafeSet);

	switch(result.status)
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
	default:  // never happen to linear systems
		printf("Fail to compute flowpipes.\n");
	}



	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	// The following code is not necessary if you do not need a plot

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting, initialSet);

	Plot_Setting plot_setting;
	plot_setting.printOn();
	plot_setting.setOutputDims(x_id, y_id);
	plot_setting.plot_2D_interval_MATLAB("LTI_Abstraction_test", result);

	return 0;
}


















