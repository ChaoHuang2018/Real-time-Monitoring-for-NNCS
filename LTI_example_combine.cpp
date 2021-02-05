#include "../flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	// declare the number of variables
	unsigned int numVars = 2;

	int x_id = stateVars.declareVar("x");
	int y_id = stateVars.declareVar("y");


	// define the matrix A
	Matrix<Real> A(numVars, numVars);

	A[0][0] = -2.458438292122803;
	A[0][1] = 1;
	A[1][0] = -4.676983857311751e-01;
	A[1][1] = -1.982884199522081e-01;

	// define the disturbance term u
	Matrix<UnivariateTaylorModel<Real> > B(numVars, 1);

	Interval u1(1.378755017343599e+00 , 1.588583383644269e+00), u2(-4.783868693373951e+00 , -2.879023190611354e+00);

//	Interval u1(1.4,1.4), u2;

	B[0][0] = u1;
	B[1][0] = u2;


	// define the dynamics
	Linear_Time_Invariant_Dynamics dynamics(A, B);


	// set the reachability parameters
	Computational_Setting setting;

	// set the stepsize and the order
	setting.setFixedStepsize(0.05, 4);

	// set the time horizon
	setting.setTime(0.15);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// print out the computation steps
	setting.printOn();


	// call this function when all of the parameters are defined
	setting.prepare();


	// define the initial set which is a box
	Interval init_x(0.814205174564876), init_y(0.735371011525061);

	vector<Interval> box(numVars);
	box[x_id] = init_x;
	box[y_id] = init_y;

	Flowpipe initialSet(box);


	// unsafe set
	vector<Constraint> unsafeSet;
	Constraint constraint("-x + 2.2"); // x >= 2.2, change the constraint to x >= 2.1 will produce UNKNOWN
//	unsafeSet.push_back(constraint);

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

	dynamics.reach(result, setting, initialSet, unsafeSet);

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


	// +++++++++++++++++++++++++++++++++++++++++++

	result.transformToTaylorModels(setting, initialSet);

	vector<vector<Interval> > flowpipe_ranges;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();

	for(; tmvIter != result.tmv_flowpipes.end() ; ++tmvIter, ++fpIter)
	{
		std::vector<Interval> box;

		tmvIter->intEval(box, fpIter->domain);

		flowpipe_ranges.push_back(box);
	}

	vector<Interval> range_union = flowpipe_ranges[0];

	for(int i=1; i<flowpipe_ranges.size(); ++i)
	{
		for(int j=0; j<numVars; ++j)
		{
			range_union[j].hull_assign(flowpipe_ranges[i][j]);
		}
	}

	for(int j=0; j<numVars; ++j)
	{
		cout << range_union[j] << endl;
	}

	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	// +++++++++++++++++++++++++++++++++++++++++++





	// The following code is not necessary if you do not need a plot

	// flowpipes should be translated to single Taylor model vectors before plotting
//	result.transformToTaylorModels(setting, initialSet);

	Plot_Setting plot_setting;
	plot_setting.printOn();
	plot_setting.setOutputDims(x_id, y_id);
	plot_setting.plot_2D_interval_MATLAB("LTI_example", result);

	return 0;
}


















