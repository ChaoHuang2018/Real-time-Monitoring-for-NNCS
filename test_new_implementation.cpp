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

    string nn_name = "systems_with_networks/Benchmark1/nn_1_sigmoid";
    string act_name = "sigmoid";
    NeuralNetwork nn(nn_name, act_name);

    vector<string> state_vars;
    unsigned int numVars = 2;
    state_vars.push_back("x0");
    state_vars.push_back("x1");

    NNTaylor nn_taylor(nn);

    /* The follow is used to test nn parser */
    // vector<Interval> network_input_box;
    // network_input_box.push_back(Interval(0.4, 0.5));
    // network_input_box.push_back(Interval(0.4, 0.5));
    // nn_taylor.set_taylor_linear(state_vars, network_input_box);
    // // nn_taylor.set_range_by_IBP(network_input_box);
    // cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
    // cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;
    /* The above is used to test nn parser */

    Matrix<Interval> remainder(numVars, 1), stateSpace(numVars, 1), u(1, 1), d(2, 1);
    Interval x0(-2, 2), x1(-2, 2), control(-4.5, 4.5), disturbance(-0.1, 0.1);

    stateSpace[0][0] = x0;
    stateSpace[1][0] = x1;

    u[0][0] = control;
    d[0][0] = d[1][0] = disturbance;

    // Define the continuous dynamics.

    // set the current state
    Matrix<Interval> x_current(2, 1);
    x_current[0][0] = 0.754166;
    x_current[1][0] = 0.745562;

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
    // cout << "Initial guess of the reachable set:" << domain << endl;

    // get the linear taylor model and the corresponding remainder of the dynamics over the primary guess
    Matrix<Interval> linearization_remainder(2, 1);
    Matrix<string> dynamics_string(2, 1);
    Matrix<double> dynamics_coeff(2, 2);
    Matrix<double> dynamics_const_term(2, 1);
    dynamics_linear_taylor_benchmark1(dynamics_string, dynamics_coeff, dynamics_const_term, domain);
    remainder_linear_taylor_benchmark1(linearization_remainder, domain);
    // cout << "linear taylor expansion of dynamics: " << dynamics_string[0][0] << ", " << dynamics_string[1][0] << endl;
    // cout << "linear taylor Remainder of dynamics: " << linearization_remainder << endl;

    // get the linear taylor model and the corresponding remainder of the nn controller over the primary guess
    vector<Interval> network_input_box;
    network_input_box.push_back(domain[0][0]);
    network_input_box.push_back(domain[1][0]);
    nn_taylor.set_taylor_linear(state_vars, network_input_box);
    // nn_taylor.set_range_by_IBP(network_input_box);
    // cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
    // cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;

    Matrix<double> nn_coeff(2, 2);
    Matrix<double> nn_const_term(2, 1);
    nn_coeff[1][0] = nn_taylor.get_jacobian()[0];
    nn_coeff[1][1] = nn_taylor.get_jacobian()[1];
    nn_const_term[1][0] = nn_taylor.get_output();

    end = clock();
    cout << "Computition time for preparing: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

    // compute reachable set for LTI
    // Declaration of the state variables.
    unsigned int numVars = 3;

    int x0_id = stateVars.declareVar("x0");
    int x1_id = stateVars.declareVar("x1");
    int u_id = stateVars.declareVar("u");

    int domainDim = numVars + 1;

    Expression_AST<Real> deriv_x0(dynamics_string[0][0] + " + " + linearization_remainder[0][0].toString() + " + " + d[0][0].toString());     // do not forget the disturbance
    Expression_AST<Real> deriv_x1(dynamics_string[1][0] + " + " + linearization_remainder[1][0].toString() + " + u + " + d[1][0].toString()); // do not forget the disturbance and approximation error
    Expression_AST<Real> deriv_u("0");

    vector<Expression_AST<Real>> ode_rhs(numVars);
    ode_rhs[x0_id] = deriv_x0;
    ode_rhs[x1_id] = deriv_x1;
    ode_rhs[u_id] = deriv_u;

    Deterministic_Continuous_Dynamics dynamics(ode_rhs);

    // Specify the parameters for reachability computation.
    Computational_Setting setting;
    unsigned int order = 4;
    // stepsize and order for reachability analysis
    setting.setFixedStepsize(0.01, order);
    // time horizon for a single control step
    setting.setTime(0.02);
    // cutoff threshold
    setting.setCutoffThreshold(1e-10);
    // queue size for the symbolic remainder
    setting.setQueueSize(1000);
    // print out the steps
    setting.printOn();
    // remainder estimation
    Interval I(-0.01, 0.01);
    vector<Interval> remainder_estimation(numVars, I);
    setting.setRemainderEstimation(remainder_estimation);
    setting.printOff();
    setting.prepare();

    // define the initial set which is a box
    Interval init_x0(x_current[0][0]), init_x1(x_current[0][1]), init_u(0);
    vector<Interval> initialSet;
    initialSet.push_back(init_x0);
    initialSet.push_back(init_x1);
    initialSet.push_back(init_u);

    // translate the initial set to a flowpipe
    Flowpipe initial_set(initialSet);

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

    // perform 10 control steps
    for (int iter = 0; iter < 10; ++iter)
    {
        vector<Interval> box;
        initial_set.intEval(box, order, setting.tm_setting.cutoff_threshold);

        string strBox = "[" + box[0].toString() + "," + box[1].toString() + "]";

        string strExpU = nn_taylor.get_taylor_expression();

        double err = nn_taylor.get_taylor_remainder().sup();

        Expression_AST<Real> exp_u(strExpU);

        TaylorModel<Real> tm_u;
        exp_u.evaluate(tm_u, initial_set.tmvPre.tms, order, initial_set.domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

        tm_u.remainder.bloat(err);

        initial_set.tmvPre.tms[u_id] = tm_u;

        dynamics.reach(result, setting, initial_set, unsafeSet);

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
    }

    end = clock();
    cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

    return 0;
}