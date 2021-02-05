#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"
#include "domain_computation.h"
#include "dynamics_linearization.h"
#include "LTI_Abstraction.h"
#include "Trajectories.h"
#include "Result_Info.h"

using namespace flowstar;
using namespace std;

int run_ex1(string nn_name, string act_name, string trajectory_file_name)
{
    clock_t begin, end;
    begin = clock();

    NeuralNetwork nn(nn_name, act_name);

    Trajectories tr(nn.get_num_of_inputs(), trajectory_file_name);
    cout << tr.get_traces().size() << endl;

    // declare the number of state variables
    unsigned int numVars = 2;
    vector<string> state_vars;
    state_vars.push_back("x0");
    state_vars.push_back("x1");

    NNTaylor nn_taylor(nn);

    Matrix<Interval> remainder(numVars, 1), stateSpace(numVars, 1), u(1, 1), d(2, 1);
    Interval x0(-2, 2), x1(-2, 2), control(-4.5, 4.5), disturbance(-0.1, 0.1);

    stateSpace[0][0] = x0;
    stateSpace[1][0] = x1;

    u[0][0] = control;
    d[0][0] = d[1][0] = disturbance;

    Matrix<Interval> x_current(numVars, 1);
    Interval initx0(0.894960333436978);
    Interval initx1(0.830269906483071);
    x_current[0][0] = initx0;
    x_current[1][0] = initx1;

    // get the primary guess of the reachable set in 1 second by 2 iterations
    Matrix<Interval> domain(numVars, 1);
    remainderEval_benchmark1(remainder, stateSpace, u, d, 0.1);
    int N = 2;
    for (int s = 0; s < N; ++s)
    {
        stateSpace = domain;
        remainderEval_benchmark1(remainder, stateSpace, u, d, 0.1);
        linearizationDomainEval_benchmark1(domain, x_current, remainder, u, d, 0.1);
    }
    // cout << "x_current: " << x_current << endl;
    cout << "Initial guess of the reachable set:" << domain << endl;

    // get the linear taylor model and the corresponding remainder of the dynamics over the primary guess
    // A is defined by dynamics_coeff
    Matrix<Interval> linearization_remainder(numVars, 1);
    Matrix<string> dynamics_string(numVars, 1);
    Matrix<double> dynamics_coeff(numVars, numVars);
    Matrix<double> dynamics_const_term(numVars, 1);
    dynamics_linear_taylor_benchmark1(dynamics_string, dynamics_coeff, dynamics_const_term, domain);
    remainder_linear_taylor_benchmark1(linearization_remainder, domain);
    Matrix<Real> A(numVars, numVars);
    A = dynamics_coeff;

    // Ux is computed based on the disturbance, sample value on the domain center, and the taylor remainder
    Matrix<Interval> Ux(numVars, 1);
    for (int s = 0; s < numVars; s++)
    {
        Ux[s][0] = d[s][0] + dynamics_const_term[s][0] + linearization_remainder[s][0];
    }
    // cout << "linear taylor expansion of dynamics: " << dynamics_string[0][0] << ", " << dynamics_string[1][0] << endl;
    // cout << "linear taylor Remainder of dynamics: " << linearization_remainder << endl;

    // B is determined by the dynamics
    Matrix<Real> B(numVars, numVars);
    B[0][0] = 0;
    B[0][1] = 0;
    B[1][0] = 0;
    B[1][1] = 1;

    // get the linear taylor model and the corresponding remainder of the nn controller over the primary guess
    vector<Interval> network_input_box;
    for (int s = 0; s < numVars; s++)
    {
        network_input_box.push_back(domain[s][0]);
    }
    nn_taylor.set_taylor_linear(state_vars, network_input_box);
    // nn_taylor.set_range_by_IBP(network_input_box);
    // cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
    // cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;

    Matrix<double> nn_coeff(numVars, numVars);
    Matrix<double> nn_const_term(numVars, 1);
    for (int s = 0; s < numVars; s++)
    {
        nn_coeff[1][s] = nn_taylor.get_jacobian()[s];
    }
    nn_const_term[1][0] = nn_taylor.get_output();
    cout << "nn_coeff: " << nn_coeff << endl;
    cout << "nn_const_term: " << nn_const_term << endl;

    // C is defined by the nn_coeff
    Matrix<Real> C(numVars, numVars);
    C = nn_coeff;
    // a is defined by nn_const_term
    Matrix<Real> a(numVars, 1);
    a = nn_const_term;
    // Uu is computed based on the disturbance, sample value on the domain center, and the taylor remainder
    Matrix<Interval> Uu(numVars, 1);
    Uu[0][0] = Interval(0);
    Uu[1][0] = nn_taylor.get_taylor_remainder();

    end = clock();
    cout << "Computition time for preparing: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

    // compute reachable set for LTI
    // Declaration of the state variables.
    int x0_id = stateVars.declareVar("x0");
    int x1_id = stateVars.declareVar("x1");

    // define the abstraction
    // LTI_Abstraction abstraction(A, B, a, C, Ux, Uu);
    // cout << "A: " << A << endl;
    // cout << "B: " << B << endl;
    // cout << "a: " << a << endl;
    // cout << "C: " << C << endl;
    Matrix<Real> AA(numVars, numVars);
    AA = A + B * C;
    Matrix<Interval> BB(numVars, 1);
    BB = B * (Uu + a) + Ux;
    cout << "AA: " << AA << endl;
    cout << "BB: " << BB << endl;
    Linear_Time_Invariant_Dynamics abstraction(AA, BB);

    // Specify the parameters for reachability computation.
    Computational_Setting setting;
    unsigned int order = 4;
    // stepsize and order for reachability analysis
    setting.setFixedStepsize(0.05, order);
    // time horizon. NOT a single control step
    setting.setTime(0.15);
    // cutoff threshold
    setting.setCutoffThreshold(1e-8);
    // print out the steps
    setting.printOn();
    setting.prepare();

    // define the initial set which is a box
    Interval init_x0(x_current[0][0]), init_x1(x_current[1][0]);
    cout << "x_current[0][0]: " << x_current[0][0] << endl;
    cout << "x_current[1][0]: " << x_current[1][0] << endl;
    vector<Interval> box;
    box.push_back(init_x0);
    box.push_back(init_x1);

    Flowpipe initialSet(box);

    // translate the initial set to a flowpipe
    // Flowpipe initial_set(initialSet);

    // unsafe set
    vector<Constraint> unsafeSet;
    // Constraint constraint1("-x0 + 0.1"); // x0 >= 0.1,
    // Constraint constraint2("x0 - 0.5");  // x0 <= 0.5,
    // Constraint constraint3("-x1 - 0.5"); // x1 >= - 0.5,
    // Constraint constraint4("x1 - 0.5");  // x1 <= 0.5,
    // unsafeSet.push_back(constraint1);
    // unsafeSet.push_back(constraint2);
    // unsafeSet.push_back(constraint3);
    // unsafeSet.push_back(constraint4);
    Constraint constraint2("x0 + 0.1"); //
    unsafeSet.push_back(constraint2);

    /*
            * The structure of the class Result_of_Reachability is defined as below:
            * nonlinear_flowpipes: the list of computed flowpipes
            * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
            * fp_end_of_time: the flowpipe at the time T
            */
    Result_of_Reachability result;

    abstraction.reach(result, setting, initialSet, unsafeSet);

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

    // flowpipes should be translated to single Taylor model vectors before plotting
    result.transformToTaylorModels(setting, initialSet);

    Plot_Setting plot_setting;
    plot_setting.setOutputDims(x0_id, x1_id);
    int mkres = mkdir("./outputs", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (mkres < 0 && errno != EEXIST)
    {
        printf("Can not create the directory for images.\n");
        exit(1);
    }

    plot_setting.plot_2D_interval_MATLAB(string("nn_1_") + act_name + string("_step_0"), result);

    end = clock();
    cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;
    return 0;
}

int main()
{
    string nn_name = "systems_with_networks/Benchmark1/nn_1_sigmoid";
    string act_name = "sigmoid";
    string trajectory_file_name = "systems_with_networks/Benchmark1/nn_1_sigmoid.txt";

    return run_ex1(nn_name, act_name, trajectory_file_name);
}
