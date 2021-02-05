#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"
#include "domain_computation.h"
#include "dynamics_linearization.h"
#include "LTI_Abstraction.h"
#include "Trajectories.h"
#include "Result_Info.h"

using namespace flowstar;
using namespace std;

vector<Result_Info> run_ex6(string nn_name, string act_name, string trajectory_file_name)
{
    clock_t begin, end;
    begin = clock();

    NeuralNetwork nn(nn_name, act_name);
    Trajectories tr(nn.get_num_of_inputs(), trajectory_file_name);

    // declare the number of state variables
    unsigned int numVars = 4;
    vector<string> state_vars;
    state_vars.push_back("x0");
    state_vars.push_back("x1");
    state_vars.push_back("x2");
    state_vars.push_back("x3");

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

    Matrix<Interval> remainder(numVars, 1), stateSpace(numVars, 1), u(1, 1), d(numVars, 1);
    Interval x0(-2, 2), x1(-2, 2), x2(-2, 2), x3(-2, 2), control(-1, 1), disturbance(-0.1, 0.1);

    stateSpace[0][0] = x0;
    stateSpace[1][0] = x1;
    stateSpace[2][0] = x2;
    stateSpace[3][0] = x3;

    u[0][0] = control;
    d[0][0] = d[1][0] = d[2][0] = d[3][0] = disturbance;

    vector<Result_Info> result_list;
    // start checking the safety along each trace
    vector<vector<vector<double>>> traces = tr.get_traces();
    for (int i = 1; i < traces.size(); i++)
    {
        // computation time list
        vector<double> verification_time;
        // safe (0) or unsafe (1) or unknown (2) or fail (3)
        vector<int> verification_result;

        for (int j = 1; j < traces[i].size(); j++)
        {
            // set the current state
            Matrix<Interval> x_current(numVars, 1);
            for (int s = 0; s < numVars; s++)
            {
                x_current[s][0] = traces[i][j][s];
            }

            // get the primary guess of the reachable set in 1 second by 2 iterations
            Matrix<Interval> domain(numVars, 1);
            remainderEval_benchmark5(remainder, stateSpace, u, d, 1);
            int N = 2;
            for (int i = 0; i < N; ++i)
            {
                stateSpace = domain;
                remainderEval_benchmark6(remainder, stateSpace, u, d, 1);
                linearizationDomainEval_benchmark6(domain, x_current, remainder, u, d, 1);
            }
            // cout << "Initial guess of the reachable set:" << domain << endl;

            // get the linear taylor model and the corresponding remainder of the dynamics over the primary guess
            // A is defined by dynamics_coeff
            Matrix<Interval> linearization_remainder(numVars, 1);
            Matrix<string> dynamics_string(numVars, 1);
            Matrix<double> dynamics_coeff(numVars, numVars);
            Matrix<double> dynamics_const_term(numVars, 1);
            dynamics_linear_taylor_benchmark5(dynamics_string, dynamics_coeff, dynamics_const_term, domain);
            remainder_linear_taylor_benchmark5(linearization_remainder, domain);
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
            B[3][3] = 1;

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
                nn_coeff[3][s] = nn_taylor.get_jacobian()[s];
            }
            nn_const_term[2][0] = nn_taylor.get_output();

            // C is defined by the nn_coeff
            Matrix<Real> C(numVars, numVars);
            C = nn_coeff;
            // a is defined by nn_const_term
            Matrix<Real> a(numVars, 1);
            a = nn_const_term;
            // Uu is the taylor remainder
            Matrix<Interval> Uu(numVars, 1);
            Uu[0][0] = Interval(0);
            Uu[1][0] = Interval(0);
            Uu[2][0] = nn_taylor.get_taylor_remainder();

            end = clock();
            cout << "Computition time for preparing: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

            // compute reachable set for LTI
            // Declaration of the state variables.
            int x0_id = stateVars.declareVar("x0");
            int x1_id = stateVars.declareVar("x1");
            int x2_id = stateVars.declareVar("x2");

            // define the abstraction
            LTI_Abstraction abstraction(A, B, a, C, Ux, Uu);

            // Specify the parameters for reachability computation.
            Computational_Setting setting;
            unsigned int order = 4;
            // stepsize and order for reachability analysis
            setting.setFixedStepsize(0.2, order);
            // time horizon. NOT a single control step
            setting.setTime(1.0);
            // cutoff threshold
            setting.setCutoffThreshold(1e-8);
            // print out the steps
            setting.printOn();
            setting.prepare();

            // define the initial set which is a box
            Interval init_x0(x_current[0][0]), init_x1(x_current[0][1]), init_x2(x_current[0][2]), init_u(0);
            vector<Interval> initialSet;
            initialSet.push_back(init_x0);
            initialSet.push_back(init_x1);
            initialSet.push_back(init_x2);
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

            abstraction.reach(result, setting, initialSet, unsafeSet);

            switch (result.status)
            {
            case COMPLETED_SAFE:
                printf("Safe\n");
                verification_result.push_back(0);
                break;
            case COMPLETED_UNSAFE:
                printf("Unsafe\n");
                verification_result.push_back(1);
                break;
            case COMPLETED_UNKNOWN:
                printf("Unknown\n");
                verification_result.push_back(2);
                break;
            default: // never happen to linear systems
                printf("Fail to compute flowpipes.\n");
                verification_result.push_back(3);
            }

            end = clock();
            cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;
            verification_time.push_back((double)(end - begin) / CLOCKS_PER_SEC);
        }
        Result_Info r(verification_time, verification_result);
        result_list.push_back(r);
    }

    return result_list;
}

int test_ex6()
{
    string nn_name = "nn/nn_6_sigmoid";
    string act_name = "sigmoid";
    string trajectory_file_name = "trajectories/nn_6_sigmoid";

    vector<Result_Info> result_list = run_ex6(nn_name, act_name, trajectory_file_name);

    return 0;
}
