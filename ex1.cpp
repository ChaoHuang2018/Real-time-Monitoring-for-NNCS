#include "ex1.h"

using namespace flowstar;
using namespace std;

vector<Result_Info> run_ex1_taylor(string nn_name, string act_name, string trajectory_file_name)
{
    clock_t begin, end;

    NeuralNetwork nn(nn_name, act_name);
    NNTaylor nn_taylor(nn);
    // vector<double> nn_range;
    // for (int i = 0; i < 10000; i++)
    // {
    //     double rand_x0 = ((double)rand() / (RAND_MAX)) * 4 - 2;
    //     double rand_x1 = ((double)rand() / (RAND_MAX)) * 4 - 2;
    //     nn_range.push_back(nn.);
    // }

    Trajectories tr(nn.get_num_of_inputs(), trajectory_file_name);
    cout << tr.get_traces().size() << endl;

    // declare the number of state variables
    unsigned int numVars = 2;
    vector<string> state_vars;
    state_vars.push_back("x0");
    state_vars.push_back("x1");

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

    vector<Result_Info> result_list;
    // start checking the safety along each trace
    vector<vector<vector<double>>> traces = tr.get_traces();
    for (int i = 0; i < traces.size(); i++)
    {
        cout << "Processing Trace: " << i << endl;

        // computation time list
        vector<double> verification_time;
        // safe (0) or unsafe (1) or unknown (2) or fail (3)
        vector<int> verification_result;

        // for (int j = 0; j < traces[i].size(); j++)
        for (int j = 0; j < 1; j++)
        {
            cout << "Processing Step: " << j << endl;
            begin = clock();

            // set the current state
            Matrix<Interval> x_current(numVars, 1);
            for (int s = 0; s < numVars; s++)
            {
                cout << traces[i][j][s] << endl;
                Interval eps(traces[i][j][s], traces[i][j][s] + 0.0000001);
                // Interval eps(traces[i][j][s]);
                x_current[s][0] = eps;
            }

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
            cout << "Initial guess of the reachable set:" << domain << endl;
            stateSpace = domain;
            // abstract_domain_by_nn_range(domain, stateSpace, nn, x_current);
            domain[0][0] = Interval(7.659786642459603e-01, 8.624801728752246e-01);
            domain[1][0] = Interval(-1.024677427504299e-01, 7.195489429823004e-01);
            // cout << "x_current: " << x_current << endl;
            cout << "Initial guess of the reachable set after refinement:" << domain << endl;

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
            cout << "1111111111111111111111111" << endl;
            nn_taylor.set_taylor_linear(state_vars, network_input_box);
            // nn_taylor.set_range_by_IBP(network_input_box);
            // cout << "Linear Taylor Expression of nn controller: " << nn_taylor.get_taylor_expression() << endl;
            // cout << "Linear Taylor Remainder of nn controller: " << nn_taylor.get_taylor_remainder() << endl;
            cout << "222222222222222" << endl;

            Matrix<double> nn_coeff(numVars, numVars);
            Matrix<double> nn_const_term(numVars, 1);
            for (int s = 0; s < numVars; s++)
            {
                nn_coeff[1][s] = nn_taylor.get_jacobian()[s];
            }
            nn_const_term[1][0] = nn_taylor.get_output();
            // cout << "nn_coeff: " << nn_coeff << endl;
            // cout << "nn_const_term: " << nn_const_term << endl;

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
            cout << "Ux for taylor: " << BB << endl;
            // cout << "AA: " << AA << endl;
            // cout << "BB: " << BB << endl;
            Linear_Time_Invariant_Dynamics abstraction(AA, BB);

            // Specify the parameters for reachability computation.
            Computational_Setting setting;
            unsigned int order = 4;
            // stepsize and order for reachability analysis
            setting.setFixedStepsize(0.02, order);
            // time horizon. NOT a single control step
            setting.setTime(0.2);
            // cutoff threshold
            setting.setCutoffThreshold(1e-8);
            // print out the steps
            setting.printOn();
            setting.prepare();

            // define the initial set which is a box
            Interval init_x0(x_current[0][0]), init_x1(x_current[1][0]);
            // cout << "x_current[0][0]: " << x_current[0][0] << endl;
            // cout << "x_current[1][0]: " << x_current[1][0] << endl;
            vector<Interval> box;
            box.push_back(init_x0);
            box.push_back(init_x1);

            Flowpipe initialSet(box);

            // translate the initial set to a flowpipe
            // Flowpipe initial_set(initialSet);

            // unsafe set
            vector<Constraint> unsafeSet;
            Constraint constraint1("-x0 + 0.2"); // x0 >= 0.2,
            Constraint constraint2("x0 - 0.5");  // x0 <= 0.5,
            Constraint constraint3("-x1 - 0.2"); // x1 >= - 0.2,
            Constraint constraint4("x1 - 0.5");  // x1 <= 0.5,
            unsafeSet.push_back(constraint1);
            unsafeSet.push_back(constraint2);
            unsafeSet.push_back(constraint3);
            unsafeSet.push_back(constraint4);
            // Constraint constraint2("x0 + 0.1");
            // unsafeSet.push_back(constraint2);

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

            plot_setting.plot_2D_interval_MATLAB(string("nn_1_") + act_name + string("_step_") + to_string(j), result);

            end = clock();
            cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;
            verification_time.push_back((double)(end - begin) / CLOCKS_PER_SEC);
        }
        Result_Info r(verification_time, verification_result);
        result_list.push_back(r);
    }

    return result_list;
}

vector<Result_Info> run_ex1_nn_range(string nn_name, string act_name, string trajectory_file_name)
{
    clock_t begin, end;

    NeuralNetwork nn(nn_name, act_name);
    // vector<double> nn_range;
    // for (int i = 0; i < 10000; i++)
    // {
    //     double rand_x0 = ((double)rand() / (RAND_MAX)) * 4 - 2;
    //     double rand_x1 = ((double)rand() / (RAND_MAX)) * 4 - 2;
    //     nn_range.push_back(nn.);
    // }

    Trajectories tr(nn.get_num_of_inputs(), trajectory_file_name);
    cout << tr.get_traces().size() << endl;

    // declare the number of state variables
    unsigned int numVars = 2;
    vector<string> state_vars;
    state_vars.push_back("x0");
    state_vars.push_back("x1");

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

    vector<Result_Info> result_list;
    // start checking the safety along each trace
    vector<vector<vector<double>>> traces = tr.get_traces();
    for (int i = 0; i < traces.size(); i++)
    {
        cout << "Processing Trace: " << i << endl;

        // computation time list
        vector<double> verification_time;
        // safe (0) or unsafe (1) or unknown (2) or fail (3)
        vector<int> verification_result;

        // for (int j = 0; j < traces[i].size(); j++)
        for (int j = 0; j < 1; j++)
        {
            cout << "Processing Step: " << j << endl;
            begin = clock();

            NNTaylor nn_taylor(nn);
            // set the current state
            Matrix<Interval> x_current(numVars, 1);
            for (int s = 0; s < numVars; s++)
            {
                cout << traces[i][j][s] << endl;
                Interval eps(traces[i][j][s], traces[i][j][s] + 0.0000001);
                // Interval eps(traces[i][j][s]);
                x_current[s][0] = eps;
            }

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
            // cout << "Initial guess of the reachable set:" << domain << endl;

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
            nn_taylor.set_range_by_IBP(network_input_box);

            // C is defined by the nn_coeff
            Matrix<Real> C(numVars, numVars);
            // a is defined by nn_const_term
            Matrix<Real> a(numVars, 1);
            // Uu is computed based on the disturbance, sample value on the domain center, and the taylor remainder
            Matrix<Interval> Uu(numVars, 1);
            Uu[0][0] = Interval(0);
            Uu[1][0] = nn_taylor.get_range_by_IBP();

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
            cout << "Ux for output range: " << BB << endl;
            // cout << "AA: " << AA << endl;
            // cout << "BB: " << BB << endl;
            Linear_Time_Invariant_Dynamics abstraction(AA, BB);

            // Specify the parameters for reachability computation.
            Computational_Setting setting;
            unsigned int order = 4;
            // stepsize and order for reachability analysis
            setting.setFixedStepsize(0.02, order);
            // time horizon. NOT a single control step
            setting.setTime(0.2);
            // cutoff threshold
            setting.setCutoffThreshold(1e-8);
            // print out the steps
            setting.printOn();
            setting.prepare();

            // define the initial set which is a box
            Interval init_x0(x_current[0][0]), init_x1(x_current[1][0]);
            // cout << "x_current[0][0]: " << x_current[0][0] << endl;
            // cout << "x_current[1][0]: " << x_current[1][0] << endl;
            vector<Interval> box;
            box.push_back(init_x0);
            box.push_back(init_x1);

            Flowpipe initialSet(box);

            // translate the initial set to a flowpipe
            // Flowpipe initial_set(initialSet);

            // unsafe set
            vector<Constraint> unsafeSet;
            Constraint constraint1("-x0 + 0.2"); // x0 >= 0.2,
            Constraint constraint2("x0 - 0.5");  // x0 <= 0.5,
            Constraint constraint3("-x1 - 0.2"); // x1 >= - 0.2,
            Constraint constraint4("x1 - 0.5");  // x1 <= 0.5,
            unsafeSet.push_back(constraint1);
            unsafeSet.push_back(constraint2);
            unsafeSet.push_back(constraint3);
            unsafeSet.push_back(constraint4);
            // Constraint constraint2("x0 + 0.1");
            // unsafeSet.push_back(constraint2);

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

            plot_setting.plot_2D_interval_MATLAB(string("nn_1_") + act_name + string("_step_") + to_string(j), result);

            end = clock();
            cout << "Totoal computition time for one point: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;
            verification_time.push_back((double)(end - begin) / CLOCKS_PER_SEC);
        }
        Result_Info r(verification_time, verification_result);
        result_list.push_back(r);
    }

    return result_list;
}

vector<Result_Info> run_ex1_control_range(string nn_name, string act_name, string trajectory_file_name)
{
}

void abstract_domain_by_nn_range(Matrix<Interval> &domain, Matrix<Interval> stateSpace, NeuralNetwork nn, Matrix<Interval> x_current)
{
    int numVars = 2;
    NNTaylor nn_taylor(nn);

    Matrix<Interval> remainder(numVars, 1), u(1, 1), d(2, 1);
    Interval x0(-2, 2), x1(-2, 2), control(-4.5, 4.5), disturbance(-0.1, 0.1);

    u[0][0] = control;
    d[0][0] = d[1][0] = disturbance;

    // cout << "x_current: " << x_current << endl;
    // cout << "Initial guess of the reachable set:" << domain << endl;

    // get the linear taylor model and the corresponding remainder of the dynamics over the primary guess
    // A is defined by dynamics_coeff
    Matrix<Interval> linearization_remainder(numVars, 1);
    Matrix<string> dynamics_string(numVars, 1);
    Matrix<double> dynamics_coeff(numVars, numVars);
    Matrix<double> dynamics_const_term(numVars, 1);
    dynamics_linear_taylor_benchmark1(dynamics_string, dynamics_coeff, dynamics_const_term, stateSpace);
    remainder_linear_taylor_benchmark1(linearization_remainder, stateSpace);
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
    nn_taylor.set_range_by_IBP(network_input_box);

    // C is defined by the nn_coeff
    Matrix<Real> C(numVars, numVars);
    // a is defined by nn_const_term
    Matrix<Real> a(numVars, 1);
    // Uu is computed based on the disturbance, sample value on the domain center, and the taylor remainder
    Matrix<Interval> Uu(numVars, 1);
    Uu[0][0] = Interval(0);
    Uu[1][0] = nn_taylor.get_range_by_IBP();

    Matrix<Real> AA(numVars, numVars);
    AA = A + B * C;
    Matrix<Interval> BB(numVars, 1);
    BB = B * (Uu + a) + Ux;

    int x0_id = stateVars.declareVar("x0");
    int x1_id = stateVars.declareVar("x1");

    Linear_Time_Invariant_Dynamics abstraction(AA, BB);

    // Specify the parameters for reachability computation.
    Computational_Setting setting;
    unsigned int order = 4;
    // stepsize and order for reachability analysis
    setting.setFixedStepsize(0.02, order);
    // time horizon. NOT a single control step
    setting.setTime(0.2);
    // cutoff threshold
    setting.setCutoffThreshold(1e-8);
    // print out the steps
    setting.printOn();
    setting.prepare();

    // define the initial set which is a box
    Interval init_x0(x_current[0][0]), init_x1(x_current[1][0]);
    // cout << "x_current[0][0]: " << x_current[0][0] << endl;
    // cout << "x_current[1][0]: " << x_current[1][0] << endl;
    vector<Interval> box;
    box.push_back(init_x0);
    box.push_back(init_x1);

    Flowpipe initialSet(box);

    // unsafe set
    vector<Constraint> unsafeSet;
    Constraint constraint1("-x0 + 0.2"); // x0 >= 0.2,
    Constraint constraint2("x0 - 0.5");  // x0 <= 0.5,
    Constraint constraint3("-x1 - 0.2"); // x1 >= - 0.2,
    Constraint constraint4("x1 - 0.5");  // x1 <= 0.5,
    unsafeSet.push_back(constraint1);
    unsafeSet.push_back(constraint2);
    unsafeSet.push_back(constraint3);
    unsafeSet.push_back(constraint4);
    // Constraint constraint2("x0 + 0.1");
    // unsafeSet.push_back(constraint2);

    /*
            * The structure of the class Result_of_Reachability is defined as below:
            * nonlinear_flowpipes: the list of computed flowpipes
            * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
            * fp_end_of_time: the flowpipe at the time T
            */
    Result_of_Reachability result;

    abstraction.reach(result, setting, initialSet, unsafeSet);

    result.transformToTaylorModels(setting, initialSet);

    vector<vector<Interval>> flowpipe_ranges;

    std::list<TaylorModelVec<Real>>::const_iterator tmvIter = result.tmv_flowpipes.begin();
    std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();

    for (; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
    {
        std::vector<Interval> box;

        tmvIter->intEval(box, fpIter->domain);

        flowpipe_ranges.push_back(box);
    }

    vector<Interval> range_union = flowpipe_ranges[0];

    for (int i = 1; i < flowpipe_ranges.size(); ++i)
    {
        for (int j = 0; j < numVars; ++j)
        {
            range_union[j].hull_assign(flowpipe_ranges[i][j]);
        }
    }

    for (int j = 0; j < numVars; ++j)
    {
        domain[j][0] = range_union[j];
    }
}

int test_ex1()
{
    string nn_name = "systems_with_networks/Benchmark1/nn_12_tanh_origin";
    string act_name = "tanh";
    string trajectory_file_name = "systems_with_networks/Benchmark1/nn_12_tanh_origin.txt";

    // use taylor model
    vector<Result_Info> result_list_taylor = run_ex1_taylor(nn_name, act_name, trajectory_file_name);
    // cout << "result_size:" << result_list.size() << endl;
    // cout << "size:" << result_list[0].get_result().size() << endl;

    double safe_rate = 0;
    double average_time = 0;
    // for (int i = 0; i < result_list_taylor.size(); i++)
    // {
    //     safe_rate += result_list_taylor[i].get_safe_rate();
    //     average_time += result_list_taylor[i].get_average_time();
    // }

    // use nn output range
    // vector<Result_Info> result_list_nn_range = run_ex1_nn_range(nn_name, act_name, trajectory_file_name);
    // safe_rate = 0;
    // average_time = 0;
    // for (int i = 0; i < result_list_nn_range.size(); i++)
    // {
    //     safe_rate += result_list_nn_range[i].get_safe_rate();
    //     average_time += result_list_nn_range[i].get_average_time();
    // }

    // cout << "False positive rate of Taylor: " << 1 - safe_rate * 1.0 / result_list_taylor.size() << endl;
    // cout << "Average time of Taylor: " << average_time * 1.0 / result_list_taylor.size() << endl;
    // cout << "False positive rate of output range: " << 1 - safe_rate * 1.0 / result_list_nn_range.size() << endl;
    // cout << "Average time of output range: " << average_time * 1.0 / result_list_nn_range.size() << endl;

    return 0;

    // declare the number of variables
    // vector<Interval> intervals;

    // Interval I1(1, 2), I2(-2, 0), I3(0, 5), I4(6, 6);
    // intervals.push_back(I1);
    // intervals.push_back(I2);
    // intervals.push_back(I3);
    // intervals.push_back(I4);

    // Interval int_union = I1;

    // for (int i = 1; i < 4; ++i)
    // {
    //     int_union.hull_assign(intervals[i]);
    // }

    // cout << int_union << endl;

    // return 0;
}
