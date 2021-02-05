#include "flowstar-template/Continuous.h"
#include "NNTaylor.h"
#include "domain_computation.h"
#include "dynamics_linearization.h"
#include "LTI_Abstraction.h"
#include "Trajectories.h"
#include "Result_Info.h"

using namespace flowstar;
using namespace std;

vector<Result_Info> run_ex1_taylor(string nn_name, string act_name, string trajectory_file_name);

vector<Result_Info> run_ex1_nn_range(string nn_name, string act_name, string trajectory_file_name);

vector<Result_Info> run_ex1_control_range(string nn_name, string act_name, string trajectory_file_name);

void abstract_domain_by_nn_range(Matrix<Interval> &domain, Matrix<Interval> stateSpace, NeuralNetwork nn, Matrix<Interval> x_current);

int test_ex1();
