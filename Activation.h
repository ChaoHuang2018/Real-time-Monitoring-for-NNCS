#include "flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;

class Activation
{

    //
protected:
    string activation;
    Interval input;
    Interval input_bound;

    Interval value;
    Interval de;
    Interval de2;
    Interval output_range;
    Interval de_range;
    Interval de2_range;

public:
    Activation();
    Activation(string act, Interval in, Interval in_bound, string approach = "taylor");

    string get_activation();
    Interval get_input();
    Interval get_input_bound();
    Interval get_value();
    Interval get_de();
    Interval get_de2();
    Interval get_output_range();
    Interval get_de_range();
    Interval get_de2_range();

    static Real softplus(Real x);
    static Real softplus_de(Real x);
    static Real softplus_de2(Real x);
    static Interval softplus(Interval x);
    static Interval softplus_de(Interval x);
    static Interval softplus_de2(Interval x);

    static Real tanh(Real x);
    static Real tanh_de(Real x);
    static Real tanh_de2(Real x);
    static Interval tanh(Interval x);
    static Interval tanh_de(Interval x);
    static Interval tanh_de2(Interval x);

    static Real sigmoid(Real x);
    static Real sigmoid_de(Real x);
    static Real sigmoid_de2(Real x);
    static Interval sigmoid(Interval x);
    static Interval sigmoid_de(Interval x);
    static Interval sigmoid_de2(Interval x);

    static Real affine(Real x);
    static Real affine_de(Real x);
    static Real affine_de2(Real x);
    static Interval affine(Interval x);
    static Interval affine_de(Interval x);
    static Interval affine_de2(Interval x);
}