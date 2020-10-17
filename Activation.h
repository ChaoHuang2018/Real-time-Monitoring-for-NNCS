#include "flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;

class Activation{

    //
    protected:
        string activation;
        Interval input;
        Interval input_bound;

    public:
        Activation(string act, Interval in, Interval in_bound);

        static Interval softplus(Interval x);
        static Interval softplus_de(Interval x);
        static Interval softplus_de2(Interval x);

        static Interval tanh(Interval x);
        static Interval tanh_de(Interval x);
        static Interval tanh_de2(Interval x);

        static Interval sigmoid(Interval x);
        static Interval sigmoid_de(Interval x);
        static Interval sigmoid_de2(Interval x);

        static Interval affine(Interval x);
        static Interval affine_de(Interval x);
        static Interval affine_de2(Interval x);

}