#include "Activation.h"

using namespace flowstar;
using namespace std;


Activation::Activation(string act, Interval in, Interval in_bound) {
    activation = act;
    input = in;
    input_bound = in_bound;

    if (activation = "softplus") {

    }

}

Interval Activation::softplus(Interval x) {
    Real inf(x.inf());
    Real sup(x.sup());

	inf.exp_assign();
    inf.log_assign(1 + inf);

	sup.exp_assign();
    sup.log_assign(1 + inf);

	Interval result(*inf, *sup);
    return result;
}

Interval Activation::softplus_de(Interval x) {
    Real inf(x.inf());
    Real sup(x.sup());

    inf = - inf;
    inf.exp_assign();
    inf = 1/(1 + inf);

	sup = - sup;
    sup.exp_assign();
    sup = 1/(1 + sup);

	Interval result(inf, sup);
    return result;
}

Interval Activation::softplus_de2(Interval x) {
    Real inf(x.inf());
    Real sup(x.sup());


    if (x.inf <= 0) && (x.sup >= 0) {

    }

}


Interval Activation::tanh(Interval x) {

}

Interval Activation::tanh_de(Interval x) {

}

Interval Activation::tanh_de2(Interval x) {

}


Interval Activation::sigmoid(Interval x) {

}

Interval Activation::sigmoid_de(Interval x) {

}

Interval Activation::sigmoid_de2(Interval x) {

}


Interval Activation::affine(Interval x) {

}

Interval Activation::affine_de(Interval x) {

}

Interval Activation::affine_de2(Interval x) {

}
