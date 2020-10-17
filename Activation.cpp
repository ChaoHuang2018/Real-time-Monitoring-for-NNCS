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
    double inf;
    double sup;

	mpfr_exp(inf, x.inf, MPFR_RNDD);
    mpfr_add(inf, 1, inf, MPFR_RNDD);
    inf = mpfr_log(inf, inf, MPFR_RNDD);

	mpfr_exp(sup, x.sup, MPFR_RNDU);
	sup = 1 + sup;
    sup = mpfr_log(sup, sup, MPFR_RNDD);

	Interval result(inf, sup);

    return result;
}

Interval Activation::softplus_de(Interval x) {
    double inf;
    double sup;

	mpfr_exp(inf, -x.inf, MPFR_RNDD);
	inf = 1 + inf;
    inf = mpfr_log(inf, inf, MPFR_RNDD);

	mpfr_exp(sup, x.sup, MPFR_RNDU);
	sup = 1 + sup;
    sup = mpfr_log(sup, sup, MPFR_RNDD);

	Interval result(inf, sup);

    return result;
}

Interval Activation::softplus_de2(Interval x) {

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
