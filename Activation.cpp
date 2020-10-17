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
    long inf;
    long sup;

	mpfr_exp(inf, x.inf, MPFR_RNDD);
    mpfr_add(inf, 1L, inf, MPFR_RNDD);
    inf = mpfr_log(inf, inf, MPFR_RNDD);

	mpfr_exp(sup, x.sup, MPFR_RNDU);
	mpfr_add(inf, 1L, sup, MPFR_RNDU);
    sup = mpfr_log(sup, sup, MPFR_RNDU);

	Interval result(inf, sup);
    return result;
}

Interval Activation::softplus_de(Interval x) {
    long inf;
    long sup;

    mpfr_mul_si(inf, x.inf, -1L, MPFR_RNDD);
	mpfr_exp(inf, inf, MPFR_RNDD);
	mpfr_add(inf, 1L, inf, MPFR_RNDD);
    mpfr_div_d(inf, 1L, inf, MPFR_RNDD);

	mpfr_exp(sup, x.sup, MPFR_RNDU);
	mpfr_exp(sup, sup, MPFR_RNDU);
	mpfr_add(sup, 1L, sup, MPFR_RNDU);
    mpfr_div_d(sup, 1L, sup, MPFR_RNDU);

	Interval result(inf, sup);
    return result;
}

Interval Activation::softplus_de2(Interval x) {
    long inf;
    long sup;

    if x.inf

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
