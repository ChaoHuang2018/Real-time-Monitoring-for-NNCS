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

}

Interval Activation::softplus_de(Interval x) {

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
