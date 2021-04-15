#include "flowstar-template/Continuous.h"

using namespace flowstar;
using namespace std;

int main(int argc, char *argv[])
{
    intervalNumPrecision = 200;

    string act = "sigmoid";
    Interval intv(2.979895446948051e+00, 2.982221967857086e+00);
    int d = stoi(argv[1]);
    int p = stoi(argv[2]);
    UnivariatePolynomial<Real> up = gen_bern_poly(act, intv, d);
    cout << up << endl;
    double err = gen_bern_err(act, intv, d);
    cout << "Naive: " << err << endl;
    double err_sample = gen_bern_err_by_sample(up, act, intv, p);
    cout << "sample: " << err_sample << endl;
    return 0;
}
