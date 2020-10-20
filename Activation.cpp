#include "Activation.h"

using namespace flowstar;
using namespace std;

Activation::Activation()
{
}

Activation::Activation(string act, Interval in, Interval in_bound, string approach)
{
    string activation = act;
    Interval input = in;
    Interval input_bound = in_bound;

    Interval value = Interval();
    Interval de = Interval();
    Interval de2 = Interval();
    Interval output_range = Interval();
    Interval de_range = Interval();
    Interval de2_range = Interval();
    if (activation == "softplus")
    {
        if (approach != "taylor")
        {
            value = softplus(input);
            output_range = softplus(input_bound);
        }
        else
        {
            value = softplus(input);
            de = softplus_de(input);
            de2 = softplus_de2(input);
            output_range = softplus(input_bound);
            de_range = softplus_de(input_bound);
            de2_range = softplus_de2(input_bound);
        }
    }
    if (activation == "sigmoid")
    {
        if (approach != "taylor")
        {
            value = sigmoid(input);
            output_range = sigmoid(input_bound);
        }
        else
        {
            value = sigmoid(input);
            de = sigmoid_de(input);
            de2 = sigmoid_de2(input);
            output_range = sigmoid(input_bound);
            de_range = sigmoid_de(input_bound);
            de2_range = sigmoid_de2(input_bound);
        }
    }
    if (activation == "tanh")
    {
        if (approach != "taylor")
        {
            value = tanh(input);
            output_range = tanh(input_bound);
        }
        else
        {
            value = tanh(input);
            de = tanh_de(input);
            de2 = tanh_de2(input);
            output_range = tanh(input_bound);
            de_range = tanh_de(input_bound);
            de2_range = tanh_de2(input_bound);
        }
    }
    if (activation == "Affine")
    {
        if (approach != "taylor")
        {
            value = affine(input);
            output_range = affine(input_bound);
        }
        else
        {
            value = affine(input);
            de = affine_de(input);
            de2 = affine_de2(input);
            output_range = affine(input_bound);
            de_range = affine_de(input_bound);
            de2_range = affine_de2(input_bound);
        }
    }
}

string Activation::get_activation()
{
    return activation;
}

Interval Activation::get_input()
{
    return Interval(input);
}

Interval Activation::get_input_bound()
{
    return Interval(input_bound);
}

Interval Activation::get_value()
{
    return Interval(value);
}

Interval Activation::get_de()
{
    return Interval(de);
}

Interval Activation::get_de2()
{
    return Interval(de2);
}

Interval Activation::get_output_range()
{
    return Interval(output_range);
}

Interval Activation::get_de_range()
{
    return Interval(de_range);
}

Interval Activation::get_de2_range()
{
    return Interval(de2_range);
}

Real Activation::softplus(Real x)
{
    Real result(x);

    result.exp_assign();
    result = result + 1;
    result.log_assign();

    return result;
}

Real Activation::softplus_de(Real x)
{
    Real result(x);

    result = -result;
    result.exp_assign();
    result = 1 / (1 + result);

    return result;
}

Real Activation::softplus_de2(Real x)
{
    return softplus_de(x) * (1 - softplus_de(x));
}

Interval Activation::softplus(Interval x)
{
    Real inf(softplus(Real(x.inf())));
    Real sup(softplus(Real(x.sup())));

    Interval result(inf, sup);
    return result;
}

Interval Activation::softplus_de(Interval x)
{
    Real inf(softplus_de(Real(x.inf())));
    Real sup(softplus_de(Real(x.sup())));

    Interval result(inf, sup);
    return result;
}

Interval Activation::softplus_de2(Interval x)
{
    vector<Real> check_list;
    check_list.push_back(Real(x.inf()));
    check_list.push_back(Real(x.sup()));

    if ((x.inf() <= 0) && (x.sup() >= 0))
    {
        check_list.push_back(Real(0));
    }

    Real inf(10000);
    Real sup(-10000);
    for (int i = 0; i < check_list.size(); i++)
    {
        if (softplus_de2(check_list[i]) <= softplus_de2(inf))
        {
            inf = check_list[i];
        }
        if (softplus_de2(check_list[i]) >= softplus_de2(sup))
        {
            sup = check_list[i];
        }
    }

    Interval result((softplus_de2(inf)), (softplus_de2(sup)));
    return result;
}

Real Activation::tanh(Real x)
{
    Real temp1(x);
    temp1.exp_assign();
    Real temp2(-x);
    temp2.exp_assign();

    Real result((temp1 - temp2) / (temp1 + temp2));

    return result;
}

Real Activation::tanh_de(Real x)
{
    Real temp1(x);
    temp1 = tanh(temp1);
    temp1.pow_assign(2);

    Real result(1 - temp1);
    return result;
}

Real Activation::tanh_de2(Real x)
{
    Real temp1(x);
    Real temp2(x);
    temp1 = tanh(temp1);
    temp2 = tanh(temp2);
    temp1.pow_assign(2);

    Real result(-2 * temp1 * (1 - temp2));

    return result;
}

Interval Activation::tanh(Interval x)
{
    Real inf(tanh(Real(x.inf())));
    Real sup(tanh(Real(x.sup())));

    Interval result(inf, sup);
    return result;
}

Interval Activation::tanh_de(Interval x)
{
    vector<Real> check_list;
    check_list.push_back(Real(x.inf()));
    check_list.push_back(Real(x.sup()));

    if ((x.inf() <= 0) && (x.sup() >= 0))
    {
        check_list.push_back(Real(0));
    }

    Real inf(10000);
    Real sup(-10000);
    for (int i = 0; i < check_list.size(); i++)
    {
        if (tanh_de(check_list[i]) <= tanh_de(inf))
        {
            inf = check_list[i];
        }
        if (tanh_de(check_list[i]) >= tanh_de(sup))
        {
            sup = check_list[i];
        }
    }

    Interval result((tanh_de(inf)), (tanh_de(sup)));
    return result;
}

Interval Activation::tanh_de2(Interval x)
{
    vector<Real> check_list;
    check_list.push_back(Real(x.inf()));
    check_list.push_back(Real(x.sup()));

    double const_temp = pow(3, 1 / 3.) / 3.;

    if ((x.inf() <= -const_temp) && (x.sup() >= -const_temp))
    {
        check_list.push_back(Real(-const_temp));
    }
    if ((x.inf() <= const_temp) && (x.sup() >= const_temp))
    {
        check_list.push_back(Real(const_temp));
    }

    Real inf(10000);
    Real sup(-10000);
    for (int i = 0; i < check_list.size(); i++)
    {
        if (tanh_de2(check_list[i]) <= tanh_de2(inf))
        {
            inf = check_list[i];
        }
        if (tanh_de2(check_list[i]) >= tanh_de2(sup))
        {
            sup = check_list[i];
        }
    }

    Interval result((tanh_de(inf)), (tanh_de(sup)));
    return result;
}

Real Activation::sigmoid(Real x)
{
    Real temp1(x);
    temp1.exp_assign();

    Real result = Real();
    result = 1.0 - 1.0 / (1.0 + temp1);

    return result;
}

Real Activation::sigmoid_de(Real x)
{
    Real temp1(x);
    temp1 = sigmoid(temp1);

    Real result(temp1 * (1 - temp1));

    return result;
}

Real Activation::sigmoid_de2(Real x)
{
    Real temp1(x);
    temp1 = sigmoid(temp1);

    Real result(2 * temp1 * temp1 * temp1 - 3 * temp1 * temp1 + temp1);

    return result;
}

Interval Activation::sigmoid(Interval x)
{
    Real inf(sigmoid(Real(x.inf())));
    Real sup(sigmoid(Real(x.sup())));

    Interval result(inf, sup);
    return result;
}

Interval Activation::sigmoid_de(Interval x)
{
    vector<Real> check_list;
    check_list.push_back(Real(x.inf()));
    check_list.push_back(Real(x.sup()));

    if ((x.inf() <= 0) && (x.sup() >= 0))
    {
        check_list.push_back(Real(0));
    }

    Real inf(10000);
    Real sup(-10000);
    for (int i = 0; i < check_list.size(); i++)
    {
        if (sigmoid_de(check_list[i]) <= sigmoid_de(inf))
        {
            inf = check_list[i];
        }
        if (sigmoid_de(check_list[i]) >= sigmoid_de(sup))
        {
            sup = check_list[i];
        }
    }

    Interval result((sigmoid_de(inf)), (sigmoid_de(sup)));
    return result;
}

Interval Activation::sigmoid_de2(Interval x)
{
    vector<Real> check_list;
    check_list.push_back(Real(x.inf()));
    check_list.push_back(Real(x.sup()));

    double const_temp1 = log(2 - pow(3, 1 / 2.));
    double const_temp2 = log(2 + pow(3, 1 / 2.));

    if ((x.inf() <= const_temp1) && (x.sup() >= const_temp1))
    {
        check_list.push_back(Real(const_temp1));
    }
    if ((x.inf() <= const_temp2) && (x.sup() >= const_temp2))
    {
        check_list.push_back(Real(const_temp2));
    }

    Real inf(10000);
    Real sup(-10000);
    for (int i = 0; i < check_list.size(); i++)
    {
        if (tanh_de2(check_list[i]) <= tanh_de2(inf))
        {
            inf = check_list[i];
        }
        if (tanh_de2(check_list[i]) >= tanh_de2(sup))
        {
            sup = check_list[i];
        }
    }

    Interval result((tanh_de(inf)), (tanh_de(sup)));
    return result;
}

Real Activation::affine(Real x)
{
    return Real(x);
}

Real Activation::affine_de(Real x)
{
    return Real(1);
}

Real Activation::affine_de2(Real x)
{
    return Real(0);
}

Interval Activation::affine(Interval x)
{
    return Interval(x);
}

Interval Activation::affine_de(Interval x)
{
    return Interval(1, 1);
}

Interval Activation::affine_de2(Interval x)
{
    return Interval(0);
}
