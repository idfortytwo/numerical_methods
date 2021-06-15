#include <iostream>
#include <iomanip>
#include <cmath>
#include <float.h>
#include <fstream>

using namespace std;

template <class T>
double f(T x) {
    return sin(x);
}

template <class T>
double f_derivative(T x) {
    return cos(x);
}

template <class T>
T backward_difference_3p(T xn_prev_prev, T xn_prev, T xn, T h) {
    return (static_cast<T>(0.5)*f(xn_prev_prev) - static_cast<T>(2.0)*f(xn_prev) + static_cast<T>(1.5)*f(xn)) / h;
}

template <class T>
T backward_difference_2p(T xn_prev, T xn, T h) {
    return (f(xn) - f(xn_prev)) / h;
}

template <class T>
T central_difference(T xn_prev, T xn_next, T h) {
    return (f(xn_next) - f(xn_prev)) / (static_cast<T>(2.0)*h);
}

template <class T>
T forward_difference_2p(T xn, T xn_next, T h) {
    return (f(xn_next) - f(xn)) / h;
}

template <class T>
T forward_difference_3p(T xn, T xn_next, T xn_next_next, T h) {
    return (-static_cast<T>(1.5)*f(xn) + static_cast<T>(2.0)*f(xn_next) - static_cast<T>(0.5)*f(xn_next_next)) / h;
}

template <class T>
void calculate(string filename) {
    T limit;
    if (is_same<T, double>::value)
        limit = DBL_EPSILON;
    else if (is_same<T, float>::value)
        limit = FLT_EPSILON;
    else
        exit(1);

    T a = static_cast<T>(0.0), b = static_cast<T>(M_PI_4), c = static_cast<T>(M_PI_2);
    T h = static_cast<T>(0.1);

    T derivative_a = f_derivative(a);
    T derivative_b = f_derivative(b);
    T derivative_c = f_derivative(c);

    ofstream out;
    out.open(filename);

    out << "h,0 progresywna trzypunktowa,0 progresywna dwupunktowa,"
           "$\\dfrac{\\pi}{4}$ wsteczna trzypunktowa,"
           "$\\dfrac{\\pi}{4}$ wsteczna dwupunktowa,"
           "$\\dfrac{\\pi}{4}$ centralna dwupunktowa,"
           "$\\dfrac{\\pi}{4}$ progresywna dwupunktowa,"
           "$\\dfrac{\\pi}{4}$ progresywna trzypunktowa,"
           "$\\dfrac{\\pi}{2}$ wsteczna trzypunktowa,"
           "$\\dfrac{\\pi}{2}$ wsteczna dwupunktowa" << endl;

    while (h > limit) {
        char delim = ',';
        out << h << delim
            << abs(derivative_a - forward_difference_2p(a, a+h, h)) << delim
            << abs(derivative_a - forward_difference_3p(a, a+h, a+h+h, h)) << delim

            << abs(derivative_b - backward_difference_3p(b-h-h, b-h, b, h)) << delim
            << abs(derivative_b - backward_difference_2p(b-h, b, h)) << delim
            << abs(derivative_b - central_difference(b-h, b+h, h)) << delim
            << abs(derivative_b - forward_difference_2p(b, b+h, h)) << delim
            << abs(derivative_b - forward_difference_3p(b, b+h, b+h+h, h)) << delim

            << abs(derivative_c - backward_difference_3p(c-h-h, c-h, c, h)) << delim
            << abs(derivative_c - backward_difference_2p(c-h, c, h)) << endl;

        h /= static_cast<T>(1.2);
    }

    out.close();
}


int main() {
    calculate<double>("../lab8/double.csv");
    calculate<float>("../lab8/float.csv");
}