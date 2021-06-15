#include <cmath>
#include <iostream>
#include <fstream>
#include <vector.h>

using namespace std;


double p(double x) { return  1.0; }
double q(double x) { return  0.0; }
double r(double x) { return -4.0; }
double s(double x) { return -x;   }

double Alpha = 0.0, Beta = 1.0, Gamma  = -1.0,
       Phi   = 0.0, Psi  = 1.0, Theta  =  0.0;
double left_bound = 0.0, right_bound = 1.0;


double f_exact(double x) {
    return (exp(2-2*x) - 4*exp(4-2*x) + 4*exp(2*x) - exp(2+2*x) - x + x*exp(4))
           / (4 - 4*exp(4));
}


void transformThomas(Vector& u, Vector& d, Vector& l, Vector &b) {
    for (int i = 1; i < u.N; i++) {
        d[i] -= l[i-1] / d[i-1] * u[i-1];
        b[i] -= l[i-1] / d[i-1] * b[i-1];
    }
}

Vector solveThomas(Vector &u, Vector &d, Vector &b) {
    int n = u.N;
    Vector x(n, 0);

    x[n-1] = b[n-1] / d[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = (b[i] - u[i] * x[i+1]) / d[i];
    }

    return x;
}


double conventional(double h, bool save=false) {
    int n = (int) ceil((right_bound - left_bound) / h) + 1;
    Vector u(n, 0), d(n, 0), l(n, 0), b(n, 0), y;

    double x_0 = left_bound, x_i, y_i, h_square = h * h;
    double maxError = 0.0, currentError;

    d[0] = Beta - Alpha / h;
    u[0] = Alpha / h;
    b[0] = -Gamma;

    for (int i = 1; i < n - 1; i++) {
        x_i = x_0 + h * i;
        h_square = h * h;

        l[i-1] = p(x_i) / h_square - q(x_i) / (2.0 * h);
        d[i]   = r(x_i) - 2.0 * p(x_i) / h_square;
        u[i]   = p(x_i) / h_square + q(x_i) / (2.0 * h);
        b[i]   = -s(x_i);
    }
    l[n-2] = -Phi / h;
    d[n-1] =  Phi / h + Psi;
    b[n-1] = -Theta;

    transformThomas(u, d, l, b);
    y = solveThomas(u, d, b);

    ofstream out;
    if (save) {
        out.open("../lab9/conventional.csv");
        out << "h,conventional" << endl;
    }

    for (int i = 0; i < n; i++) {
        x_i = x_0 + h * i;
        y_i = y[i];

        currentError = abs(y_i - f_exact(x_i));
        if (maxError < currentError)
            maxError = currentError;

        if (save) {
            out << x_i << "," << y_i << endl;
        }
    }

    out.close();
    return maxError;
}

double numerov(double h, bool save=false) {
    int n = (int) ceil((right_bound - left_bound) / h) + 1;
    Vector u(n, 0), d(n, 0), l(n, 0), b(n, 0), y;

    double x_0 = left_bound, x_i, y_i, h_square = h * h;
    double maxError = 0.0, currentError;

    d[0] = Beta - Alpha / h;
    u[0] = Alpha / h;
    b[0] = -Gamma;

    for (int i = 1; i < n - 1; i++) {
        x_i = x_0 + h * i;

        l[i-1] = p(x_i) / h_square - q(x_i) / (2.0 * h) + r(x_i) / 12.0;
        d[i]   = 10.0 / 12.0 * r(x_i) - 2.0 * p(x_i) / h_square;
        u[i]   = p(x_i) / h_square - q(x_i) / (2.0 * h) + r(x_i) / 12.0;
        b[i]   = -s(x_i);
    }

    l[n-2] = -Phi / h;
    d[n-1] =  Phi / h + Psi;
    b[n-1] = -Theta;

    transformThomas(u, d, l, b);
    y = solveThomas(u, d, b);

    ofstream out;
    if (save) {
        out.open("../lab9/numerov.csv");
        out << "h,Numerov" << endl;
    }

    for (int i = 0; i < n; i++) {
        x_i = x_0 + h * i;
        y_i = y[i];

        currentError = abs(y_i - f_exact(x_i));
        if (maxError < currentError)
            maxError = currentError;

        if (save) {
            out << x_i << "," << y_i << endl;
        }
    }

    out.close();
    return maxError;
}

double order(double h_i, double h_j, double y_i, double y_j) {
    return (y_i - y_j) / (h_i - h_j);
}

int main() {
    double h = 0.01;
    conventional(h, true);
    numerov(h, true);

    double conventionalError, numerovError;
    double h_0, conventional_0, numerov_0,
           h_5, conventional_5, numerov_5;

    ofstream errors;
    errors.open("../lab9/errors.csv");
    errors << "h,Numerov,conventional" << endl;

    int count = 0;
    h = 0.125;
    while (h > 1e-6) {
        conventionalError = conventional(h);
        numerovError = numerov(h);
        errors << log10(h) << ","
               << log10(conventionalError) << ","
               << log10(numerovError) << endl;

        if (count == 0) {
            h_0 = log10(h);
            conventional_0 = log10(conventionalError);
            numerov_0 = log10(numerovError);
        } else if (count == 5) {
            h_5 = log10(h);
            conventional_5 = log10(conventionalError);
            numerov_5 = log10(numerovError);
        }

        count++;

        h /= 2.0;
    }

    cout << "rząd metody z dyskretyzacją trzypunktową konwencjonalną: "
         << order(h_0, h_5, conventional_0, conventional_5) << endl;
    cout << "rząd metody z dyskretyzacją Numerowa: "
         << order(h_0, h_5, numerov_0, numerov_5) << endl;

    errors.close();
}