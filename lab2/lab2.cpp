#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;


double f(double x) {
    return (1 - exp(-x)) / x;
}

double recalculated(double x, double exact) {
    double result = 0, iter_part = 1;
    int sign = 1;

    for(int n = 1; abs(exact - result) > DBL_EPSILON; n++) {
        result += sign * iter_part;

        iter_part *= x / (n + 1);
        if (iter_part == 0)
            break;

        sign *= -1;
    }

    return result;
}


int main() {
    ifstream plik;
    ofstream output;

    plik.open("dane.txt");
    output.open("output.dat");

    double log10x, x, exact, calculated;
    while(plik.good()) {
        plik >> log10x >> x >> exact;

        calculated = log10x < -0.5 ? recalculated(x, exact) : f(x);

        output << log10x << ' ' << log10(abs((exact - calculated) / exact)) << endl;
    }
}
