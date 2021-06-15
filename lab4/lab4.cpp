#pragma ide diagnostic ignored "modernize-use-auto"
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cfloat>

using namespace std;

const int NMAX = 30;
const double TOLX = DBL_EPSILON;
const double TOLF = DBL_EPSILON;


double f1(vector<double> xn) {
    double x = xn[0], y = xn[1], z = xn[2];
    return x*x + y*y + z*z - 2.0;
}

double f2(vector<double> xn) {
    double x = xn[0], y = xn[1];
    return x*x + y*y - 1.0;
}

double f3(vector<double> xn) {
    double x = xn[0], y = xn[1];
    return x*x - y;
}


double delta1(vector<double> xn) {
    double x = xn[0], y = xn[1];
    return (2.0*x*x*y + x*x - y*y - 1.0) / (2.0*x*(2.0*y + 1.0));
}

double delta2(vector<double> xn) {
    double y = xn[1];
    return (y*y + y - 1.0) / (2.0*y + 1.0);
}

double delta3(vector<double> xn) {
    double z = xn[2];
    return z/2.0 - 1.0/(2.0*z);
}

// pomocnicze funkcje do printowania
void printHeader(int vector_size) {
    cout << "n   ";
    for (int i = 0; i < vector_size; i++)
        cout << "xn_" << setw(17) << left << i << " ";
    cout << setw(26) << "max estymator" << "max residuum" << endl;
}

void printIteration(int vector_size, const vector<double> &xn1, double max_residuum, double max_estimator, int n) {
    cout << setw(2) << left << n << "  ";
    for (int i = 0; i < vector_size; i++)
        cout << setw(19) << setprecision(17) << left << xn1[i] << "  ";
    cout << setw(24) << max_estimator << "  " << setw(24) << max_residuum << endl;
}

void printTOLX(int n, double estimator) {
    cout << endl << "przekroczono TOLX po " << n
         << " iteracjach - max estymator = " << estimator << endl << endl;
}

void printTOLF(int n, double residuum) {
    cout << endl << "przekroczono TOLF po " << n
         << " iteracjach - max residuum = " << residuum << endl << endl;
}

void printTOLXAndTOLF(int n, double estimator, double residuum) {
    cout << endl << "przekroczono TOLX oraz TOLF po " << n
         << " iteracjach - max estymator = " << estimator
         << ", max residuum = " << residuum << endl << endl;
}

void printNMAX() {
    cout << "przekroczono maksymalną liczbę iteracji - " << NMAX << endl << endl;
}

vector<double> solve_system(vector<double> x0, vector<double (*)(vector<double>)> delta,
                            vector<double (*)(vector<double>)> f, int vector_size) {

    // wektor obliczonych miejsc zerowych układu równań f
    vector<double> xn1(vector_size);
    vector<double> xn = x0;

    vector<double> residuum(vector_size);
    vector<double> estimator(vector_size);
    double max_residuum;
    double max_estimator;

    printHeader(vector_size);

    // arbitralne ograniczenie na liczbę iteracji
    for (int n = 1; n <= NMAX; n++) {
        // iterowanie po elementach wektora xn
        for (int i = 0; i < vector_size; i++) {
            xn1[i] = xn[i] - delta[i](xn);

            estimator[i] = abs(xn1[i] - xn[i]);
            residuum[i]  = abs(f[i](xn1));
        }

        // wyznaczanie największego estymatora i residuum
        max_estimator = *max_element(estimator.begin(), estimator.end());
        max_residuum  = *max_element(residuum.begin(), residuum.end());

        printIteration(vector_size, xn1, max_residuum, max_estimator, n);

        // kryterium dokładności wyznaczenia xn1
        // oraz wiarygodności xn1 jako przybliżenia pierwiastków
        if (max_estimator <= TOLX && max_residuum <= TOLF) {
            printTOLXAndTOLF(n, max_estimator, max_residuum);
            return xn1;

        } else {
            // kryterium dokładności wyznaczenia xn1
            if (max_estimator <= TOLX) {
                printTOLX(n, max_estimator);
                return xn1;
            }

            // kryterium wiarygodności xn1 jako przybliżenia pierwiastków
            if (max_residuum <= TOLF) {
                printTOLF(n, max_residuum);
                return xn1;
            }
        }

        xn = xn1;
    }

    printNMAX();
    return xn1;
}


int main() {
    int vector_size = 3;

    vector<double> xn = {7.5, 5.0, 2.5};
    vector<double (*)(vector<double>)> delta = {delta1, delta2, delta3};
    vector<double (*)(vector<double>)> f = {f1, f2, f3};

    vector<double> roots = solve_system(xn, delta, f, vector_size);

    for (int i = 0; i < vector_size; i++)
        cout << "x_" << i << " = " << roots[i] << endl;
}