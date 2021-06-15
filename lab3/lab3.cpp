#include <iostream>
#include <cmath>
#include <cfloat>
#include <iomanip>

using namespace std;


const double NMAX = 30;
const double TOLX = DBL_EPSILON;
const double TOLF = DBL_EPSILON;


// pomocnicze funkcje do printowania
void printHeader() {
    cout << endl << "n   x_n1                    "
                    "en                        residuum" << endl;
}

void printHeaderSecant() {
    cout << endl << "n   x_n2                    "
                    "en                        residuum" << endl;
}

void printIteration(int n, double xn1, double en, double residuum) {
    cout << setw(2) << left << n << "  "
         << setprecision(17) << setw(24) << left << xn1
         << setprecision(17) << setw(24) << left << en  << "  "
         << setprecision(17) << setw(24) << left << residuum << endl;
}

void printTOLX(int n, double estimator) {
    cout << "przekroczono TOLX po " << n
         << " iteracjach - estymator = " << estimator << endl;
}

void printTOLF(int n, double residuum) {
    cout << "przekroczono TOLF po " << n
         << " iteracjach - residuum = " << residuum << endl;
}

void printTOLXAndTOLF(int n, double estimator, double residuum) {
    cout << "przekroczono TOLX oraz TOLF po " << n
         << " iteracjach - estymator = " << estimator
         << ", residuum = " << residuum << endl;
}

void printNMAX() {
    cout << "przekroczono maksymalną liczbę iteracji - " << NMAX << endl;
}


double f1(double x) {
    return sin(x/4.0) * sin(x/4.0) - x;
}

double f1_derivative(double x) {
    return 0.25 * sin(x/2.0) - 1;
}

double f2(double x) {
    return tan(2.0*x) - x - 1.0;
}

double f2_derivative(double x) {
    return 2.0 / (cos(2.0*x) * cos(2.0*x));

}


double phi1_Picard(double xn) {
    return sin(xn / 4.0) * sin(xn / 4.0);
}

double phi2_Picard(double xn) {
    return atan(xn + 1.0) / 2.0;
}

double phi2_Picard_divergent(double xn) {
    return 2.0 / (cos(2*xn) * cos(2*xn));
}

double phi1_Newton(double xn) {
    return xn - f1(xn) / f1_derivative(xn);
}

double phi2_Newton(double xn) {
    return xn - f2(xn) / f2_derivative(xn);
}

double phi1_Secant(double xn, double xn1) {
    return xn1 - f1(xn1) / ((f1(xn1) - f1(xn)) / (xn1 - xn));
}

double phi2_Secant(double xn, double xn1) {
    return xn1 - f2(xn1) / (f2(xn1 - f2(xn)) / (xn1 - xn));
}

double solvePhi(double (*phi)(double), double (*f)(double), double x0) {
    double xn = x0, xn1, en, residuum;

    printHeader();

    // arbitralne ograniczenie na liczbę iteracji
    for (int n = 1; n <= NMAX; n++) {
        xn1 = phi(xn);

        en = xn1 - xn;
        residuum = f(xn1);

        printIteration(n, xn1, en, residuum);

        // kryterium dokładności wyznaczenia xn1
        // oraz wiarygodności xn1 jako przybliżenia pierwiastków
        if (abs(en) <= TOLX && abs(residuum) <= TOLF) {
            printTOLXAndTOLF(n, en, residuum);
            return xn1;

        } else {
            // kryterium dokładności wyznaczenia xn
            if (abs(en) <= TOLX) {
                printTOLX(n, en);
                cout << "res: " << residuum << ", f(x): " << f(xn1) << endl;
                return xn1;
            }

            // kryterium wiarygodności xn jako przybliżenia pierwiastka
            if (abs(residuum) <= TOLF) {
                printTOLF(n, residuum);
                return xn1;
            }
        }

        xn = xn1;
    }

    printNMAX();
    return xn1;
}

double solveSecant(double (*f)(double), double x0, double x1) {
    double xn = x0, xn1 = x1, xn2, en, residuum;

    printHeaderSecant();

    // arbitralne ograniczenie na liczbę iteracji
    for (int n = 0; n < NMAX; n++) {
        xn2 = xn1 - f(xn1) / (((f(xn1) - f(xn)) / (xn1 - xn)));

        en = xn2 - xn1;
        residuum = f(xn2);

        printIteration(n, xn2, en, residuum);

        // kryterium dokładności wyznaczenia xn1
        // oraz wiarygodności xn1 jako przybliżenia pierwiastków
        if (abs(en) <= TOLX && abs(residuum) <= TOLF) {
            printTOLXAndTOLF(n, en, residuum);
            return xn2;

        } else {
            // kryterium dokładności wyznaczenia xn
            if (abs(en) <= TOLX) {
                printTOLX(n, en);
                return xn2;
            }

            // kryterium wiarygodności xn jako przybliżenia pierwiastka
            if (abs(residuum) <= TOLF) {
                printTOLF(n, residuum);
                return xn2;
            }
        }

        xn = xn1;
        xn1 = xn2;
    }

    printNMAX();
    return xn2;
}

double solveBisection(double (*f)(double), double a0, double b0) {
    double xn, an = a0, bn = b0, f_xn, f_an, f_bn, en, residuum;

    // sprawdzenie poprawności podanego przedziału [an, bn]
    double f_a0 = f(a0), f_b0 = f(b0);
    if ((f_a0 < 0.0 && f_b0 < 0.0) || (f_a0 >= 0.0 && f_b0 >= 0.0)) {
        cout << "podany przedział jest zły - f(a0) oraz f(b0) są tego samego znaku"
             << endl;
        return -1.0;
    }

    printHeader();

    // arbitralne ograniczenie na liczbę iteracji
    for (int n = 0; n < NMAX; n++) {
        xn = (an + bn) / 2.0;

        f_xn = f(xn);
        f_an = f(an);
        f_bn = f(bn);

        if ((f_an <= 0.0 && f_xn >= 0.0) || (f_an >= 0.0 && f_xn <= 0.0))
            bn = xn;
        else if ((f_xn <= 0.0 && f_bn >= 0.0) || (f_xn >= 0.0 && f_bn <= 0.0))
            an = xn;

        en = (bn - an) / 2.0;
        residuum = f_xn;

        printIteration(n, xn, en, residuum);

        // kryterium dokładności wyznaczenia xn
        // oraz wiarygodności xn jako przybliżenia pierwiastków
        if (abs(en) <= TOLX && abs(residuum) <= TOLF) {
            printTOLXAndTOLF(n, en, residuum);
            return xn;

        } else {
            // kryterium dokładności wyznaczenia xn
            if (abs(en) <= TOLX) {
                printTOLX(n, en);
                return xn;
            }

            // kryterium wiarygodności xn jako przybliżenia pierwiastka
            if (abs(residuum) <= TOLF) {
                printTOLF(n, residuum);
                return xn;
            }
        }
    }

    printNMAX();
    return xn;
}

//int main() {
//    cout << endl << "Metoda Picarda dla f1: ";
//    double root1_Picard = solvePhi(phi1_Picard, f1, 1.0);
//    cout << "x0 = " << root1_Picard << ", f1(x0) = " << f1(root1_Picard) << endl;
//
//    cout << endl << "Metoda Picarda dla f2: ";
//    double root2_Picard = solvePhi(phi2_Picard, f2, 1.0);
//    cout << "x0 = " << root2_Picard << ", f2(x0) = " << f2(root2_Picard) << endl;
//};

//int main() {
//    cout << endl << "Metoda Picarda dla f2: ";
//    double root2_Picard = solvePhi(phi2_Picard_divergent, f1, 1.0);
//    cout << "x0 = " << root2_Picard << ", f2(x0) = " << f2(root2_Picard) << endl;
//};

//int main() {
//    cout << endl << "Metoda Newtona dla f1: ";
//    double root1_Newton = solvePhi(phi1_Newton, f1, -1.0);
//    cout << "x0 = " << root1_Newton << ", f1(x0) = " << f1(root1_Newton) << endl;
//
//    cout << endl << "Metoda Newtona dla f2: ";
//    double root2_Newton = solvePhi(phi2_Newton, f2, -1.5);
//    cout << "x0 = " << root2_Newton << ", f2(x0) = " << f2(root2_Newton) << endl;
//};

//int main() {
//    cout << endl << "Metoda siecznych dla f1: ";
//    double root1_Secant = solveSecant(f1, 2.0, 1.0);
//    cout << "x0 = " << root1_Secant    << ", f1(x0) = " << f1(root1_Secant) << endl;
//
//    cout << endl << "Metoda siecznych dla f2: ";
//    double root2_Secant = solveSecant(f2, 2.0, 1.0);
//    cout << "x0 = " << root2_Secant    << ", f2(x0) = " << f2(root2_Secant) << endl;
//}

//int main() {
//    cout << endl << "Metoda bisekcji dla f1: ";
//    double root1_Bisection = solveBisection(f1, -1.5, 2.0);
//    cout << "x0 = " << root1_Bisection << ", f1(x0) = " << f1(root1_Bisection) << endl;
//
//    cout << endl << "Metoda bisekcji dla f2: ";
//    double root2_Bisection = solveBisection(f2, -0.6, 0.7);
//    cout << "x0 = " << root2_Bisection << ", f2(x0) = " << f2(root2_Bisection) << endl;
//}

int main() {
    cout << endl << "Metoda Picarda dla f1: ";
    double root1_Picard = solvePhi(phi1_Picard, f1, 1.0);
    cout << endl << "Metoda Picarda dla f2: ";
    double root2_Picard = solvePhi(phi2_Picard, f2, 1.0);
    cout << endl << "Metoda Newtona dla f1: ";
    double root1_Newton = solvePhi(phi1_Newton, f1, -1.0);
    cout << endl << "Metoda Newtona dla f2: ";
    double root2_Newton = solvePhi(phi2_Newton, f2, -1.0);
    cout << endl << "Metoda siecznych dla f1: ";
    double root1_Secant = solveSecant(f1, 2.0, 1.0);
    cout << endl << "Metoda siecznych dla f2: ";
    double root2_Secant = solveSecant(f2, 2.0, 1.0);
    cout << endl << "Metoda bisekcji dla f1: ";
    double root1_Bisection = solveBisection(f1, -1.5, 2.0);
    cout << endl << "Metoda bisekcji dla f2: ";
    double root2_Bisection = solveBisection(f2, -0.6, 0.7);

    cout << endl << "Dla funkcji f1(x) = sin(x/4)^2 - x: " << endl;
    cout << "Picard)    x0 = " << root1_Picard    << ", f1(x0) = " << f1(root1_Picard) << endl;
    cout << "Bisekji)   x0 = " << root1_Bisection << ", f1(x0) = " << f1(root1_Bisection) << endl;
    cout << "Newtona)   x0 = " << root1_Newton    << ", f1(x0) = " << f1(root1_Newton) << endl;
    cout << "Siecznych) x0 = " << root1_Secant    << ", f1(x0) = " << f1(root1_Secant) << endl;

    cout << endl << "Dla funkcji f2(x) = tan(2x) - x - 1: " << endl;
    cout << "Picard)    x0 = " << root2_Picard    << ", f2(x0) = " << f2(root2_Picard) << endl;
    cout << "Bisekji)   x0 = " << root2_Bisection << ", f2(x0) = " << f2(root2_Bisection) << endl;
    cout << "Newtona)   x0 = " << root2_Newton    << ", f2(x0) = " << f2(root2_Newton) << endl;
    cout << "Siecznych) x0 = " << root2_Secant    << ", f2(x0) = " << f2(root2_Secant) << endl;
};

