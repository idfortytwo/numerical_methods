#include <iostream>
#include <vector.h>

using namespace std;


 void calculateEta(Vector &l, Vector &d, Vector &u);
 void calculateR(Vector &b, Vector &l, Vector &eta);
 Vector solve(Vector &eta, Vector &r, Vector &u);


int main() {
    Vector b = {31.0, 165.0/4.0, 917.0/30, 851.0/28.0, 3637.0/90, 332.0/11.0};
    cout << "Rozwiązujemy układ równań Ax = b używając algorytm Thomasa" << endl;
    b.print("Wektor b:");

    cout << endl << "Macierz trójdiagonalną A przedstawiamy jako trzy wektory l, d oraz u" << endl;
    Vector l = {1.0/3.0, 1.0/5.0, 1.0/7.0, 1.0/9.0, 1.0/11.0};
    Vector d = {10.0, 20.0, 30.0, 30.0, 20.0, 10};
    Vector u = {1.0/2.0, 1.0/4.0, 1.0/6.0, 1.0/8.0, 1.0/10.0};
    l.print("Wektor l:");
    d.print("Wektor d:");
    u.print("Wektor u:");


    cout << endl << "Najpierw przekształcamy macierz A - zamiast wektora d obliczamy nowy wektor eta" << endl;
    calculateEta(l, d, u);
    Vector eta = d;
    eta.print("Wektor eta: ");

    cout << endl << "Następnie z wektor b przekształcamy na r" << endl;
    calculateR(b, l, eta);
    Vector r = b;
    b.print("Wektor r: ");

    cout << endl << "Ostatecznie obliczamy wektor rozwiązań x" << endl;
    Vector x = solve(eta, r, u);
    x.print("Wektor x:");


    cout << endl << "Sprawdźmy wynik mnożenia Ax" << endl;
    Vector expected_b;
    expected_b[0] = d[0]*x[0] + u[0]*x[1];
    for (int i = 1; i < 6; i++) {
        expected_b[i] = l[i-1]*x[i-1] + d[i]*x[i] + u[i]*x[i+1];
    }

    expected_b.print("Wektor Ax:");
    cout << endl << "Otrzymaliśmy wyniki bardzo bliskie do wektora b";
}


void calculateEta(Vector &l, Vector &d, Vector &u) {
    int N =  d.N;
    for (int i = 1; i < N; i++)
        d[i] -= l[i-1] / d[i-1] * u[i-1];
}

void calculateR(Vector &b, Vector &l, Vector &eta) {
    int N =  b.N;
    for (int i = 1; i < N; i++)
        b[i] -= l[i-1] / eta[i-1] * b[i-1];
}

Vector solve(Vector &eta, Vector &r, Vector &u) {
    Vector x;
    int N =  x.N-1;

    x[N-1] = r[N-1] / eta[N-1];
    for (int i = N-2; i >= 0; i--)
        x[i] = (r[i] - u[i] * x[i+1]) / eta[i];

    return x;
}
