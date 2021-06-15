#include <iostream>
#include <old2/vector_old.h>
#include <old2/matrix_old.h>
#include <fstream>

using namespace std;


#define NMAX 30
#define TOLF 1e-8
#define TOLX 1e-8

template <int N>
void LDU_decomposition(const Matrix<N, N> &A, Matrix<N, N> &L, Matrix<N, N> &D, Matrix<N, N> &U) {
    L = {{0}};
    D = {{0}};
    U = {{0}};

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double value = A[i][j];
            if (i > j)
                L[i][j] = value;
            else if (i == j)
                D[i][j] = value;
            else
                U[i][j] = value;
        }
    }
}
template <int N, int M>
Vector<N> solveLowerTriangular(const Matrix<N, M> &L, const Vector<N> &y) {
    Vector<N> x;
    for (int n = 0; n < N; n++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += L[n][j] * x[j];
        x[n] = (y[n] - sum) / L[n][n];
    }
    return x;
}

template <int N>
Vector<N> solveSOR(const Matrix<N, N> &A, const Vector<N> &b, const Vector<N> &x0, double omega, int &iterations) {
    Vector<N> x = x0, xnext, y;
    Matrix<N, N> L, D, U;
    double residuum, estimator;

    LDU_decomposition(A, L, D, U);

    // L + 1/omega*D
    Matrix L_omega_D =  L + (D * (1/omega));
    // (1 - 1/omega)*D + U
    Matrix omega_DU = (D * (1 - 1/omega)) + U;

    int n = 1;
    for (; n <= NMAX; n++) {
        // (L + 1/omega*D) * xnext = -[(1 - 1/omega)*D + U] * x + b
        Vector<N> L_omega_D_xnext  = b - (omega_DU * x);

        xnext = solveLowerTriangular(L_omega_D, L_omega_D_xnext);

        residuum  = (A * xnext - b).normMax();
        estimator = (xnext - x).normMax();

        x = xnext;

        if (residuum < TOLF && estimator < TOLX) {
            iterations = n;
            return x;
        }
    }

    iterations = n;
    return x;
}

int main() {
    Matrix<4, 4> A = {{
        100.0, -1.0, 2.0, -3.0,
        1.0, 200.0, -4.0, 5.0,
        -2.0, 4.0, 300.0, -6.0,
        3.0, -5.0, 6.0, 400.0
    }};

    Vector<4> b = {{116.0, -226.0, 912.0, -1174.0}};
    Vector<4> x0 = {{2.0, 2.0, 2.0, 2.0}};

    ofstream out;
    out.open("../lab7/sor.csv");

    for (double omega = 0.001; omega < 2.0; omega += 0.001) {
        int iterations;
        solveSOR(A, b, x0, omega, iterations);
        out << omega << "," << iterations << endl;
    }

    out.close();
}