#include <iostream>
#include <old2/vector_old.h>
#include <old2/matrix_old.h>

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

template <int N>
Vector<N> solveDiagonal(const Matrix<N, N> &M, const Vector<N> &y) {
    Vector<N> x;
    for (int i = 0; i < N; i++)
        x[i] = y[i] / M[i][i];
    return x;
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

void printHeader(int N) {
    cout << " n  " << setw(14*N) << "xn   " << "    residuum      estimator" << endl;
}

template <int N>
void printIteration(int n, const Vector<N> &xn, double residuum, double estimator) {
    cout << setw(2) << n << "   " << xn << " " << setw(14)
         << residuum << " " << setw(14) << estimator << endl;
}

template <int N>
Vector<N> solveJacobi(Matrix<N, N> &A, Vector<N> &b, Vector<N> &x0) {
    Vector<N> x = x0, xnext, y;
    Matrix<N, N> L, D, U;
    double residuum, estimator;

    LDU_decomposition(A, L, D, U);

    Matrix LU = L + U;

    printHeader(N);
    cout << " 0   " << x << endl;

    for (int n = 1; n <= NMAX; n++) {
        // D * xnext = -(L + U) * x + b
        Vector<N> D_xnext = b - (LU * x);

        xnext = solveDiagonal(D, D_xnext);

        residuum  = (A * xnext - b).normMax();
        estimator = (xnext - x).normMax();

        x = xnext;

        printIteration(n, x, residuum, estimator);

        if (residuum < TOLF && estimator < TOLX) {
            cout << "TOLF i TOLX po " << n << " iteracjach" << endl;
            return x;
        }
    }

    cout << "przekroczono NMAX" << endl;
    return x;
}

template <int N>
Vector<N> solveGaussSeidel(const Matrix<N, N> &A, const Vector<N> &b, const Vector<N> &x0) {
    Vector<N> x = x0, xnext, y;
    Matrix<N, N> L, D, U;
    double residuum, estimator;

    LDU_decomposition(A, L, D, U);

    Matrix LD = L + D;

    printHeader(N);
    cout << " 0   " << x << endl;

    for (int n = 1; n <= NMAX; n++) {
        // (L + D) * xnext = -U * x + b
        Vector<N> LD_xnext = b - U * x;

        xnext = solveLowerTriangular(LD, LD_xnext);

        residuum  = (A * xnext - b).normMax();
        estimator = (xnext - x).normMax();

        x = xnext;

        printIteration(n, x, residuum, estimator);

        if (residuum < TOLF && estimator < TOLX) {
            cout << "TOLF i TOLX po " << n << " iteracjach" << endl;
            return x;
        }
    }

    cout << "przekroczono NMAX" << endl;
    return x;
}

template <int N>
Vector<N> solveSOR(const Matrix<N, N> &A, const Vector<N> &b, const Vector<N> &x0, double omega) {
    Vector<N> x = x0, xnext, y;
    Matrix<N, N> L, D, U;
    double residuum, estimator;

    LDU_decomposition(A, L, D, U);

    // L + 1/omega*D
    Matrix L_omega_D =  L + (D * (1/omega));
    // (1 - 1/omega)*D + U
    Matrix omega_DU = (D * (1 - 1/omega)) + U;

    printHeader(N);
    cout << " 0   " << x << endl;

    for (int n = 1; n <= NMAX; n++) {
        // (L + 1/omega*D) * xnext = -[(1 - 1/omega)*D + U] * x + b
        Vector<N> L_omega_D_xnext  = b - (omega_DU * x);

        xnext = solveLowerTriangular(L_omega_D, L_omega_D_xnext);

        residuum  = (A * xnext - b).normMax();
        estimator = (xnext - x).normMax();

        x = xnext;

        printIteration(n, x, residuum, estimator);

        if (residuum < TOLF && estimator < TOLX) {
            cout << "TOLF i TOLX po " << n << " iteracjach" << endl;
            return x;
        }
    }

    cout << "przekroczono NMAX" << endl;
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

    cout << "Metoda Jacobi:" << endl;
    solveJacobi(A,b, x0);

    cout << endl << "Metoda Gaussa-Seidela:" << endl;
    solveGaussSeidel(A, b, x0);

    cout << endl << "Metoda SOR z parametrem 0.75:" << endl;
    solveSOR(A, b, x0, 0.75);

    cout << endl << "Metoda SOR z parametrem 0.5:" << endl;
    solveSOR(A, b, x0, 0.5);
}