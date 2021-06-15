#include "old2/vector_old.h"
#include "old2/matrix_old.h"

using namespace std;


template <int N, int M> int getPartialPivot(const Matrix<N, M> &A, int from);
template <int N, int M> pair<Matrix<N, M>, Matrix<N, M>> LU_decomposition(const Matrix<N, M> &A, Vector<N> &indexes);
template <int N> void rearrangeIndexes(Vector<N> &vect, const Vector<N> &indexes);
template <int N, int M> Vector<N> solveL(Matrix<N, M> &L, Vector<N> &b);
template <int N, int M> Vector<N> solveU(Matrix<N, M> &U, Vector<N> &y);
template <int N, int M> Vector<N> solveAxb(const Matrix<N, M> &A, Vector<N> &b);


int main() {
    Matrix<3, 2> l;

    Matrix<4, 4> A = {{
        1.0, -20.0, 30.0, -4.0,
        2.0, -40.0, -6.0, 50.0,
        9.0, -180.0, 11.0, -12.0,
        -16.0, 15.0, -140.0, 13.0
    }};

    Vector<4> b = {{35.0, 104.0, -366.0, -354.0}};

    Vector<4> x = solveAxb(A, b);

    cout << endl << "Na koniec sprawdźmy czy otrzymany wektor x faktycznie zawiera rozwiązania do pierwotnego równania Ax = b" << endl;

//    Vector<4> expected_b = matrixTimesVector(A, x);
//    expected_b.print();
//    cout << "Po wymnożeniu macierzy A i obliczonego wektora x ewidentnie otrzymaliśmy wektor b" << endl;
}

template <int N, int M>
Vector<N> solveAxb(const Matrix<N, M> &A, Vector<N> &b) {
    Vector<N> x;
    Vector<N> indexes;
    for (int i = 0; i < 4; i++)
        indexes[i] = i;

    cout << "Rozwiązujemy układ równań w postaci Ax = b metodą dekompozycji" << endl;
    A.print("Macierz A: ", 5);
    b.print("Wektor b: ", 5);


    cout << endl << endl << "Zaczynamy dekompozycję LU" << endl;
    auto [L, U] = LU_decomposition(A, indexes);

    cout << "Po dekompozycji otrzymaliśmy macierze L oraz U" << endl;
    L.print("Macierz L: ");
    U.print("Macierz U: ");
    indexes.print("Nowa kolejność wierszy w macierzach L, U oraz wektorze b,\nktóra powstała w wyniku zamian wierszy po wyborze częściowym:", 2);

    rearrangeIndexes(b, indexes);
    b.print("Wektor b po zmianie kolejności wierszy: ", 8);


    cout << endl << "Rozwiązujemy układ równań Ly = b, gdzie L - macierz dolnotrójkątna" << endl;
    Vector<N> y = solveL(L, b);
    y.print("Wektor y: ", 8);


    cout << endl << "Rozwiązujemy układ równań Ux = y, gdzie U - macierz górnotrójkątna" << endl;
    x = solveU(U, y);
    x.print("Wektor x:", 8);
    cout << "Znaleźliśmy rozwiązanie dla układu równań Ax = b" << endl;

    return x;
}

template <int N, int M>
int getPartialPivot(const Matrix<N, M> &A, int from) {
    int pivot = from;
    double max = A[from][from];

    for (int i = from; i < N; i++) {
        double value = std::abs(A[i][from]);
        if (value > max) {
            max = value;
            pivot = i;
        }
    }

    return pivot;
}

template <int N, int M>
pair<Matrix<N, M>, Matrix<N, M>> LU_decomposition(const Matrix<N, M> &A, Vector<N> &indexes) {
    Matrix<N, M> L;
    Matrix<N, M> U = A;

    for (int k = 0; k < N-1; k++) {
        double divisor = U[k][k];
        if (divisor <= 1.0e-6 && divisor >= -1.0e-6) {
            int pivot = getPartialPivot(U, k);
            U.swapRows(k, pivot);
            L.swapRows(k, pivot);
            indexes.swapElements(k, pivot);
            divisor = U[k][k];

            cout << "Zamiana kolumn " << k << " oraz " << pivot << ": " << endl;
            U.print();
        }

        L[k][k] = 1;

        Vector<N> firstRow = U[k];
        for (int row = k+1; row < N; row++) {
            double firstInCol = U[row][k];
            double multiplier = firstInCol / divisor;

            L[row][k] = multiplier;

            for (int col = k; col < N; col++)
                U[row][col] -= firstRow[col] * multiplier;
        }

        cout << "Po iteracji " << k+1 << ": " << endl;
        U.print();
    }

    L[N-1][N-1] = 1;

    return make_pair(L, U);
}

template <int N>
void rearrangeIndexes(Vector<N> &vect, const Vector<N> &indexes) {
    Vector tmp = vect;
    for (int i = 0; i < N; i++) {
        double val = tmp[i];
        int index = indexes[i];
        vect[index] = val;
    }
};

template <int N, int M>
Vector<N> solveL(Matrix<N, M> &L, Vector<N> &b) {
    Vector<N> x;
    for (int n = 0; n < N; n++) {
        double sum = 0.0;
        for (int m = 0; m <= n-1; m++)
            sum += L[n][m] * x[m];
        x[n] = b[n] - sum;
    }
    return x;
}

template <int N, int M>
Vector<N> solveU(Matrix<N, M> &U, Vector<N> &y) {
    Vector<N> x;
    int sumElementsCount = 0;
    for (int n = N-1; n >= 0; n--) {
        double sum = 0.0;
        int m = N-1;
        for (int _ = 0; _ < sumElementsCount; _++) {
            sum += U[n][m] * x[m];
            m--;
        }
        sumElementsCount++;
        x[n] = (y[n] - sum) / U[n][n];
    }
    return x;
}