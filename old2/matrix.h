#include "vector.h"

#ifndef MATRIX_H
#define MATRIX_H
using namespace std;

template<int N, int M>
class Matrix : array<Vector<N>, M> {
public:
    Matrix<N, M>() = default;
    Matrix<N, M>(array<Vector<N>, M> a) : array<Vector<N>, M> (a) {}

    using array<Vector<N>, M>::operator[];

    void swapRows(int a, int b) {
        swap((*this)[a], (*this)[b]);
    };

    void print(int width=8) const {
        for (Vector row : *this)
            row.print(width);
        std::cout << std::endl;
    }

    void print(string text, int width=8) const {
        std::cout << text << std::endl;
        this->print(width);
    }

    Matrix<N, M> operator*(double scalar) {
        Matrix<N, M> result;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                result[i][j] = (*this)[i][j] * scalar;
        return result;
    }

    template <int, int>
    friend Vector<N> operator*(const Matrix<N, M> m, const Vector<M> &x);

    Matrix<N, M> operator+(const Matrix<N, M> &other) {
        Matrix<N, M> result;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                result[i][j] = (*this)[i][j] + other[i][j];
        return result;
    }

    Matrix<N, M> operator*(const Matrix<N, M> &other) {
        Matrix<N, M> result;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                double sum = 0.0;
                for (int k = 0; k < M; k++) {
                    sum += (*this)[i][k] * other[k][j];
                }
                result[i][j] = sum;
            }
        }

        return result;
    }
};

template <int N, int M>
Vector<N> operator*(const Matrix<N, M> m, const Vector<M> &x) {
    Vector<N> result;
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < M; j++)
            sum += m[i][j] * x[j];
        result[i] = sum;
    }
    return result;
}

#endif