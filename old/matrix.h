#include <vector>
#include "vector.h"

#ifndef LAB11_MATRIX_H
#define LAB11_MATRIX_H

class Matrix {
    std::vector<Vector> data;

public:
    int N, M;
    Matrix(const Matrix &matrix) : N(matrix.N), M(matrix.M), data(matrix.data) {};
    Matrix(int N, int M) : N(N), M(M), data(std::vector<Vector>(N, Vector(M))) {};
    Matrix(int N) : N(N), M(N), data(std::vector<Vector>(N, Vector(N))) {};
    Matrix(std::initializer_list<Vector> matrix) : N(matrix.size()), M((*matrix.begin()).getN()),
                                                   data(std::vector<Vector>(matrix)) {};

    Matrix(int N, int M, double value) : N(N), M(M), data(std::vector<Vector>(N, Vector(M, value))) {};

    Vector& operator[](int i) { return data[i]; };
    Vector operator[](int i) const { return data[i]; };
    std::vector<Vector>::iterator begin() { return data.begin(); }
    std::vector<Vector>::iterator end() { return data.end(); }

    void swapRows(int a, int b) { std::swap(data[a], data[b]); };

    void print(int width = 8) const;
};

#endif //LAB11_MATRIX_H
