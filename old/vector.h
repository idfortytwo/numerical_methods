#include <vector>
#include <iostream>

#ifndef LAB11_VECTOR_H
#define LAB11_VECTOR_H

class Vector {
    std::vector<double> data;

public:
    int N;

    Vector(const Vector &vect) : N(vect.N), data(vect.data) {};
    Vector(std::initializer_list<double> data) : N(data.size()), data(std::vector<double>(data)) {};
    Vector(const std::vector<double> &data) : N(data.size()), data(data) {};
    Vector(int size) : N(size), data(std::vector<double>(size, 0)) {};
    Vector(int size, double value) : N(size), data(std::vector<double>(size, value)) {};

    Vector& operator=(const Vector &vect);
    double& operator[](int i);
    double operator[](int i) const;

    void swapElements(int a, int b) { std::swap(data[a], data[b]); };
    void print(int width = 10);
};


#endif //LAB11_VECTOR_H
