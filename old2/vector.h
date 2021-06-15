#include <vector>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <iomanip>

#ifndef VECTOR_H
#define VECTOR_H

class Vector: public std::vector<double> {
public:
    int N;

    Vector(std::initializer_list<double> input) :
        N(input.size()), vector<double>::vector(input) {};
    Vector(int size, double value) : N(size), Vector::vector(size, value) {}
    Vector(int size) : N(size) {
        this->reserve(N);
    }

    friend Vector operator-(const Vector &left, const Vector &right);
    friend Vector operator+(Vector &left, const Vector &right);

    Vector operator*(double scalar) {
        Vector result(N);
        for(int i = 0; i < N; ++i)
            result[i] = (*this)[i] * scalar;
        return result;
    }

    void print(int width=8) {
        for (auto value : (*this))
            std::cout << std::setw(width) << value << " ";
        std::cout << std::endl;
    }

    void print(const std::string &text, int width=10) {
        std::cout << text << std::endl;
        print(width);
    }

    void printVertical(int width=8) {
        for (auto value : (*this))
            std::cout << value << std::endl;
    }

    friend std::ostream& operator<<(std::ostream &os, const Vector &vect);

    double normMax() {
        double maximum = std::abs((*this)[0]);
        for (int i = 1; i < N; i++) {
            double current;
            if ((current = std::abs((*this)[i])) > maximum) {
                maximum = current;
            }
        }
        return maximum;
    }
};

std::ostream &operator<<(std::ostream &os, const Vector &vect) {
    os << '[';
    for (double value : vect)
        os << std::setw(11) << std::setprecision(8) << value << ", ";
    os << "\b\b]";

    return os;
}

Vector operator-(const Vector &left, const Vector &right) {
    int n = left.N;
    Vector result(n);
    for(int i = 0; i < n; i++)
        result[i] = left[i] - right[i];
    return result;
}

Vector operator+(Vector &left, const Vector &right) {
    int n = left.N;
    Vector result(n);
    for(int i = 0; i < n; i++)
        result[i] = left[i] + right[i];
    return result;
}

#endif