#include <array>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <iomanip>

#ifndef VECTOR_H
#define VECTOR_H
using namespace std;

template <int N>
class Vector: public array<double,N> {
public:
    template <int>
    friend Vector<N> operator-(const Vector<N> &left, const Vector<N> &right);

    template <int>
    friend Vector<N> &operator+(Vector<N> &left, const Vector<N> &right);

    Vector<N> operator*(double scalar) {
        Vector<N> result;
        for(int i = 0; i < N; ++i)
            result[i] = (*this)[i] + scalar;
        return result;
    }

    void swapElements(int a, int b) {
        std::swap((*this)[a], (*this)[b]);
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

    template <int>
    friend ostream& operator<<(ostream &os, const Vector<N> &vect);

    double normMax() {
        double maximum = (*this)[0];
        for (int i = 1; i < N; i++) {
            if ((*this)[i] > maximum)
                maximum = (*this)[i];
        }
        return maximum;
    }
};

template <int N>
ostream &operator<<(ostream &os, const Vector<N> &vect) {
    os << '[';
    for (double value : vect)
        os << setw(11) << setprecision(8) << value << ", ";
    os << "\b\b]";

    return os;
}

template <int N>
Vector<N> operator-(const Vector<N> &left, const Vector<N> &right) {
    Vector<N> result;
    for(int i = 0; i < N; i++)
        result[i] = left[i] - right[i];
    return result;
}

template <int N>
Vector<N> &operator+(Vector<N> &left, const Vector<N> &right) {
    Vector<N> result;
    for(int i = 0; i < N; i++)
        result[i] = left[i] + right[i];
    return result;
}

#endif