#include "vector.h"
#include <iomanip>


Vector& Vector::operator=(const Vector &vect) {
    N = vect.N;
    data = vect.data;
    return *this;
};

double &Vector::operator[](int i) {
    return data[i];
}

double Vector::operator[](int i) const {
    return data[i];
}

void Vector::print(int width) {
    for (auto value : data) {
        std::cout << std::setw(width) << value << " ";
    }
    std::cout << std::endl;
}