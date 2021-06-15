#include "matrix.h"

void Matrix::print(int width) const {
    for (Vector row : data)
        row.print(width);
    std::cout << std::endl;
}

