//#include <vector.h>
//#include <algorithm>
//#include <ostream>
//#include <iostream>
//#include <iomanip>
//
//#ifndef MATRIX_H
//#define MATRIX_H
//
//using namespace std;
//
//
////class Matrix {
////    int N, M;
////    double** _arrayofarrays;
////
////public:
////    Matrix(int N, int M) : N(N), M(M) {
////        _arrayofarrays = new double*[N];
////        for(int i = 0; i < N; ++i)
////            _arrayofarrays[i] = new double[M];
////    }
////
////    class Proxy {
////    public:
////        Proxy(double* _array) : _array(_array) { }
////
////        double operator[](int index) {
////            return _array[index];
////        }
////    private:
////        double* _array;
////    };
////
////    Proxy operator[](int index) {
////        return Proxy(_arrayofarrays[index]);
////    }
////};
//
////class Matrix {
////    double** data;
////
////public:
////    const int N, M;
////
////    Matrix(int N, int M) : N(N), M(M) {
////        data = (double**) new double[N * sizeof(int*)];
////        for(int i = 0; i < N; i++)
////            data[i] = new double[M * sizeof(int)];
////    }
////
////    Matrix(int N, int M, double value) : N(N), M(M) {
////        data = (double**) new double[N * sizeof(int*)];
////        for (int i = 0; i < N; i++) {
////            data[i] = new double[M * sizeof(int)];
////            for (int j = 0; j < M; j++) {
////                data[i][j] = value;
////            }
////        }
////    }
////
////    ~Matrix() {
////        for (int i = 0; i < N; i++)
////            delete[] data[i];
////        delete[] data;
////    }
////
////    class Proxy {
////    public:
////        Proxy(double* temp) : temp(temp) {}
////        double& operator[](int i) {
////            return temp[i];
////        }
////        double* temp;
////    };
////
////    Proxy operator[] (int i) {
////        return Proxy(data[i]);
////    }
////};
//
//class Matrix {
//    Vector* data;
//
//public:
//    const int N, M;
//
//    Matrix(int N, int M) : N(N), M(M) {
////        data = (Vector*) new double[N * sizeof(int*)];
//        data = (Vector*) new Vector(N);
//        for(int i = 0; i < N; i++)
//            data[i] = Vector(N);
//    }
//
//    Matrix(int N, int M, double value) : N(N), M(M) {
//        data = (Vector*) new Vector(N);;
//        for (int i = 0; i < N; i++) {
//            data[i] = Vector(N);
//            for (int j = 0; j < M; j++) {
//                data[i][j] = value;
//            }
//        }
//    }
//
//    ~Matrix() {
////        for (int i = 0; i < N; i++)
////            delete data[i];
//        delete[] data;
//    }
//
////    class Proxy {
////    public:
////        Proxy(double* temp) : temp(temp) {}
////        double& operator[](int i) {
////            return temp[i];
////        }
////        double* temp;
////    };
//
//    Vector &operator[] (int i) {
//        return data[i];
//    }
//};
//
//
////class Matrix {
////    int N, M;
////    Vector* m_points;
////
////public:
////    Matrix(int N, int M, double default_value) : N(N), M(M) {
////        m_points = new Vector(M);
////
////        for (int i = 0; i < N; i++) {
////            m_points[i] = default_value;
////        }
////    }
////
////    Matrix(int N, int M) : N(N), M(M) {
////        m_points = new Vector(M);
////    }
////
////    Matrix(std::initializer_list<double> input) : size(input.size()) {
////        m_points = new double [size];
////        auto begin = input.begin();
////        for (int i = 0; i < size; i++) {
////            m_points[i] = *begin;
////            begin++;
////        }
////    }
////
////    Matrix(const Vector &vect) : size(vect.size) {
////        m_points = new double[size];
////        for (int i = 0; i < size; i++) {
////            m_points[i] = vect.m_points[i];
////        }
////    }
////
////    ~Matrix() {
////        delete[] m_points;
////        m_points = nullptr;
////    }
////
////    double &operator[](int i) const {
////        return m_points[i];
////    }
////
////    double &operator[](int i) {
////        return m_points[i];
////    }
////
////    void operator=(Vector &vect) {
////        delete[] m_points;
////
////        size = vect.size;
////        m_points = new double[vect.size];
////        for (int i = 0; i < size; i++) {
////            m_points[i] = vect.m_points[i];
////        }
////    }
////
////    void print(int width = 8) {
////        cout << left;
////        for (int i = 0; i < size; i++) {
////            cout << setw(width) << m_points[i] << endl;
////        }
////        cout << endl << endl;
////    }
////
////    int n() const {
////        return int(size);
////    }
////};
//
//#endif