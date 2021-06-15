//#include <vector>
//#include <algorithm>
//#include <ostream>
//#include <iostream>
//#include <iomanip>
//
//#ifndef VECTOR_H
//#define VECTOR_H
//
//using namespace std;
//
//
////class Vector {
////    double* data;
////
////public:
////    const int N;
////
////    Vector(int N, double default_value) : N(N) {
////        data = new double[N];
////        for (int i = 0; i < N; i++) {
////            data[i] = default_value;
////        }
////    }
////
////    Vector(int N) : N(N) {
////        data = new double[N];
////    }
////
////    Vector(std::initializer_list<double> input) : N(input.size()) {
////        data = new double [N];
////        auto begin = input.begin();
////        for (int i = 0; i < N; i++) {
////            data[i] = *begin;
////            begin++;
////        }
////    }
////
////    Vector(const Vector &vect) : N(vect.N) {
////        data = new double[N];
////        for (int i = 0; i < N; i++) {
////            data[i] = vect.data[i];
////        }
////    }
////
////    ~Vector() {
////        delete[] data;
////    }
////
////    double &operator[](int i) const {
////        return data[i];
////    }
////
////    double &operator[](int i) {
////        return data[i];
////    }
////
////    void operator=(const Vector &vect) {
////        delete[] data;
////
////        data = new double[vect.N];
////        for (int i = 0; i < vect.N; i++) {
////            data[i] = vect.data[i];
////        }
////    }
////
////    void operator=(Vector &vect) {
////        delete[] data;
////
////        data = new double[vect.N];
////        for (int i = 0; i < vect.N; i++) {
////            data[i] = vect.data[i];
////        }
////    }
////
////    void print(int width = 8) {
////        cout << left;
////        for (int i = 0; i < N; i++) {
////            cout << i << ": " << setw(width) << data[i] << endl;
////        }
////        cout << endl << endl;
////    }
////};
////
////#endif
//
//class Vector {
//    vector<double> data;
//
//public:
//    const int N;
//
//    Vector(int N, double default_value) : N(N) {
//        data = new double[N];
//        for (int i = 0; i < N; i++) {
//            data[i] = default_value;
//        }
//    }
//
//    Vector(int N) : N(N) {
//        data = new double[N];
//    }
//
//    Vector(std::initializer_list<double> input) : N(input.size()) {
//        data = new double [N];
//        auto begin = input.begin();
//        for (int i = 0; i < N; i++) {
//            data[i] = *begin;
//            begin++;
//        }
//    }
//
//    Vector(const Vector &vect) : N(vect.N) {
//        data = new double[N];
//        for (int i = 0; i < N; i++) {
//            data[i] = vect.data[i];
//        }
//    }
//
//    ~Vector() {
//        delete[] data;
//    }
//
//    double &operator[](int i) const {
//        return data[i];
//    }
//
//    double &operator[](int i) {
//        return data[i];
//    }
//
//    void operator=(const Vector &vect) {
//        delete[] data;
//
//        data = new double[vect.N];
//        for (int i = 0; i < vect.N; i++) {
//            data[i] = vect.data[i];
//        }
//    }
//
//    void operator=(Vector &vect) {
//        delete[] data;
//
//        data = new double[vect.N];
//        for (int i = 0; i < vect.N; i++) {
//            data[i] = vect.data[i];
//        }
//    }
//
//    void print(int width = 8) {
//        cout << left;
//        for (int i = 0; i < N; i++) {
//            cout << i << ": " << setw(width) << data[i] << endl;
//        }
//        cout << endl << endl;
//    }
//};
//
//#endif