#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;


double f_exact(double t) {
    return 1 - exp(-10*(t + atan(t)));
}

double solveMethod(double left_bound, double right_bound, double h,
                   double (*calculate_y_k1)(double t_k, double y_k, double h),
                   string filename, string header) {

    double y_k1 = 0.0, y_k = 0.0;
    double maxError = 0, currentError;
    double t_k = left_bound;

    ofstream out;
    bool shouldSave = (!filename.empty()) && !(header.empty());
    if (shouldSave) {
        out.open("../lab10/" + filename);
        out << header << endl;
    }

    while (t_k < right_bound) {
        y_k1 = calculate_y_k1(t_k, y_k, h);

        currentError = abs(y_k1 - f_exact(t_k+h));
        if (maxError < currentError) {
            maxError = currentError;
        }

        t_k += h;
        y_k = y_k1;

        if (shouldSave)
            out << t_k << ',' << y_k << endl;
    }

    out.close();

    return maxError;
}


double f_temp(double t){
    return (10.0 * t * t + 20.0) / (t * t + 1.0);
}

double calculate_y_k1_BME(double t_k, double y_k, double h) {
    return y_k + (1.0 - y_k) * f_temp(t_k) * h;
}

double calculate_y_k1_PME(double t_k, double y_k, double h) {
    double temp_k = f_temp(t_k);
    return (y_k + temp_k * h) / (1.0 + temp_k * h);
}

double calculate_y_k1_PMT(double t_k, double y_k, double h) {
    double temp_k = f_temp(t_k);
    double temp_k1 = f_temp(t_k + h);
    return (y_k + h / 2.0 * temp_k + h / 2.0 * temp_k1 - y_k * h / 2.0 * temp_k)
           / (h / 2.0 * temp_k1 + 1.0);
}


double solveBME(double left_bound, double right_bound, double h,
                string filename="", string header="s") {
    return solveMethod(left_bound, right_bound, h, calculate_y_k1_BME,
                       filename, header);
}

double solvePME(double left_bound, double right_bound, double h,
                string filename="", string header="") {
    return solveMethod(left_bound, right_bound, h, calculate_y_k1_PME,
                       filename, header);
}

double solvePMT(double left_bound, double right_bound, double h,
                string filename="", string header="") {
    return solveMethod(left_bound, right_bound, h, calculate_y_k1_PMT,
                       filename, header);
}

double order(double h_i, double h_j, double y_i, double y_j) {
    return (y_i - y_j) / (h_i - h_j);
}

void BME_PME_PMT_errors(double left_bound, double right_bound) {
    ofstream out;
    out.open("../lab10/BME_PME_PMT_errors.csv");
    out << "h,BME,PME,PMT" << endl;

    double error_BME, error_PME, error_PMT;
    double h_0,  BME_0,  PME_0,  PMT_0,
           h_10, BME_10, PME_10, PMT_10;

    int count = 0;
    double h = 0.001;
    while (h > 1e-7) {
        error_BME = solveBME(left_bound, right_bound, h);
        error_PME = solvePME(left_bound, right_bound, h);
        error_PMT = solvePMT(left_bound, right_bound, h);

        out << log10(h) << ','
            << log10(error_BME) << ','
            << log10(error_PME) << ','
            << log10(error_PMT) << endl;

        if (count == 0) {
            h_0 = log10(h);
            BME_0 = log10(error_BME);
            PME_0 = log10(error_PME);
            PMT_0 = log10(error_PMT);
        } else if (count == 10) {
            h_10 = log10(h);
            BME_10 = log10(error_BME);
            PME_10 = log10(error_PME);
            PMT_10 = log10(error_PMT);
        }

        h /= 2.0;
        count++;
    }

    cout << "rząd BME: " << order(h_0, h_10, BME_0, BME_10) << endl;
    cout << "rząd PME: " << order(h_0, h_10, PME_0, PME_10) << endl;
    cout << "rząd PMT: " << order(h_0, h_10, PMT_0, PMT_10) << endl;

    out.close();
}

int main() {
    double left_bound = 0.0, right_bound = 0.75;
    double h = 0.001, h_unstable = 0.25;

    solveBME(0.0, 3.0, h,
             "BME_stable.csv", "t,BME_stable h=0.01");
    solveBME(0.0, 3.0, h_unstable,
             "BME_unstable.csv", "t,BME_unstable h=0.25");
    solvePME(left_bound, right_bound, h,
             "PME.csv", "t,PME");
    solvePMT(left_bound, right_bound, h,
             "PMT.csv", "t,PMT");

    BME_PME_PMT_errors(left_bound, right_bound);
}