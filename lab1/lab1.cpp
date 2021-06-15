#include <iostream>

using namespace std;

int main() {
    float epsilon_float = 1;
    int t = -1;
    while (epsilon_float + 1 > 1) {
        epsilon_float /= 2.0f;
        t++;
    }

    cout << "mantysa float: " << t << endl;
    cout << "epsilon float: " << epsilon_float << endl;


    double epsilon_double = 1;
    t = -1;
    while (epsilon_double + 1 > 1) {
        epsilon_double /= 2.0;
        t++;
    }

    cout << "mantysa double: " << t << endl;
    cout << "epsilon double: " << epsilon_double << endl;
}