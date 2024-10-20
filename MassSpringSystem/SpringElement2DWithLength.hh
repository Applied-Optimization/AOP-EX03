#pragma once
using namespace std;

#include <Downloads/aopt-exercise3/include/FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {
class ParametricFunctionBase {
public:
    virtual int n_unknowns() = 0;
    virtual double eval_f(const vector<double> &_x, const vector<double> &_coeffs) = 0;
    virtual void eval_gradient(const vector<double> &_x, const vector<double> &_coeffs, vector<double> &_g) = 0;
    virtual void eval_hessian(const vector<double> &_x, const vector<double> &_coeffs, vector<vector<double>> &_H) = 0;
};

class SpringElement2DWithLength : public ParametricFunctionBase {
public:
    SpringElement2DWithLength() : ParametricFunctionBase() {}

    inline virtual int n_unknowns() override { return 4; }

    inline virtual double eval_f(const vector<double> &_x, const vector<double> &_coeffs) override {
        double k = _coeffs[0];
        double l = _coeffs[1];
        double dx = _x[0] - _x[2];
        double dy = _x[1] - _x[3];
        double dist_squared = dx * dx + dy * dy;
        return 0.5 * k * ((dist_squared - l * l) * (dist_squared - l * l));
    }

    inline virtual void eval_gradient(const vector<double> &_x, const vector<double> &_coeffs, vector<double> &_g) override {
        double k = _coeffs[0];
        double l = _coeffs[1];
        double dx = _x[0] - _x[2];
        double dy = _x[1] - _x[3];
        double dist_squared = dx * dx + dy * dy;
        double diff = dist_squared - l * l;
        _g[0] = 2 * k * diff * dx;
        _g[1] = 2 * k * diff * dy;
        _g[2] = -2 * k * diff * dx;
        _g[3] = -2 * k * diff * dy;
    }

    inline virtual void eval_hessian(const vector<double> &_x, const vector<double> &_coeffs, vector<vector<double>> &_H) override {
        double k = _coeffs[0];
        double l = _coeffs[1];
        double dx = _x[0] - _x[2];
        double dy = _x[1] - _x[3];
        double dist_squared = dx * dx + dy * dy;
        double diff = dist_squared - l * l;

        _H.resize(4, vector<double>(4, 0.0));
        _H[0][0] = 2 * k * (2 * dx * dx + diff);
        _H[1][1] = 2 * k * (2 * dy * dy + diff);
        _H[2][2] = 2 * k * (2 * dx * dx + diff);
        _H[3][3] = 2 * k * (2 * dy * dy + diff);
        _H[0][2] = _H[2][0] = -2 * k * (2 * dx * dx + diff);
        _H[1][3] = _H[3][1] = -2 * k * (2 * dy * dy + diff);
    }
};

int main() {
    SpringElement2DWithLength spring;

    // Define the points positions
    vector<double> x = {1.0, 1.0, 3.0, 3.0};
    vector<double> coeffs = {10.0, 2.0};  // k = 10, l = 2 

    // Calcul of energy
    double energy = spring.eval_f(x, coeffs);
    cout << "Energy with length : " << energy << endl;

    // Calcul of gradient
    vector<double> gradient(4);
    spring.eval_gradient(x, coeffs, gradient);
    cout << "Gradient : ";
    for (double val : gradient) cout << val << " ";
    cout << endl;

    // Calcul of Hessian matrice
    vector<vector<double>> hessian;
    spring.eval_hessian(x, coeffs, hessian);
    cout << "Hessian matrice:" << endl;
    for (const auto &row : hessian) {
        for (double val : row) cout << val << " ";
        cout << endl;
    }

    return 0;
}
} // namespace AOPT



