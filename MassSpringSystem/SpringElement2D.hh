
#pragma once
using namespace std;

class ParametricFunctionBase {
public:
    virtual int n_unknowns() = 0;
    virtual double eval_f(const vector<double> &_x, const vector<double> &_coeffs) = 0;
    virtual void eval_gradient(const vector<double> &_x, const vector<double> &_coeffs, vector<double> &_g) = 0;
    virtual void eval_hessian(const vector<double> &_x, const vector<double> &_coeffs, vector<vector<double>> &_H) = 0;
};

class SpringElement2D : public ParametricFunctionBase {
public:
    SpringElement2D() : ParametricFunctionBase() {}

    inline virtual int n_unknowns() override { return 4; }

    inline virtual double eval_f(const vector<double> &_x, const vector<double> &_coeffs) override {
        double k = _coeffs[0];
        double dx = _x[0] - _x[2];
        double dy = _x[1] - _x[3];
        return 0.5 * k * (dx * dx + dy * dy);
    }

    inline virtual void eval_gradient(const vector<double> &_x, const vector<double> &_coeffs, vector<double> &_g) override {
        double k = _coeffs[0];
        _g[0] = k * (_x[0] - _x[2]);
        _g[1] = k * (_x[1] - _x[3]);
        _g[2] = k * (_x[2] - _x[0]);
        _g[3] = k * (_x[3] - _x[1]);
    }

    inline virtual void eval_hessian(const vector<double> &_x, const vector<double> &_coeffs, vector<vector<double>> &_H) override {
        double k = _coeffs[0];
        _H.resize(4, vector<double>(4, 0.0));
        _H[0][0] = k; _H[1][1] = k; _H[2][2] = k; _H[3][3] = k;
        _H[0][2] = -k; _H[1][3] = -k;
        _H[2][0] = -k; _H[3][1] = -k;
    }
};

int main() {
    SpringElement2D spring;

    // Define the points positions
    vector<double> x = {1.0, 1.0, 3.0, 3.0};
    vector<double> coeffs = {10.0};

    // Calcul of energy
    double energy = spring.eval_f(x, coeffs);
    cout << "Springs energy : " << energy << endl;

    // Calcul of gradient
    vector<double> gradient(4);
    spring.eval_gradient(x, coeffs, gradient);
    cout << "Gradient : ";
    for (double val : gradient) cout << val << " ";
    cout << endl;

    // Calcul of Hessian matrice
    vector<vector<double>> hessian;
    spring.eval_hessian(x, coeffs, hessian);
    cout << "Hessian Matrice:" << endl;
    for (const auto &row : hessian) {
        for (double val : row) cout << val << " ";
        cout << endl;
    }

    return 0;
}