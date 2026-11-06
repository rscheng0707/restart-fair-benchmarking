// benchmarks.hpp — test functions for time-fair benchmarking harness
//
// NOTE ON SCOPE: these are hand-implemented standard functions used for
// harness development and pilot runs only. The main matrix should use the
// official BBOB/COCO suite (via IOHexperimenter) so that problem definitions,
// instance transformations, and target values are the community-standard ones.
// Reviewers will expect BBOB, not hand-rolled Rastrigin.

#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <functional>
#include <atomic>

namespace tfb {

// Function-evaluation counter. Every objective call increments it.
// Made explicit rather than global so that FE accounting is auditable —
// FE traces must be reported alongside time traces (see paper §3.7).
struct FECounter {
    long long count = 0;
    void reset() { count = 0; }
};

struct Problem {
    std::string name;
    int dim;
    double lower, upper;
    std::function<double(const double*, int)> f;
    double optimum;   // known global optimum value (0 for all below)
};

// ---- landscape: unimodal, separable ----
inline double sphere(const double* x, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i] * x[i];
    return s;
}

// ---- landscape: unimodal, non-separable, ill-conditioned ----
inline double rosenbrock(const double* x, int n) {
    double s = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        double a = x[i + 1] - x[i] * x[i];
        double b = 1.0 - x[i];
        s += 100.0 * a * a + b * b;
    }
    return s;
}

// ---- landscape: multimodal, separable, regular structure ----
inline double rastrigin(const double* x, int n) {
    double s = 10.0 * n;
    for (int i = 0; i < n; ++i)
        s += x[i] * x[i] - 10.0 * std::cos(2.0 * M_PI * x[i]);
    return s;
}

// ---- landscape: multimodal, non-separable ----
inline double ackley(const double* x, int n) {
    double s1 = 0.0, s2 = 0.0;
    for (int i = 0; i < n; ++i) {
        s1 += x[i] * x[i];
        s2 += std::cos(2.0 * M_PI * x[i]);
    }
    return -20.0 * std::exp(-0.2 * std::sqrt(s1 / n))
           - std::exp(s2 / n) + 20.0 + M_E;
}

// ---- landscape: multimodal, weak global structure ----
inline double griewank(const double* x, int n) {
    double s = 0.0, p = 1.0;
    for (int i = 0; i < n; ++i) {
        s += x[i] * x[i] / 4000.0;
        p *= std::cos(x[i] / std::sqrt((double)(i + 1)));
    }
    return s - p + 1.0;
}

inline std::vector<Problem> standard_suite(int dim) {
    return {
        {"sphere",     dim, -5.0,   5.0,   sphere,     0.0},
        {"rosenbrock", dim, -5.0,   10.0,  rosenbrock, 0.0},
        {"rastrigin",  dim, -5.12,  5.12,  rastrigin,  0.0},
        {"ackley",     dim, -32.77, 32.77, ackley,     0.0},
        {"griewank",   dim, -600.0, 600.0, griewank,   0.0},
    };
}

} // namespace tfb
