// sapso.hpp — surrogate-assisted PSO (high per-iteration cost slot)
//
// ROLE IN THE STUDY (paper §3.4)
// This fills the "high per-iteration cost" portfolio slot and anchors coupling
// Rule B: T is defined as the time this algorithm takes to exhaust its FE
// allocation. It must therefore be expensive PER ITERATION, not merely slower
// overall — an algorithm that is slow because it runs more iterations does not
// test the mechanism under study.
//
// It also mirrors the anchor paper's own construction: its illustrative example
// contrasts a plain PSO against an "advanced" PSO with heavier per-iteration
// cost. Using the same base algorithm (PSO) for both the cheap and expensive
// slots isolates per-iteration cost as the treatment, exactly as the RQ3 design
// isolates parallelism by holding the algorithm fixed.
//
// HONESTY REQUIREMENT
// This is a straightforward RBF-surrogate-assisted PSO built from published
// components (Kriging/RBF-assisted EC, cf. Jin 2011 survey). It is NOT presented
// as a novel algorithm and must not be described as one. Its purpose is to be a
// realistic, defensible instance of an expensive-per-iteration method. Report it
// as such.
//
// WHERE THE COST COMES FROM
//   1. RBF model fit over an archive of evaluated points: O(m^3) Cholesky on an
//      m x m matrix, m = archive size (capped).
//   2. Surrogate screening: n_cand candidate positions per particle are
//      predicted by the model, and only the most promising is truly evaluated.
// Both are real surrogate-assisted practice, not artificial padding. Cost is
// tunable via archive_cap and n_cand so that the ratio to plain PSO can be
// SET and REPORTED rather than left to chance.
//
// FE ACCOUNTING
// Surrogate predictions are NOT function evaluations and are not counted. Only
// true objective calls increment the FE counter. This is the standard
// convention and is what makes the FE-fair vs time-fair contrast meaningful:
// under FE-matched budgets the surrogate work is invisible, which is precisely
// the accounting gap the study examines.

#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include "benchmarks.hpp"
#include "pso.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace tfb {

struct SAPSOParams {
    int    swarm_size  = 120;
    double w           = 0.7298;
    double c1          = 1.49618;
    double c2          = 1.49618;
    double vmax_frac   = 0.2;
    int    archive_cap = 100;  // m in the O(m^3) fit; the main cost knob
    int    n_cand      = 3;    // candidates screened per particle per iteration
    int    refit_every = 5;    // refit surrogate every k iterations
    double rbf_eps     = 1.0;   // multiquadric shape parameter
};

// ---- RBF surrogate: multiquadric phi(r) = sqrt(r^2 + eps^2) ----
// Fit solves (Phi + lambda I) c = y by Cholesky. Ridge term keeps the system
// well conditioned when archive points cluster, which they do late in a run.
class RBFSurrogate {
public:
    bool fit(const std::vector<double>& X, const std::vector<double>& y,
             int n, double eps, double ridge = 1e-8)
    {
        m_ = (int)y.size();
        n_ = n;
        eps_ = eps;
        if (m_ < 3) return false;
        X_ = X; y_ = y;

        std::vector<double> A((size_t)m_ * m_);
        for (int i = 0; i < m_; ++i)
            for (int j = 0; j <= i; ++j) {
                double r2 = 0.0;
                for (int d = 0; d < n_; ++d) {
                    double t = X_[(size_t)i * n_ + d] - X_[(size_t)j * n_ + d];
                    r2 += t * t;
                }
                double v = std::sqrt(r2 + eps_ * eps_);
                A[(size_t)i * m_ + j] = A[(size_t)j * m_ + i] = v;
            }
        for (int i = 0; i < m_; ++i) A[(size_t)i * m_ + i] += ridge;

        // Cholesky, in place -> lower triangle. O(m^3): the dominant cost.
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j <= i; ++j) {
                double s = A[(size_t)i * m_ + j];
                for (int k = 0; k < j; ++k)
                    s -= A[(size_t)i * m_ + k] * A[(size_t)j * m_ + k];
                if (i == j) {
                    if (s <= 0.0) return false;   // not PD; skip this refit
                    A[(size_t)i * m_ + i] = std::sqrt(s);
                } else {
                    A[(size_t)i * m_ + j] = s / A[(size_t)j * m_ + j];
                }
            }
        }
        // Solve L L^T c = y
        c_.assign(m_, 0.0);
        for (int i = 0; i < m_; ++i) {
            double s = y_[i];
            for (int k = 0; k < i; ++k) s -= A[(size_t)i * m_ + k] * c_[k];
            c_[i] = s / A[(size_t)i * m_ + i];
        }
        for (int i = m_ - 1; i >= 0; --i) {
            double s = c_[i];
            for (int k = i + 1; k < m_; ++k) s -= A[(size_t)k * m_ + i] * c_[k];
            c_[i] = s / A[(size_t)i * m_ + i];
        }
        ok_ = true;
        return true;
    }

    // Prediction is NOT a function evaluation and is not counted as one.
    double predict(const double* x) const {
        if (!ok_) return std::numeric_limits<double>::infinity();
        double s = 0.0;
        for (int i = 0; i < m_; ++i) {
            double r2 = 0.0;
            for (int d = 0; d < n_; ++d) {
                double t = x[d] - X_[(size_t)i * n_ + d];
                r2 += t * t;
            }
            s += c_[i] * std::sqrt(r2 + eps_ * eps_);
        }
        return s;
    }

    bool ok() const { return ok_; }

private:
    int m_ = 0, n_ = 0;
    double eps_ = 1.0;
    bool ok_ = false;
    std::vector<double> X_, y_, c_;
};

// Surrogate-assisted PSO. Signature mirrors pso_run so the harness can drive
// either without special-casing.
inline RunResult sapso_run(const Problem& prob,
                           const SAPSOParams& p,
                           long long fe_budget,
                           uint64_t seed,
                           bool use_omp,
                           int nthreads)
{
    const int n = prob.dim;
    const int S = p.swarm_size;
    const double range = prob.upper - prob.lower;
    const double vmax = p.vmax_frac * range;

    std::vector<double> pos((size_t)S * n), vel((size_t)S * n), pbest((size_t)S * n);
    std::vector<double> pbest_f(S, std::numeric_limits<double>::infinity());
    std::vector<double> cur_f(S, std::numeric_limits<double>::infinity());
    std::vector<std::mt19937_64> rng(S);
    for (int i = 0; i < S; ++i) rng[i].seed(mix_seed(seed, i));
    std::uniform_real_distribution<double> U(0.0, 1.0);

    // Archive of truly-evaluated points, used to fit the surrogate.
    std::vector<double> archX;
    std::vector<double> archY;

    for (int i = 0; i < S; ++i)
        for (int d = 0; d < n; ++d) {
            pos[(size_t)i * n + d] = prob.lower + U(rng[i]) * range;
            vel[(size_t)i * n + d] = (U(rng[i]) * 2.0 - 1.0) * vmax;
            pbest[(size_t)i * n + d] = pos[(size_t)i * n + d];
        }

    long long fes = 0;
    double gbest_f = std::numeric_limits<double>::infinity();
    std::vector<double> gbest(n, 0.0);

    auto evaluate_all = [&]() {
#ifdef _OPENMP
        if (use_omp) {
            #pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int i = 0; i < S; ++i) cur_f[i] = prob.f(&pos[(size_t)i * n], n);
        } else
#endif
        { for (int i = 0; i < S; ++i) cur_f[i] = prob.f(&pos[(size_t)i * n], n); }
        fes += S;
    };

    auto archive_push = [&](const double* x, double y) {
        if ((int)archY.size() >= p.archive_cap) {
            // Drop oldest — keeps the fit cost bounded and the model local.
            archX.erase(archX.begin(), archX.begin() + n);
            archY.erase(archY.begin());
        }
        archX.insert(archX.end(), x, x + n);
        archY.push_back(y);
    };

    evaluate_all();
    for (int i = 0; i < S; ++i) {
        pbest_f[i] = cur_f[i];
        archive_push(&pos[(size_t)i * n], cur_f[i]);
        if (cur_f[i] < gbest_f) {
            gbest_f = cur_f[i];
            for (int d = 0; d < n; ++d) gbest[d] = pos[(size_t)i * n + d];
        }
    }

    RBFSurrogate model;
    std::vector<double> cand((size_t)p.n_cand * n);
    int iter = 0;

    while (fes + S <= fe_budget) {
        // ---- COST DRIVER 1: surrogate refit, O(m^3) ----
        if (iter % p.refit_every == 0)
            model.fit(archX, archY, n, p.rbf_eps);

        // ---- COST DRIVER 2: candidate screening ----
        // Each particle generates n_cand velocity/position proposals; the
        // surrogate ranks them and only the best becomes the new position.
        // Serial in fixed index order, matching pso.hpp's discipline: RNG draws
        // come from each particle's own stream, so results do not depend on
        // thread scheduling.
        for (int i = 0; i < S; ++i) {
            double best_pred = std::numeric_limits<double>::infinity();
            int best_k = 0;
            for (int k = 0; k < p.n_cand; ++k) {
                for (int d = 0; d < n; ++d) {
                    double r1 = U(rng[i]), r2 = U(rng[i]);
                    double v = p.w * vel[(size_t)i * n + d]
                             + p.c1 * r1 * (pbest[(size_t)i * n + d] - pos[(size_t)i * n + d])
                             + p.c2 * r2 * (gbest[d] - pos[(size_t)i * n + d]);
                    if (v >  vmax) v =  vmax;
                    if (v < -vmax) v = -vmax;
                    double x = pos[(size_t)i * n + d] + v;
                    if (x < prob.lower) x = prob.lower;
                    if (x > prob.upper) x = prob.upper;
                    cand[(size_t)k * n + d] = x;
                }
                double pred = model.ok() ? model.predict(&cand[(size_t)k * n])
                                         : 0.0;
                if (k == 0 || pred < best_pred) { best_pred = pred; best_k = k; }
            }
            for (int d = 0; d < n; ++d) {
                double newx = cand[(size_t)best_k * n + d];
                vel[(size_t)i * n + d] = newx - pos[(size_t)i * n + d];
                pos[(size_t)i * n + d] = newx;
            }
        }

        evaluate_all();

        for (int i = 0; i < S; ++i) {
            archive_push(&pos[(size_t)i * n], cur_f[i]);
            if (cur_f[i] < pbest_f[i]) {
                pbest_f[i] = cur_f[i];
                for (int d = 0; d < n; ++d) pbest[(size_t)i * n + d] = pos[(size_t)i * n + d];
            }
            if (cur_f[i] < gbest_f) {
                gbest_f = cur_f[i];
                for (int d = 0; d < n; ++d) gbest[d] = pos[(size_t)i * n + d];
            }
        }
        ++iter;
    }

    RunResult r;
    r.best_fitness = gbest_f;
    r.fes = fes;
    r.iterations = iter;
    return r;
}

} // namespace tfb
