// pso.hpp — canonical PSO, serial and OpenMP paths
//
// DESIGN CONSTRAINT (paper §3.4): the parallel and serial implementations must
// produce equivalent search behaviour under matched seeds at thread count 1,
// so that any divergence at higher thread counts is attributable to parallel
// execution rather than to differing search dynamics.
//
// This is achieved by giving each particle its OWN RNG stream, seeded
// deterministically from (base_seed, particle_index). A single shared RNG would
// make results depend on thread scheduling order and destroy reproducibility.
// This is the single most important implementation decision in the harness.
//
// Reference: Kennedy & Eberhart (1995), with inertia weight (Shi & Eberhart 1998)
// and constriction-style coefficients. Cite the canonical specification in the
// paper; do not present this as a novel variant.

#pragma once
#include <vector>
#include <random>
#include <limits>
#include <cstdint>
#include "benchmarks.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace tfb {

struct PSOParams {
    int swarm_size   = 40;
    double w         = 0.7298;   // inertia
    double c1        = 1.49618;  // cognitive
    double c2        = 1.49618;  // social
    double vmax_frac = 0.2;      // velocity clamp as fraction of range
};

struct RunResult {
    double best_fitness;
    long long fes;
    int iterations;
    std::vector<std::pair<long long, double>> fe_trace;  // (FE count, best-so-far)
};

// Deterministic per-particle seeding. Splitmix64 finalizer gives good
// dispersion so that adjacent particle indices get unrelated streams.
inline uint64_t mix_seed(uint64_t base, uint64_t idx) {
    uint64_t z = base + 0x9E3779B97F4A7C15ULL * (idx + 1);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// Core PSO. `nthreads` == 1 with use_omp=false is the serial reference path.
// Termination: FE budget. (Time-budgeted termination is handled by the harness
// wrapping this in the Algorithm 1 restart loop.)
inline RunResult pso_run(const Problem& prob,
                         const PSOParams& p,
                         long long fe_budget,
                         uint64_t seed,
                         bool use_omp,
                         int nthreads,
                         bool trace = false,
                         int trace_every = 100)
{
    const int n = prob.dim;
    const int S = p.swarm_size;
    const double range = prob.upper - prob.lower;
    const double vmax = p.vmax_frac * range;

    std::vector<double> pos(S * n), vel(S * n), pbest(S * n);
    std::vector<double> pbest_f(S, std::numeric_limits<double>::infinity());
    std::vector<double> cur_f(S, std::numeric_limits<double>::infinity());
    std::vector<std::mt19937_64> rng(S);

    // Per-particle streams — identical regardless of thread count or schedule.
    for (int i = 0; i < S; ++i) rng[i].seed(mix_seed(seed, i));

    std::uniform_real_distribution<double> U(0.0, 1.0);

    // Initialisation (serial, deterministic order)
    for (int i = 0; i < S; ++i) {
        for (int d = 0; d < n; ++d) {
            pos[i * n + d] = prob.lower + U(rng[i]) * range;
            vel[i * n + d] = (U(rng[i]) * 2.0 - 1.0) * vmax;
            pbest[i * n + d] = pos[i * n + d];
        }
    }

    long long fes = 0;
    double gbest_f = std::numeric_limits<double>::infinity();
    std::vector<double> gbest(n, 0.0);
    RunResult res;
    int iter = 0;

    // Initial evaluation
    auto evaluate_all = [&]() {
#ifdef _OPENMP
        if (use_omp) {
            #pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int i = 0; i < S; ++i)
                cur_f[i] = prob.f(&pos[i * n], n);
        } else
#endif
        {
            for (int i = 0; i < S; ++i)
                cur_f[i] = prob.f(&pos[i * n], n);
        }
        fes += S;
    };

    evaluate_all();
    for (int i = 0; i < S; ++i) {
        pbest_f[i] = cur_f[i];
        if (cur_f[i] < gbest_f) {
            gbest_f = cur_f[i];
            for (int d = 0; d < n; ++d) gbest[d] = pos[i * n + d];
        }
    }

    while (fes + S <= fe_budget) {
        // Position/velocity update.
        // RNG draws happen per particle from that particle's own stream, so the
        // sequence of random numbers each particle sees is independent of how
        // the loop is scheduled across threads.
#ifdef _OPENMP
        if (use_omp) {
            #pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int i = 0; i < S; ++i) {
                for (int d = 0; d < n; ++d) {
                    double r1 = U(rng[i]), r2 = U(rng[i]);
                    double v = p.w * vel[i * n + d]
                             + p.c1 * r1 * (pbest[i * n + d] - pos[i * n + d])
                             + p.c2 * r2 * (gbest[d]        - pos[i * n + d]);
                    if (v >  vmax) v =  vmax;
                    if (v < -vmax) v = -vmax;
                    vel[i * n + d] = v;
                    double x = pos[i * n + d] + v;
                    if (x < prob.lower) { x = prob.lower; vel[i * n + d] = 0.0; }
                    if (x > prob.upper) { x = prob.upper; vel[i * n + d] = 0.0; }
                    pos[i * n + d] = x;
                }
            }
        } else
#endif
        {
            for (int i = 0; i < S; ++i) {
                for (int d = 0; d < n; ++d) {
                    double r1 = U(rng[i]), r2 = U(rng[i]);
                    double v = p.w * vel[i * n + d]
                             + p.c1 * r1 * (pbest[i * n + d] - pos[i * n + d])
                             + p.c2 * r2 * (gbest[d]        - pos[i * n + d]);
                    if (v >  vmax) v =  vmax;
                    if (v < -vmax) v = -vmax;
                    vel[i * n + d] = v;
                    double x = pos[i * n + d] + v;
                    if (x < prob.lower) { x = prob.lower; vel[i * n + d] = 0.0; }
                    if (x > prob.upper) { x = prob.upper; vel[i * n + d] = 0.0; }
                    pos[i * n + d] = x;
                }
            }
        }

        evaluate_all();

        // Reduction done serially in fixed index order so that ties resolve
        // identically regardless of thread count. A parallel reduction with
        // nondeterministic ordering would break serial/parallel equivalence.
        for (int i = 0; i < S; ++i) {
            if (cur_f[i] < pbest_f[i]) {
                pbest_f[i] = cur_f[i];
                for (int d = 0; d < n; ++d) pbest[i * n + d] = pos[i * n + d];
            }
            if (cur_f[i] < gbest_f) {
                gbest_f = cur_f[i];
                for (int d = 0; d < n; ++d) gbest[d] = pos[i * n + d];
            }
        }

        ++iter;
        if (trace && (iter % trace_every == 0))
            res.fe_trace.emplace_back(fes, gbest_f);
    }

    res.best_fitness = gbest_f;
    res.fes = fes;
    res.iterations = iter;
    return res;
}

} // namespace tfb
