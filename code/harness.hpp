// harness.hpp — implementation of the restart-fair protocol
//
// Implements Algorithm 1 (time-fair evaluation, per instance) and Algorithm 2
// (baseline with independent restarts) from Lian (2025), arXiv:2509.08986.
//
// DEVIATION NOTE (paper §3.2): Algorithm 1 as written exhausts T and counts
// completed runs post hoc; Definition 2 pre-computes k = floor(T/tau). These can
// disagree when run durations vary. We implement BOTH and report both, since
// the protocol is ambiguous on which is authoritative.
//
// TIMING DISCIPLINE:
//   - steady_clock (monotonic); never system_clock, which can jump.
//   - The timer brackets algorithm execution only. Problem setup, allocation of
//     result structures, logging, and I/O are OUTSIDE the timed region.
//   - Allocation inside pso_run IS inside the timed region, deliberately: it is
//     part of what a user pays to run the algorithm, and excluding it would
//     flatter implementations with expensive setup.
//   - Timer resolution is reported; runs shorter than 1000x resolution are
//     flagged, as wall-clock comparison is meaningless at that scale.

#pragma once
#include <chrono>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "benchmarks.hpp"
#include "pso.hpp"

namespace tfb {

using Clock = std::chrono::steady_clock;

inline double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

struct SingleRun {
    double best_fitness;
    double wall_seconds;
    long long fes;
};

struct TimeFairResult {
    double best_of_restarts;      // best over all completed runs
    int completed_runs;           // runs actually finished within T
    int k_predicted;              // floor(T/tau) from pilot-estimated tau
    double total_wall;
    long long total_fes;
    std::vector<SingleRun> runs;
    bool overran;                 // last run exceeded T (see policy note below)
};

// tau estimation from pilot runs. MUST be called per configuration
// (platform x threads x toolchain) — see paper §3.2. A tau measured under one
// configuration yields an incorrect k under another.
struct TauEstimate {
    double mean, sd, median;
    int n;
};

inline TauEstimate estimate_tau(const Problem& prob, const PSOParams& p,
                                long long fe_budget, int pilot_runs,
                                uint64_t seed_base, bool use_omp, int nthreads)
{
    std::vector<double> t;
    t.reserve(pilot_runs);
    for (int i = 0; i < pilot_runs; ++i) {
        auto t0 = Clock::now();
        volatile double sink = pso_run(prob, p, fe_budget,
                                       seed_base + 100000 + i,
                                       use_omp, nthreads).best_fitness;
        (void)sink;
        t.push_back(seconds_since(t0));
    }
    TauEstimate e;
    e.n = pilot_runs;
    e.mean = std::accumulate(t.begin(), t.end(), 0.0) / t.size();
    double v = 0.0;
    for (double x : t) v += (x - e.mean) * (x - e.mean);
    e.sd = std::sqrt(v / std::max<size_t>(1, t.size() - 1));
    std::sort(t.begin(), t.end());
    e.median = t[t.size() / 2];
    return e;
}

// Algorithm 1 / 2: run repeatedly until the wall-clock budget T is exhausted.
//
// RESTART ADMISSION POLICY: a new run is started only if the elapsed time plus
// the estimated tau still fits within T. This prevents starting a run that
// cannot finish, which would waste budget and bias against slow algorithms.
// The policy is stated because the protocol does not specify one, and the
// choice measurably affects results. Reported explicitly in the paper.
inline TimeFairResult run_time_fair(const Problem& prob, const PSOParams& p,
                                    long long fe_budget_per_run,
                                    double T, double tau,
                                    uint64_t seed_base,
                                    bool use_omp, int nthreads,
                                    bool admit_only_if_fits = true)
{
    TimeFairResult r;
    r.best_of_restarts = std::numeric_limits<double>::infinity();
    r.completed_runs = 0;
    r.total_fes = 0;
    r.overran = false;
    r.k_predicted = (tau > 0.0) ? (int)std::floor(T / tau) : 0;

    auto batch_start = Clock::now();
    int j = 0;
    while (true) {
        double elapsed = seconds_since(batch_start);
        if (elapsed >= T) break;
        if (admit_only_if_fits && (elapsed + tau > T) && j > 0) break;

        auto t0 = Clock::now();
        RunResult rr = pso_run(prob, p, fe_budget_per_run,
                               seed_base + j, use_omp, nthreads);
        double dt = seconds_since(t0);

        r.runs.push_back({rr.best_fitness, dt, rr.fes});
        r.best_of_restarts = std::min(r.best_of_restarts, rr.best_fitness);
        r.total_fes += rr.fes;
        ++r.completed_runs;
        ++j;

        if (seconds_since(batch_start) > T) r.overran = true;
    }
    r.total_wall = seconds_since(batch_start);
    return r;
}

// FE-matched regime: one run at the conventional FE budget, no restarts.
inline SingleRun run_fe_matched(const Problem& prob, const PSOParams& p,
                                long long fe_budget, uint64_t seed,
                                bool use_omp, int nthreads)
{
    auto t0 = Clock::now();
    RunResult rr = pso_run(prob, p, fe_budget, seed, use_omp, nthreads);
    return {rr.best_fitness, seconds_since(t0), rr.fes};
}

// Timer resolution check — reported in the paper's environment section.
inline double timer_resolution_seconds() {
    auto a = Clock::now(), b = a;
    while (b == a) b = Clock::now();
    return std::chrono::duration<double>(b - a).count();
}

} // namespace tfb
