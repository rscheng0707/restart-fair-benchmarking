// calibrate.cpp — find archive_cap values achieving target SAPSO/PSO cost ratios
//
// PURPOSE (freeze §4.4, D2)
// The cost ratio is a controlled experimental factor at levels 2x, 5x, 10x.
// Because surrogate overhead is near-constant while PSO base cost spans 86x
// across the portfolio, a single archive_cap cannot hold the ratio constant.
// This program measures the ratio as a function of archive_cap, per problem,
// so the parameter can be set to hold the TREATMENT constant.
//
// This is controlled design, not tuning toward an outcome: the quantity being
// fixed is the independent variable, and calibration happens before any
// main-matrix run. The resulting table is reported in full, including the
// ACHIEVED ratio alongside the target — discrete archive_cap values will not
// land on 2.00x exactly, and the paper reports what was obtained, not what was
// aimed for.
//
// METHOD
// For each problem: measure PSO tau once (n runs), then sweep archive_cap over
// a geometric ladder, measuring SAPSO tau at each. Cost is dominated by the
// O(m^3) Cholesky refit, so ratio rises steeply with archive_cap; a geometric
// ladder covers the useful range in few points. Selection of the best cap per
// target is done by the companion script from the raw sweep, so nothing is
// aggregated before it is recorded.
//
// COST CONTROL
// Large archive_cap values are expensive by construction. --maxsec caps the
// sweep: once a configuration exceeds it, larger caps for that problem are
// skipped and reported as out of range. Better to record that a target is
// unreachable than to spend an hour proving it.
//
// CALIBRATION PHASE. Exploratory. Fixes a methodology parameter; not part of
// the confirmatory matrix.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include "bbob_adapter.hpp"
#include "pso.hpp"
#include "sapso.hpp"
#include "harness.hpp"

using namespace tfb;

static double mean_of(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
}
static double sd_of(const std::vector<double>& v, double m) {
    if (v.size() < 2) return 0.0;
    double s = 0.0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (double)(v.size() - 1));
}

int main(int argc, char** argv) {
    int dim      = 30;
    int swarm    = 120;
    int runs     = 3;        // per sweep point; the companion script reports CV
    int instance = 1;
    long long femult = 10000;
    double maxsec = 30.0;    // skip sweep points slower than this
    std::string label = "calib";
    std::vector<int> fids;
    // Ladder spans a wide range because the achievable ratio differs by problem:
    // cheap problems reach 10x at small caps, expensive ones need large caps.
    // O(m^3) fit cost means doubling the cap multiplies surrogate work by ~8,
    // so --maxsec bounds the sweep rather than the ladder being conservative.
    std::vector<int> caps = {25, 50, 100, 200, 400, 600, 900, 1400, 2000};

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char* s) { return !std::strcmp(argv[i], s); };
        if      (eq("--dim")      && i+1 < argc) dim      = std::atoi(argv[++i]);
        else if (eq("--swarm")    && i+1 < argc) swarm    = std::atoi(argv[++i]);
        else if (eq("--runs")     && i+1 < argc) runs     = std::atoi(argv[++i]);
        else if (eq("--femult")   && i+1 < argc) femult   = std::atoll(argv[++i]);
        else if (eq("--instance") && i+1 < argc) instance = std::atoi(argv[++i]);
        else if (eq("--maxsec")   && i+1 < argc) maxsec   = std::atof(argv[++i]);
        else if (eq("--label")    && i+1 < argc) label    = argv[++i];
        else if (eq("--fids")     && i+1 < argc) {
            fids.clear();
            char* t = std::strtok(argv[++i], ",");
            while (t) { fids.push_back(std::atoi(t)); t = std::strtok(nullptr, ","); }
        }
        else if (eq("--caps")     && i+1 < argc) {
            caps.clear();
            char* t = std::strtok(argv[++i], ",");
            while (t) { caps.push_back(std::atoi(t)); t = std::strtok(nullptr, ","); }
        }
    }
    if (fids.empty()) fids = {1, 8, 10, 15, 24};   // portfolio after D1

    const long long FE = femult * dim;
    PSOParams pp; pp.swarm_size = swarm;

    std::fprintf(stderr, "# label=%s dim=%d swarm=%d runs=%d femult=%lld\n",
                 label.c_str(), dim, swarm, runs, femult);
    std::fprintf(stderr, "# maxsec=%.1f  caps=", maxsec);
    for (size_t i = 0; i < caps.size(); ++i)
        std::fprintf(stderr, "%d%s", caps[i], i + 1 < caps.size() ? "," : "\n");
    std::fprintf(stderr, "# timer_resolution_s=%.3e\n", timer_resolution_seconds());

    std::printf("label,fid,problem,dim,swarm,archive_cap,n_cand,refit_every,"
                "pso_tau_s,sapso_tau_s,ratio,pso_cv,sapso_cv,n_runs\n");

    for (int fid : fids) {
        BBOBProblem bp(fid, instance, dim, 1);
        auto prob = bp.as_problem();

        // Baseline PSO tau, measured once per problem.
        std::vector<double> tp;
        for (int r = 0; r < runs; ++r) {
            auto t0 = Clock::now();
            pso_run(prob, pp, FE, 800000ULL + r, false, 1);
            tp.push_back(seconds_since(t0));
        }
        double pso_m = mean_of(tp), pso_s = sd_of(tp, pso_m);
        std::fprintf(stderr, "# %-20s PSO tau=%.4fs cv=%.2f%%\n",
                     prob.name.c_str(), pso_m, pso_s / pso_m * 100.0);

        for (int cap : caps) {
            SAPSOParams sp;
            sp.swarm_size  = swarm;
            sp.archive_cap = cap;
            // n_cand and refit_every stay at frozen values; archive_cap is the
            // single calibration knob, so the cost mechanism stays interpretable.

            std::vector<double> ts;
            bool aborted = false;
            for (int r = 0; r < runs; ++r) {
                auto t0 = Clock::now();
                sapso_run(prob, sp, FE, 800000ULL + r, false, 1);
                double dt = seconds_since(t0);
                ts.push_back(dt);
                if (dt > maxsec) { aborted = true; break; }
            }
            double sa_m = mean_of(ts), sa_s = sd_of(ts, sa_m);
            std::printf("%s,%d,%s,%d,%d,%d,%d,%d,%.9f,%.9f,%.6f,%.4f,%.4f,%d\n",
                        label.c_str(), fid, prob.name.c_str(), dim, swarm,
                        cap, sp.n_cand, sp.refit_every,
                        pso_m, sa_m, sa_m / pso_m,
                        pso_s / pso_m * 100.0,
                        sa_m > 0 ? sa_s / sa_m * 100.0 : 0.0,
                        (int)ts.size());
            std::fflush(stdout);
            std::fprintf(stderr, "#   cap=%-5d sapso=%8.3fs ratio=%6.2fx%s\n",
                         cap, sa_m, sa_m / pso_m, aborted ? "  [exceeded maxsec]" : "");
            if (aborted) {
                std::fprintf(stderr, "#   -> skipping larger caps for %s\n",
                             prob.name.c_str());
                break;
            }
        }
    }
    return 0;
}
