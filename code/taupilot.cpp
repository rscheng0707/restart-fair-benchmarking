// taupilot.cpp — estimate tau for every main-matrix configuration
//
// Definition 2 requires tau (average single-run duration) to compute
// k = floor(T/tau). Freeze section 4.2 fixes three things about it:
//   - estimated from SEPARATE pilot runs, never from the runs under analysis
//   - estimated PER CONFIGURATION (problem x algorithm/ratio x threads x impl)
//   - fixed at freeze and held constant through the main matrix
//
// Reusing one tau across thread counts or problems would misallocate restarts
// and manufacture a platform effect — the exact artifact RQ2 exists to detect.
//
// Output is tau.csv in the format mainmatrix.cpp expects:
//     fid,ratio,threads,omp,tau_s
// with ratio=0 denoting plain PSO and ratio in {3,6,12} denoting SAPSO at the
// calibrated archive_cap for that level.
//
// A companion diagnostics stream (stderr) reports CV and the 95% CI half-width
// per configuration. Freeze section 4.2 requires the half-width below 5% of tau;
// configurations above that must get more runs or a recorded exclusion.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "bbob_adapter.hpp"
#include "pso.hpp"
#include "sapso.hpp"
#include "harness.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace tfb;

struct ProblemSpec { int fid; long long femult; std::map<int,int> cap; };

// Must match mainmatrix.cpp exactly (freeze 3, 4.5b).
static const std::vector<ProblemSpec> PROBLEMS = {
    {  1, 30000, {{3,   17}, {6,   83}, {12,  199}} },   // femult raised, freeze 3.1; caps recalibrated at femult=30000 (D8)
    {  8, 10000, {{3,  430}, {6,  812}, {12, 1285}} },
    { 10, 10000, {{3,  580}, {6, 1027}, {12, 1604}} },
    { 15, 10000, {{3,  705}, {6, 1242}, {12, 1908}} },
    { 24, 10000, {{3,  473}, {6,  893}, {12, 1401}} },
};
static const std::vector<int> RATIOS  = {0, 3, 6, 12};   // 0 = plain PSO
static const std::vector<int> THREADS = {1, 2, 4, 6, 8, 12, 16};
static const int DIM = 30, SWARM = 120, INSTANCE = 1;
static const int SAPSO_NCAND = 3, SAPSO_REFIT = 5;

// Student t 97.5th percentile, small samples.
static double tcrit(int df) {
    static const double t[] = {0,12.706,4.303,3.182,2.776,2.571,2.447,2.365,
                               2.306,2.262,2.228,2.201,2.179,2.160,2.145,2.131};
    return (df <= 15) ? t[df] : 1.96;
}

int main(int argc, char** argv) {
    int runs = 7;
    std::string out = "tau.csv";
    uint64_t seed_base = 300000ULL;
    std::set<std::string> only;   // empty => run everything

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char* s){ return !std::strcmp(argv[i], s); };
        if      (eq("--runs") && i+1 < argc) runs = std::atoi(argv[++i]);
        else if (eq("--out")  && i+1 < argc) out  = argv[++i];
        else if (eq("--only") && i+1 < argc) {
            // format: fid:ratio:threads:omp  e.g. --only 1:6:16:1,10:0:8:1
            char* t = std::strtok(argv[++i], ",");
            while (t) { only.insert(t); t = std::strtok(nullptr, ","); }
        }
    }

#ifdef _OPENMP
    const int hw = omp_get_max_threads();
#else
    const int hw = 1;
#endif
    const int cores = std::max(1, hw / 2);

    std::fprintf(stderr, "# taupilot: runs=%d dim=%d swarm=%d cores=%d\n",
                 runs, DIM, SWARM, cores);
    std::fprintf(stderr, "# timer_resolution_s=%.3e\n", timer_resolution_seconds());

    std::FILE* fp = std::fopen(out.c_str(), "w");
    if (!fp) { std::fprintf(stderr, "# FATAL: cannot write %s\n", out.c_str()); return 2; }
    std::fprintf(fp, "fid,ratio,threads,omp,tau_s\n");

    int unstable = 0, total = 0;
    auto t_all = Clock::now();

    for (const auto& ps : PROBLEMS) {
        BBOBProblem bp(ps.fid, INSTANCE, DIM, cores + 1);
        auto prob = bp.as_problem();
        const long long FE = ps.femult * DIM;

        for (int ratio : RATIOS) {
            for (int threads : THREADS) {
                if (threads > cores) continue;
                // omp at every thread level; serial only at 1 thread
                std::vector<std::pair<bool,int>> variants = {{true, threads}};
                if (threads == 1) variants.push_back({false, 1});

                for (auto& v : variants) {
                    if (!only.empty()) {
                        char key[64];
                        std::snprintf(key, sizeof key, "%d:%d:%d:%d",
                                      ps.fid, ratio, v.second, v.first ? 1 : 0);
                        if (!only.count(key)) continue;
                    }
                    std::vector<double> t;
                    for (int r = 0; r < runs; ++r) {
                        auto t0 = Clock::now();
                        if (ratio == 0) {
                            PSOParams pp; pp.swarm_size = SWARM;
                            pso_run(prob, pp, FE, seed_base + r, v.first, v.second);
                        } else {
                            SAPSOParams sp; sp.swarm_size = SWARM;
                            sp.archive_cap = ps.cap.at(ratio);
                            sp.n_cand = SAPSO_NCAND; sp.refit_every = SAPSO_REFIT;
                            sapso_run(prob, sp, FE, seed_base + r, v.first, v.second);
                        }
                        t.push_back(seconds_since(t0));
                    }
                    double m = std::accumulate(t.begin(), t.end(), 0.0) / t.size();
                    double sd = 0.0;
                    for (double x : t) sd += (x - m) * (x - m);
                    sd = std::sqrt(sd / (t.size() - 1));
                    double hw_ci = tcrit((int)t.size() - 1) * sd / std::sqrt((double)t.size());
                    double pct = hw_ci / m * 100.0;

                    std::fprintf(fp, "%d,%d,%d,%d,%.9f\n",
                                 ps.fid, ratio, v.second, v.first ? 1 : 0, m);
                    std::fflush(fp);
                    ++total;
                    if (pct > 5.0) {
                        ++unstable;
                        std::fprintf(stderr,
                            "# UNSTABLE fid=%2d ratio=%2d t=%2d %-6s tau=%8.4fs "
                            "CI=%5.2f%% <- needs more runs or a recorded exclusion\n",
                            ps.fid, ratio, v.second, v.first ? "omp" : "serial", m, pct);
                    }
                }
            }
        }
        std::fprintf(stderr, "# done fid=%d %-20s elapsed=%.1f min\n",
                     ps.fid, prob.name.c_str(), seconds_since(t_all) / 60.0);
    }

    std::fclose(fp);
    std::fprintf(stderr, "\n# %d configurations, %d above the 5%% CI threshold\n",
                 total, unstable);
    if (unstable)
        std::fprintf(stderr,
            "# Freeze 4.2: no configuration proceeds with an imprecise tau.\n"
            "# Increase --runs, or record exclusions in the deviations log.\n");
    else
        std::fprintf(stderr, "# All tau estimates within threshold. tau.csv is ready.\n");
    std::fprintf(stderr, "# total %.1f min -> %s\n", seconds_since(t_all)/60.0, out.c_str());
    return 0;
}
