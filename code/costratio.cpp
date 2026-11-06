// costratio.cpp — measure PSO vs SAPSO per-problem cost ratio on target hardware
//
// PURPOSE (freeze §10, pre-flight item 4)
// Coupling Rule B defines T as the time SAPSO takes to exhaust its FE
// allocation. T therefore cannot be set until SAPSO's runtime is measured ON
// THE MACHINE THE MAIN MATRIX WILL RUN ON. Pilot ratios were obtained on
// different hardware and do not transfer: they ranged 1.34x-10.64x there, and
// that spread is itself hardware-dependent.
//
// The ratio is a REPORTED QUANTITY, not a tuning target. SAPSO parameters are
// fixed at the freeze values and are not adjusted to hit a preferred ratio.
// What this program produces is the table that goes in the paper, plus the T
// values for each coupling rule.
//
// Output: one row per (problem, algorithm, run). Ratios and T values are
// computed by the companion script from the raw rows, so nothing is aggregated
// before it is recorded.
//
// EXPLORATORY-CALIBRATION. This measurement fixes a methodology parameter; it
// is not part of the confirmatory matrix and is reported separately.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include "bbob_adapter.hpp"
#include "pso.hpp"
#include "sapso.hpp"
#include "harness.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace tfb;

int main(int argc, char** argv) {
    int dim      = 30;
    int swarm    = 120;
    int runs     = 10;      // enough for a stable mean; this is calibration
    int instance = 1;
    long long femult = 10000;
    int sphere_femult = 0;  // 0 = same as femult; set to raise Sphere's budget
    std::string label = "costratio";
    std::vector<int> fids;   // empty => default_fids()

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char* s) { return !std::strcmp(argv[i], s); };
        if      (eq("--dim")      && i+1 < argc) dim      = std::atoi(argv[++i]);
        else if (eq("--swarm")    && i+1 < argc) swarm    = std::atoi(argv[++i]);
        else if (eq("--runs")     && i+1 < argc) runs     = std::atoi(argv[++i]);
        else if (eq("--femult")   && i+1 < argc) femult   = std::atoll(argv[++i]);
        else if (eq("--spherefemult") && i+1 < argc) sphere_femult = std::atoi(argv[++i]);
        else if (eq("--instance") && i+1 < argc) instance = std::atoi(argv[++i]);
        else if (eq("--label")    && i+1 < argc) label    = argv[++i];
        else if (eq("--fids")     && i+1 < argc) {
            // Comma-separated BBOB function ids, e.g. --fids 20,22,23,24
            fids.clear();
            char* tok = std::strtok(argv[++i], ",");
            while (tok) { fids.push_back(std::atoi(tok)); tok = std::strtok(nullptr, ","); }
        }
    }
    if (fids.empty()) fids = default_fids();

    PSOParams   pp; pp.swarm_size = swarm;
    SAPSOParams sp; sp.swarm_size = swarm;   // freeze values: cap=100, cand=3, refit=5

    std::fprintf(stderr, "# label=%s dim=%d swarm=%d runs=%d instance=%d\n",
                 label.c_str(), dim, swarm, runs, instance);
    std::fprintf(stderr, "# femult=%lld sphere_femult=%d\n", femult,
                 sphere_femult ? sphere_femult : (int)femult);
    std::fprintf(stderr, "# SAPSO archive_cap=%d n_cand=%d refit_every=%d\n",
                 sp.archive_cap, sp.n_cand, sp.refit_every);
    std::fprintf(stderr, "# timer_resolution_s=%.3e\n", timer_resolution_seconds());

    std::printf("label,fid,problem,dim,swarm,algo,fe_budget,run_idx,wall_s,"
                "ns_per_eval,best_fitness,fes\n");

    // Single-threaded throughout: Rule B's T is defined by a single run, and
    // introducing thread count here would confound the calibration with the
    // parallelism factor the main matrix exists to study.
    for (int fid : fids) {
        long long FE = (fid == 1 && sphere_femult)
                         ? (long long)sphere_femult * dim
                         : femult * dim;
        BBOBProblem bp(fid, instance, dim, 1);
        auto prob = bp.as_problem();

        for (int r = 0; r < runs; ++r) {
            auto t0 = Clock::now();
            RunResult a = pso_run(prob, pp, FE, 700000ULL + r, false, 1);
            double ta = seconds_since(t0);
            std::printf("%s,%d,%s,%d,%d,PSO,%lld,%d,%.9f,%.3f,%.17g,%lld\n",
                        label.c_str(), fid, prob.name.c_str(), dim, swarm,
                        FE, r, ta, ta*1e9/(double)a.fes, a.best_fitness, a.fes);

            t0 = Clock::now();
            RunResult b = sapso_run(prob, sp, FE, 700000ULL + r, false, 1);
            double tb = seconds_since(t0);
            std::printf("%s,%d,%s,%d,%d,SAPSO,%lld,%d,%.9f,%.3f,%.17g,%lld\n",
                        label.c_str(), fid, prob.name.c_str(), dim, swarm,
                        FE, r, tb, tb*1e9/(double)b.fes, b.best_fitness, b.fes);
            std::fflush(stdout);
        }
        std::fprintf(stderr, "# done fid=%d %s\n", fid, prob.name.c_str());
    }
    return 0;
}
