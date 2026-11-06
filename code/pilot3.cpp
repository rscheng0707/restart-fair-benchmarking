// pilot3.cpp — Phase 1c: real BBOB suite, corrected affinity levels
//
// Changes from pilot2.cpp, each traceable to a pilot-2 finding:
//
// 1. REAL BBOB PROBLEMS via IOHexperimenter. Pilot 2's hand-written functions
//    gave runs of 20-140 ms, at which scale (a) OpenMP overhead dominated,
//    (b) no thermal drift was detectable in 90/90 configurations, and (c) 29/90
//    configurations had a tau CI half-width above 5%, making k = floor(T/tau)
//    unusable. BBOB applies affine/rotation transforms (O(d^2)), raising
//    per-evaluation cost 4-27x and putting runs in the 0.3-2.3 s range at
//    dim 30 / 300k FE.
//
// 2. SMT LEVEL SEPARATED FROM CORE LEVEL. Pilot 2 ran 32 threads under
//    OMP_PLACES=cores on a 16-core part, forcing two threads onto each core.
//    The resulting variance (pinned CV 38% vs unpinned 10%) measured a policy
//    mismatch, not a placement effect, and must not be reported as one.
//    Thread levels here stop at physical core count by default; the SMT level
//    is a separate labelled pass using OMP_PLACES=threads.
//
// 3. PER-EVALUATION COST REPORTED. Wall time divided by FE count, so the
//    per-problem cost spread (the mechanism the whole study rests on) is
//    measured rather than assumed.
//
// AFFINITY PASSES — run all three, then analyse together:
//    OMP_PROC_BIND=close OMP_PLACES=cores   ./pilot3 --label pinned_cores
//    OMP_PROC_BIND=false                    ./pilot3 --label unpinned
//    OMP_PROC_BIND=close OMP_PLACES=threads ./pilot3 --label pinned_smt --smt
//
// EXPLORATORY. Informs the frozen methodology; not confirmatory evidence.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <random>
#include "bbob_adapter.hpp"
#include "pso.hpp"
#include "harness.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace tfb;

int main(int argc, char** argv) {
    int dim          = 30;
    int swarm        = 120;
    int runs         = 15;
    int instance     = 1;
    long long femult = 10000;
    bool smt         = false;
    uint64_t order_seed = 20260916ULL;
    std::string label = "default";

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char* s) { return !std::strcmp(argv[i], s); };
        if      (eq("--dim")      && i+1 < argc) dim      = std::atoi(argv[++i]);
        else if (eq("--swarm")    && i+1 < argc) swarm    = std::atoi(argv[++i]);
        else if (eq("--runs")     && i+1 < argc) runs     = std::atoi(argv[++i]);
        else if (eq("--femult")   && i+1 < argc) femult   = std::atoll(argv[++i]);
        else if (eq("--instance") && i+1 < argc) instance = std::atoi(argv[++i]);
        else if (eq("--label")    && i+1 < argc) label    = argv[++i];
        else if (eq("--orderseed")&& i+1 < argc) order_seed = std::atoll(argv[++i]);
        else if (eq("--smt"))                    smt = true;
    }

    const long long FE = femult * dim;
    PSOParams p; p.swarm_size = swarm;
    auto fids = default_fids();

#ifdef _OPENMP
    const int hw = omp_get_max_threads();
#else
    const int hw = 1;
#endif
    // Default thread levels stop at physical cores. omp_get_max_threads()
    // reports logical processors, so halve it when SMT is present. The SMT
    // level is opt-in via --smt and must be paired with OMP_PLACES=threads.
    const int cores = smt ? hw : std::max(1, hw / 2);
    std::vector<int> thread_levels;
    for (int t : {1, 2, 4, 6, 8, 12, 16, 24, 32})
        if (t <= cores) thread_levels.push_back(t);
    if (thread_levels.empty()) thread_levels.push_back(1);

    std::fprintf(stderr, "# label=%s smt=%d\n", label.c_str(), (int)smt);
    std::fprintf(stderr, "# timer_resolution_s=%.3e\n", timer_resolution_seconds());
    std::fprintf(stderr, "# dim=%d swarm=%d instance=%d fe_budget=%lld runs=%d\n",
                 dim, swarm, instance, FE, runs);
    std::fprintf(stderr, "# omp_max_threads=%d  thread_ceiling=%d\n", hw, cores);
    std::fprintf(stderr, "# order_seed=%llu\n", (unsigned long long)order_seed);
    const char* pb = std::getenv("OMP_PROC_BIND");
    const char* pl = std::getenv("OMP_PLACES");
    std::fprintf(stderr, "# OMP_PROC_BIND=%s OMP_PLACES=%s\n",
                 pb ? pb : "(unset)", pl ? pl : "(unset)");
    if (smt && (!pl || std::string(pl) != "threads"))
        std::fprintf(stderr, "# WARNING: --smt without OMP_PLACES=threads. "
                             "Pinning to cores while oversubscribing measures a "
                             "policy mismatch, not a placement effect.\n");

    // Construct BBOB problems ONCE, outside any timed region. BBOB construction
    // builds rotation matrices; charging that to the algorithm's wall-clock
    // budget would inflate tau and corrupt k = floor(T/tau).
    std::vector<std::unique_ptr<BBOBProblem>> bbob;
    std::vector<Problem> probs;
    for (int fid : fids) {
        bbob.push_back(std::make_unique<BBOBProblem>(fid, instance, dim, cores + 1));
        probs.push_back(bbob.back()->as_problem());
    }

    struct Cfg { const char* impl; bool omp; int threads; };
    std::vector<Cfg> cfgs = {{"serial", false, 1}};
    for (int t : thread_levels) cfgs.push_back({"omp", true, t});

    struct Task { int pi, ci; };
    std::vector<Task> tasks;
    for (size_t i = 0; i < probs.size(); ++i)
        for (size_t j = 0; j < cfgs.size(); ++j)
            tasks.push_back({(int)i, (int)j});
    std::mt19937_64 rng(order_seed);
    std::shuffle(tasks.begin(), tasks.end(), rng);

    std::printf("label,fid,problem,instance,dim,swarm,impl,threads,fe_budget,"
                "run_idx,task_order,wall_s,ns_per_eval,best_fitness,fes\n");

    auto t_start = Clock::now();
    int order = 0;
    for (const auto& tk : tasks) {
        const auto& prob = probs[tk.pi];
        const auto& c    = cfgs[tk.ci];
        for (int r = 0; r < runs; ++r) {
            auto t0 = Clock::now();
            RunResult rr = pso_run(prob, p, FE, 900000ULL + r, c.omp, c.threads);
            double dt = seconds_since(t0);
            std::printf("%s,%d,%s,%d,%d,%d,%s,%d,%lld,%d,%d,%.9f,%.3f,%.17g,%lld\n",
                        label.c_str(), fids[tk.pi], prob.name.c_str(), instance,
                        dim, swarm, c.impl, c.threads, FE, r, order,
                        dt, dt * 1e9 / (double)rr.fes, rr.best_fitness, rr.fes);
        }
        std::fflush(stdout);
        std::fprintf(stderr, "# [%3d/%3zu] %-18s %-6s t=%-2d elapsed=%.1fs\n",
                     order + 1, tasks.size(), prob.name.c_str(),
                     c.impl, c.threads, seconds_since(t_start));
        ++order;
    }
    std::fprintf(stderr, "# total_wall_s=%.1f\n", seconds_since(t_start));
    return 0;
}
