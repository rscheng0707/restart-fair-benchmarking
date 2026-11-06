// mainmatrix.cpp — Phase 2 confirmatory experiment
//
// Implements the frozen design (ICMLSC2027_methodology_freeze.md).
// EVERYTHING here follows the freeze. Any change to a factor level, coupling
// rule, tau value, or calibration constant is a DEVIATION and must be logged
// in freeze section 11 with a reason. Do not silently edit constants.
//
// ---------------------------------------------------------------------------
// DESIGN SUMMARY
//   Budget regimes : FE-matched | time-fair (restart-fair) | multi-budget FE
//   Coupling rules : B (primary, T = SAPSO run time) | A (T = PSO) | C (mean)
//   Cost ratio     : 3x | 6x | 12x   (per-problem archive_cap, freeze 4.5b)
//   Algorithms     : PSO | SAPSO
//   Implementation : serial | OpenMP
//   Threads        : 1 2 4 6 8 12 16
//   Problems       : BBOB f1 f8 f10 f15 f24, instance 1, dim 30
//   Seeds          : 15, run as 5 sessions x 3 seeds (freeze 6.0)
//
// SESSION BLOCKING (freeze 6.0)
// Sessions split by SEED, not by factor. Each session runs the COMPLETE
// condition set at its seed subset, so session is a balanced block: a hot
// evening or a background process inflates variance but cannot bias any factor
// comparison. Splitting by coupling rule or problem would confound session with
// that factor irrecoverably.
//
//   mainmatrix.exe --session 1 --seeds-from 1  --seeds-to 3
//   mainmatrix.exe --session 2 --seeds-from 4  --seeds-to 6
//   ... through session 5, seeds 13-15.
//
// CHECKPOINTING
// One CSV row per measured run, flushed after every condition. An interrupted
// session costs the remaining conditions, not the completed ones. --resume
// reads an existing output file and skips conditions already present.
//
// TAU
// Read from a tau file (--taufile), NOT recomputed here. Definition 2 requires
// tau estimated from separate pilot runs and held fixed; estimating it from the
// runs under analysis would fit the restart count to the data it evaluates.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include "bbob_adapter.hpp"
#include "pso.hpp"
#include "sapso.hpp"
#include "harness.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace tfb;

// ---------------------------------------------------------------------------
// FROZEN CONSTANTS — freeze sections 3, 4.5b, 5.  Do not edit without a
// deviation entry.
// ---------------------------------------------------------------------------

struct ProblemSpec {
    int fid;
    long long femult;                 // FE budget multiplier (x dim)
    std::map<int, int> cap_for_ratio; // target ratio -> archive_cap (freeze 4.5b)
};

// archive_cap values verified 16 Aug 2026; achieved ratios in freeze 4.5b.
static const std::vector<ProblemSpec> PROBLEMS = {
    // Sphere: femult raised to 30000 (freeze 3.1). At the standard 10000 its
    // single-run time is 0.292 s, short enough that timer noise inflates the
    // tau CI (8-12% observed during calibration). 30000 gives ~0.88 s.
    // NOTE: archive_cap values were calibrated at femult=10000. Because the
    // cost ratio is a ratio of two times that both scale with FE count, it is
    // approximately invariant to the multiplier - but this MUST be verified
    // (see freeze 3.1 verification requirement) before Session 1.
    {  1, 30000, {{3,   17}, {6,   83}, {12,  199}} },  // Sphere (caps recalibrated at femult=30000, D8)
    {  8, 10000, {{3,  430}, {6,  812}, {12, 1285}} },  // Rosenbrock
    { 10, 10000, {{3,  580}, {6, 1027}, {12, 1604}} },  // EllipsoidRotated
    { 15, 10000, {{3,  705}, {6, 1242}, {12, 1908}} },  // RastriginRotated
    { 24, 10000, {{3,  473}, {6,  893}, {12, 1401}} },  // LunacekBiRastrigin
};

static const std::vector<int> RATIO_LEVELS  = {3, 6, 12};
static const std::vector<int> THREAD_LEVELS = {1, 2, 4, 6, 8, 12, 16};
// Multi-budget FE regime: Piotrowski et al. (2024) argue no single FE budget is
// privileged. Multipliers relative to the conventional 10,000 x D.
//
// DEVIATION D4 (16 Aug 2026): the 10x level is dropped. At full enumeration it
// accounted for 71.6% of total compute (22.3 h -> 8.0 h without it), pushing
// sessions from ~1.6 h to ~4.5 h and defeating the session-blocking design
// (freeze 6.0). The regime is retained at 0.1x and 1.0x, spanning one order of
// magnitude rather than Piotrowski et al.'s three. The narrower span is a
// stated limitation: it tests budget sensitivity but cannot claim the same
// range as the study it responds to.
static const std::vector<double> MULTIBUDGET_FACTORS = {0.1, 1.0};

static const int    DIM        = 30;
static const int    SWARM      = 120;
static const int    INSTANCE   = 1;
static const int    SAPSO_NCAND = 3;
static const int    SAPSO_REFIT = 5;

// ---------------------------------------------------------------------------

struct TauKey {
    int fid, ratio, threads;   // ratio 0 => PSO (no surrogate)
    bool omp;
    bool operator<(const TauKey& o) const {
        return std::tie(fid, ratio, threads, omp) <
               std::tie(o.fid, o.ratio, o.threads, o.omp);
    }
};

// tau file format (CSV, header required):
//   fid,ratio,threads,omp,tau_s
// ratio=0 denotes plain PSO. omp=0/1.
static std::map<TauKey, double> load_tau(const std::string& path) {
    std::map<TauKey, double> m;
    std::ifstream f(path);
    if (!f) {
        std::fprintf(stderr, "# FATAL: cannot open tau file '%s'\n", path.c_str());
        std::fprintf(stderr, "# tau must come from separate pilot runs (freeze 4.2);\n"
                             "# it is not estimated from the runs under analysis.\n");
        std::exit(2);
    }
    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok;
        int v[4]; double tau;
        for (int i = 0; i < 4; ++i) { std::getline(ss, tok, ','); v[i] = std::atoi(tok.c_str()); }
        std::getline(ss, tok, ','); tau = std::atof(tok.c_str());
        m[TauKey{v[0], v[1], v[2], v[3] != 0}] = tau;
    }
    std::fprintf(stderr, "# loaded %zu tau entries from %s\n", m.size(), path.c_str());
    return m;
}

static double tau_lookup(const std::map<TauKey, double>& m,
                         int fid, int ratio, int threads, bool omp)
{
    auto it = m.find(TauKey{fid, ratio, threads, omp});
    if (it == m.end()) {
        std::fprintf(stderr,
            "# FATAL: no tau for fid=%d ratio=%d threads=%d omp=%d.\n"
            "# Every configuration needs a pilot-estimated tau (freeze 4.2).\n",
            fid, ratio, threads, (int)omp);
        std::exit(2);
    }
    return it->second;
}

static std::set<std::string> load_done(const std::string& path) {
    std::set<std::string> done;
    std::ifstream f(path);
    if (!f) return done;
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> c;
        std::stringstream ss(line);
        std::string t;
        while (std::getline(ss, t, ',')) c.push_back(t);
        if (c.size() < 12) continue;
        // columns: session,regime,rule,ratio,fid,problem,algo,impl,threads,mb,seed,...
        done.insert(c[0]+"|"+c[1]+"|"+c[2]+"|"+c[3]+"|"+c[4]+"|"+c[6]+"|"+c[7]+
                    "|"+c[8]+"|"+c[9]+"|"+c[10]);
    }
    return done;
}

int main(int argc, char** argv) {
    int session = 0, seed_from = 1, seed_to = 15;
    uint64_t order_seed = 20260916ULL;
    std::string taufile = "tau.csv", out = "", label = "main";
    bool resume = false, dry = false;
    double warmup_sec = 600.0;   // freeze 6.2; 0 disables (logged as deviation)

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char* s) { return !std::strcmp(argv[i], s); };
        if      (eq("--session")    && i+1 < argc) session   = std::atoi(argv[++i]);
        else if (eq("--seeds-from") && i+1 < argc) seed_from = std::atoi(argv[++i]);
        else if (eq("--seeds-to")   && i+1 < argc) seed_to   = std::atoi(argv[++i]);
        else if (eq("--taufile")    && i+1 < argc) taufile   = argv[++i];
        else if (eq("--out")        && i+1 < argc) out       = argv[++i];
        else if (eq("--label")      && i+1 < argc) label     = argv[++i];
        else if (eq("--orderseed")  && i+1 < argc) order_seed= std::atoll(argv[++i]);
        else if (eq("--warmup")     && i+1 < argc) warmup_sec= std::atof(argv[++i]);
        else if (eq("--resume"))                   resume    = true;
        else if (eq("--dry-run"))                  dry       = true;
    }
    if (session == 0) {
        std::fprintf(stderr, "usage: mainmatrix --session N --seeds-from A --seeds-to B\n"
                             "                  [--taufile tau.csv] [--out file.csv]\n"
                             "                  [--warmup SEC] [--resume] [--dry-run]\n");
        return 1;
    }
    if (out.empty()) out = "session" + std::to_string(session) + ".csv";

    auto TAU = dry ? std::map<TauKey,double>{} : load_tau(taufile);
    std::set<std::string> done;
    if (resume) {
        done = load_done(out);
        std::fprintf(stderr, "# resume: %zu conditions already recorded in %s\n",
                     done.size(), out.c_str());
    }

#ifdef _OPENMP
    const int hw = omp_get_max_threads();
#else
    const int hw = 1;
#endif
    const int cores = std::max(1, hw / 2);   // physical cores; SMT excluded (freeze 5)

    // ---- environment metadata (freeze 6.7 / protocol checklist item 6) ----
    std::fprintf(stderr, "# label=%s session=%d seeds=%d-%d\n",
                 label.c_str(), session, seed_from, seed_to);
    std::fprintf(stderr, "# dim=%d swarm=%d instance=%d\n", DIM, SWARM, INSTANCE);
    std::fprintf(stderr, "# omp_max_threads=%d physical_cores=%d\n", hw, cores);
    std::fprintf(stderr, "# timer_resolution_s=%.3e\n", timer_resolution_seconds());
    const char* pb = std::getenv("OMP_PROC_BIND");
    const char* pl = std::getenv("OMP_PLACES");
    std::fprintf(stderr, "# OMP_PROC_BIND=%s OMP_PLACES=%s\n",
                 pb ? pb : "(unset)", pl ? pl : "(unset)");
    std::fprintf(stderr, "# order_seed=%llu warmup_sec=%.0f\n",
                 (unsigned long long)order_seed, warmup_sec);
    // Seed-space audit: pilots use 300000 (taupilot), 700000 (costratio),
    // 800000 (calibrate); the main matrix uses 500000 + seed*1000 + restart.
    // Ranges do not overlap, so tau is estimated on different streams than the
    // runs it is applied to (freeze 4.2).
    std::fprintf(stderr, "# seed_space=main:500000+seed*1000+restart\n");

    // ---- build problems once, outside any timed region ----
    // BBOB construction generates rotation matrices; charging that to an
    // algorithm's wall-clock budget would inflate tau and corrupt k.
    std::vector<std::unique_ptr<BBOBProblem>> holders;
    std::map<int, Problem> probs;
    for (const auto& ps : PROBLEMS) {
        holders.push_back(std::make_unique<BBOBProblem>(ps.fid, INSTANCE, DIM, cores + 1));
        probs[ps.fid] = holders.back()->as_problem();
    }

    // ---- enumerate the condition list, then permute ----
    // Interleaving prevents any factor level from occupying a contiguous
    // thermal window (freeze 6.2).
    struct Cond {
        const char* regime;   // "fe" | "timefair" | "multife"
        char rule;            // 'B' | 'A' | 'C' | '-'
        int ratio;            // 3 | 6 | 12
        int fid;
        const char* algo;     // "PSO" | "SAPSO"
        const char* impl;     // "serial" | "omp"
        int threads;
        double mb;            // multi-budget factor; 1.0 otherwise
        int seed;
    };
    std::vector<Cond> conds;

    for (int seed = seed_from; seed <= seed_to; ++seed)
    for (const auto& ps : PROBLEMS)
    for (int ratio : RATIO_LEVELS)
    for (int threads : THREAD_LEVELS) {
        if (threads > cores) continue;
        // OpenMP at every thread level; serial reference only at 1 thread.
        // The serial/OpenMP pair at 1 thread is the RQ3 control (freeze 3.4):
        // both must produce bit-identical results, verified by test_equiv.
        std::vector<std::pair<const char*, int>> impls = {{"omp", threads}};
        if (threads == 1) impls.push_back({"serial", 1});

        for (auto& ip : impls) {
            for (const char* algo : {"PSO", "SAPSO"}) {
                const bool is_sapso = std::strcmp(algo, "SAPSO") == 0;

                // PSO does not use archive_cap, so under the FE-matched and
                // multi-budget regimes its run is IDENTICAL across cost-ratio
                // levels: same parameters, same FE budget, same seed. Emitting
                // it once per level would write three duplicate rows (2,400
                // across the experiment) that analysis would have to detect and
                // discard. Under the time-fair regime PSO DOES differ by level,
                // because T is derived from SAPSO's runtime at that level and
                // therefore sets a different restart count k.
                //
                // So: budget-only regimes are emitted once for PSO (tagged with
                // the sentinel ratio 0); time-fair is emitted per level for both
                // algorithms.
                const bool emit_budget_regimes = is_sapso || (ratio == RATIO_LEVELS.front());
                const int  budget_ratio_tag    = is_sapso ? ratio : 0;

                if (emit_budget_regimes) {
                    conds.push_back({"fe", '-', budget_ratio_tag, ps.fid, algo,
                                     ip.first, ip.second, 1.0, seed});
                    for (double mb : MULTIBUDGET_FACTORS) {
                        if (mb == 1.0) continue;   // already covered by "fe"
                        conds.push_back({"multife", '-', budget_ratio_tag, ps.fid,
                                         algo, ip.first, ip.second, mb, seed});
                    }
                }

                for (char rule : {'B', 'A', 'C'})
                    conds.push_back({"timefair", rule, ratio, ps.fid, algo,
                                     ip.first, ip.second, 1.0, seed});
            }
        }
    }

    std::mt19937_64 orng(order_seed + session);
    std::shuffle(conds.begin(), conds.end(), orng);
    std::fprintf(stderr, "# %zu conditions this session\n", conds.size());

    if (dry) {
        std::fprintf(stderr, "# dry run: no execution. Estimated conditions above.\n");
        return 0;
    }

    // ---- output ----
    bool fresh = !resume || done.empty();
    std::FILE* fp = std::fopen(out.c_str(), fresh ? "w" : "a");
    if (!fp) { std::fprintf(stderr, "# FATAL: cannot open %s\n", out.c_str()); return 2; }
    if (fresh)
        std::fprintf(fp,
            "session,regime,rule,ratio,fid,problem,algo,impl,threads,mb,seed,"
            "archive_cap,fe_budget,T_s,tau_s,k_predicted,runs_completed,"
            "best_fitness,total_wall_s,total_fes,order_idx\n");

    // ---- thermal warm-up (freeze 6.2) ----
    // Interleaving distributes thermal bias evenly across conditions but does
    // not remove it: a cold machine is systematically faster than one at
    // equilibrium, so measurements taken during the thermal transient are not
    // comparable to those taken after it. We drive a sustained load at the
    // session's full thread count until the machine reaches steady state, and
    // discard everything computed during it.
    if (warmup_sec > 0.0) {
        std::fprintf(stderr, "# warm-up: %.0f s at %d threads (results discarded)\n",
                     warmup_sec, cores);
        const Problem& wp = probs.at(PROBLEMS.back().fid);
        PSOParams wpp; wpp.swarm_size = SWARM;
        auto tw = Clock::now();
        long long reps = 0;
        while (seconds_since(tw) < warmup_sec) {
            pso_run(wp, wpp, PROBLEMS.back().femult * DIM, 999000ULL + reps,
                    true, cores);
            ++reps;
        }
        std::fprintf(stderr, "# warm-up complete: %lld runs, %.1f s\n",
                     reps, seconds_since(tw));
        std::fprintf(stderr, "# RECORD the steady-state CPU temperature now, "
                             "and the ambient temperature (freeze 6.2, 6.7).\n");
    } else {
        std::fprintf(stderr, "# WARNING: warm-up skipped (--warmup 0). Early "
                             "conditions will run on a cold machine and be "
                             "systematically faster. This is a deviation from "
                             "freeze 6.2 and must be logged.\n");
    }

    auto t_session = Clock::now();
    int idx = 0, skipped = 0;

    for (const auto& c : conds) {
        ++idx;
        std::string k2 = std::to_string(session) + "|" + c.regime + "|" + c.rule + "|" +
                         std::to_string(c.ratio) + "|" + std::to_string(c.fid) + "|" +
                         c.algo + "|" + c.impl + "|" + std::to_string(c.threads) + "|" +
                         [&]{ char b[16]; std::snprintf(b,sizeof b,"%.3f",c.mb); return std::string(b);}() +
                         "|" + std::to_string(c.seed);
        if (resume && done.count(k2)) { ++skipped; continue; }

        const ProblemSpec* ps = nullptr;
        for (const auto& p : PROBLEMS) if (p.fid == c.fid) ps = &p;
        const Problem& prob = probs[c.fid];
        const bool is_sapso = std::strcmp(c.algo, "SAPSO") == 0;
        const bool use_omp  = std::strcmp(c.impl, "omp") == 0;
        // ratio 0 is the sentinel for PSO under budget-only regimes, where no
        // surrogate exists and archive_cap is undefined.
        const int  cap      = (c.ratio == 0) ? 0 : ps->cap_for_ratio.at(c.ratio);
        const long long FE  = (long long)(ps->femult * DIM * c.mb);

        PSOParams pp;   pp.swarm_size = SWARM;
        SAPSOParams sp; sp.swarm_size = SWARM;
        sp.archive_cap = cap; sp.n_cand = SAPSO_NCAND; sp.refit_every = SAPSO_REFIT;

        // tau for THIS configuration (freeze 4.2: per problem x impl x threads)
        const int tau_ratio = is_sapso ? c.ratio : 0;   // 0 => plain PSO
        double tau = tau_lookup(TAU, c.fid, tau_ratio, c.threads, use_omp);

        double T = 0.0; int kpred = 0;
        if (std::strcmp(c.regime, "timefair") == 0) {
            double tau_pso   = tau_lookup(TAU, c.fid, 0,       c.threads, use_omp);
            double tau_sapso = tau_lookup(TAU, c.fid, c.ratio, c.threads, use_omp);
            switch (c.rule) {
                case 'B': T = tau_sapso; break;                       // primary
                case 'A': T = tau_pso;   break;
                case 'C': T = 0.5 * (tau_pso + tau_sapso); break;
            }
            kpred = (tau > 0.0) ? (int)std::floor(T / tau) : 0;
        }

        double best = 0.0, wall = 0.0; long long fes = 0; int completed = 1;

        if (std::strcmp(c.regime, "timefair") == 0) {
            // Algorithm 1 / 2: repeat runs until the wall-clock budget is spent.
            TimeFairResult r;
            if (is_sapso) {
                // run_time_fair drives pso_run; for SAPSO we inline the same loop
                r.best_of_restarts = std::numeric_limits<double>::infinity();
                r.completed_runs = 0; r.total_fes = 0;
                auto t0 = Clock::now(); int j = 0;
                while (true) {
                    double el = seconds_since(t0);
                    if (el >= T) break;
                    if ((el + tau > T) && j > 0) break;   // admission policy (freeze 4.3)
                    RunResult rr = sapso_run(prob, sp, FE, 500000ULL + c.seed * 1000 + j,
                                             use_omp, c.threads);
                    r.best_of_restarts = std::min(r.best_of_restarts, rr.best_fitness);
                    r.total_fes += rr.fes; ++r.completed_runs; ++j;
                }
                r.total_wall = seconds_since(t0);
            } else {
                r = run_time_fair(prob, pp, FE, T, tau,
                                  500000ULL + c.seed * 1000, use_omp, c.threads, true);
            }
            best = r.best_of_restarts; wall = r.total_wall;
            fes = r.total_fes; completed = r.completed_runs;
        } else {
            // FE-matched and multi-budget FE: a single run at the FE budget.
            auto t0 = Clock::now();
            RunResult rr = is_sapso
                ? sapso_run(prob, sp, FE, 500000ULL + c.seed * 1000, use_omp, c.threads)
                : pso_run  (prob, pp, FE, 500000ULL + c.seed * 1000, use_omp, c.threads);
            wall = seconds_since(t0);
            best = rr.best_fitness; fes = rr.fes; completed = 1;
        }

        std::fprintf(fp,
            "%d,%s,%c,%d,%d,%s,%s,%s,%d,%.3f,%d,%d,%lld,%.9f,%.9f,%d,%d,%.17g,%.9f,%lld,%d\n",
            session, c.regime, c.rule, c.ratio, c.fid, prob.name.c_str(),
            c.algo, c.impl, c.threads, c.mb, c.seed,
            is_sapso ? cap : 0, FE, T, tau, kpred, completed,
            best, wall, fes, idx);
        std::fflush(fp);   // checkpoint after every condition

        if (idx % 25 == 0)
            std::fprintf(stderr, "# [%5d/%5zu] elapsed=%.1f min\n",
                         idx, conds.size(), seconds_since(t_session) / 60.0);
    }

    std::fclose(fp);
    std::fprintf(stderr, "# session %d complete: %d run, %d skipped, %.1f min\n",
                 session, idx - skipped, skipped, seconds_since(t_session) / 60.0);
    return 0;
}
