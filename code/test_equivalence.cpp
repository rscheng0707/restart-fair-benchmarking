// test_equivalence.cpp — §3.4 verification gate
//
// Requirement: the OpenMP implementation must produce results equivalent to the
// serial implementation under matched seeds. We test the strongest form —
// BIT-IDENTICAL output — at thread counts 1, 2, 4, 8.
//
// Why bit-identical and not "close": if results merely agree to some tolerance,
// then any ranking difference observed at higher thread counts could be
// floating-point drift rather than a genuine parallel effect. Exact equality
// removes that alternative explanation entirely.
//
// This test MUST pass before any pilot run. If it fails, the RQ3 comparison is
// confounded and the result is uninterpretable.

#include <cstdio>
#include <cstdint>
#include <cstring>
#include "benchmarks.hpp"
#include "pso.hpp"

using namespace tfb;

static bool bit_equal(double a, double b) {
    uint64_t ua, ub;
    std::memcpy(&ua, &a, 8);
    std::memcpy(&ub, &b, 8);
    return ua == ub;
}

int main() {
    PSOParams p;
    auto suite = standard_suite(10);
    const long long FE = 10000;
    const int threads[] = {1, 2, 4, 8};
    int failures = 0, checks = 0;

    std::printf("%-12s %6s %8s  %-24s %-24s %s\n",
                "problem", "seed", "threads", "serial", "openmp", "result");
    std::printf("--------------------------------------------------------"
                "----------------------------------\n");

    for (const auto& prob : suite) {
        for (uint64_t seed : {1ULL, 42ULL, 12345ULL}) {
            RunResult ser = pso_run(prob, p, FE, seed, false, 1);
            for (int t : threads) {
                RunResult omp = pso_run(prob, p, FE, seed, true, t);
                bool ok = bit_equal(ser.best_fitness, omp.best_fitness)
                          && ser.fes == omp.fes
                          && ser.iterations == omp.iterations;
                ++checks;
                if (!ok) ++failures;
                if (!ok || (seed == 1ULL && t == 8)) {
                    std::printf("%-12s %6llu %8d  %-24.17g %-24.17g %s\n",
                                prob.name.c_str(), (unsigned long long)seed, t,
                                ser.best_fitness, omp.best_fitness,
                                ok ? "OK" : "*** MISMATCH ***");
                }
            }
        }
    }

    std::printf("\n%d checks, %d failures\n", checks, failures);
    if (failures == 0)
        std::printf("PASS - serial and OpenMP paths are bit-identical.\n"
                    "RQ3 comparison is unconfounded by search-dynamics differences.\n");
    else
        std::printf("FAIL - DO NOT PROCEED TO PILOTS.\n"
                    "Parallel path diverges from serial; RQ3 would be uninterpretable.\n");
    return failures ? 1 : 0;
}
